#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# H4 — Genetic plausibility (cis-eQTL × AF GWAS coloc) for the four
#       pre-registered cis-eQTL × GWAS anchors (small frozen set).
#
# Sources (locked tissue priority for primary vs blood sensitivity):
#   AF GWAS  : GCST90624411 (Miyazawa 2023 multi-ancestry, harmonised tabix)
#   eQTL     : eQTL Catalogue dataset QTD000251
#              (GTEx, heart_atrial_appendage, n=372, gene-level)
#
# Anchors (locked, §6):
#   matrix  : PRELP, OMD
#   metabolic: IDH2, NDUFA8
#
# Pipeline per anchor:
#   1. Tabix-slice eQTL within ±1 Mb of the gene; filter to gene_id.
#   2. Test cis-eQTL gate: lead variant p_eQTL < 5e-8.
#   3. Tabix-slice AF GWAS in same window.
#   4. Inner-join variants on chrom/pos + ref/alt (handles flips).
#   5. Run coloc.abf with default priors; record PP4.
#   6. Verdict tier (locked): supported / inconclusive / not supported.
#
# Cached tabix slices live under Project_P6/cache/ to avoid re-downloads.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Rsamtools)
  library(data.table)
  library(coloc)
})

find_root <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("--file=", "", grep("--file=", ca, value = TRUE)[1])
  if (length(f) && !is.na(f) && nzchar(f)) {
    d <- normalizePath(dirname(f), mustWork = FALSE)
    if (basename(d) == "scripts") return(dirname(d))
    return(d)
  }
  cwd <- getwd(); if (basename(cwd) == "scripts") return(dirname(cwd)); cwd
}
ROOT  <- find_root()
CACHE <- file.path(ROOT, "Project_P6", "cache")
OUT_DIR <- file.path(ROOT, "results")
dir.create(CACHE,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GWAS_URL <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000/GCST90624411/harmonised/GCST90624411.h.tsv.gz"
EQTL_URL <- "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000251/QTD000251.all.tsv.gz"

EQTL_NAMES <- c("molecular_trait_id", "chromosome", "position", "ref", "alt",
                "variant", "ma_samples", "maf", "pvalue", "beta", "se",
                "type", "ac", "an", "r2", "molecular_trait_object_id",
                "gene_id", "median_tpm", "rsid", "nlog10p")
GWAS_NAMES <- c("chromosome", "base_pair_location", "effect_allele",
                "other_allele", "beta", "standard_error",
                "effect_allele_frequency", "p_value", "rsid", "rs_id",
                "n", "hm_coordinate_conversion", "hm_code", "variant_id")

ANCHORS <- list(
  list(symbol = "PRELP",  axis = "matrix",
       gene_id = "ENSG00000188783",  chrom = 1L,
       start = 203475828L, end = 203491352L),
  list(symbol = "OMD",    axis = "matrix",
       gene_id = "ENSG00000127083",  chrom = 9L,
       start = 92414245L,  end = 92424461L),
  list(symbol = "IDH2",   axis = "metabolic",
       gene_id = "ENSG00000182054",  chrom = 15L,
       start = 90083045L,  end = 90102504L),
  list(symbol = "NDUFA8", axis = "metabolic",
       gene_id = "ENSG00000119421",  chrom = 9L,
       start = 122144058L, end = 122159819L)
)

WINDOW_BP    <- 1e6     # ±1 Mb cis window (locked, §6)
EQTL_GATE_P  <- 5e-8    # cis-eQTL gate (locked, §6)
COLOC_VERDICT_HI <- 0.7
COLOC_VERDICT_LO <- 0.3

# ---------------------------------------------------------------------------
# Tabix slice with on-disk cache.  Each (source, chrom, start, end) cached as
# a .tsv.gz under Project_P6/cache/.
# ---------------------------------------------------------------------------
slice_tabix <- function(url, chrom, start, end, source_tag, retries = 3) {
  key <- sprintf("%s_chr%s_%s_%s.tsv.gz", source_tag, chrom, start, end)
  cache <- file.path(CACHE, key)
  if (file.exists(cache) && file.size(cache) > 0) {
    cat(sprintf("  cached: %s\n", basename(cache)))
    return(fread(cache, sep = "\t", header = FALSE, showProgress = FALSE))
  }
  rng <- GRanges(as.character(chrom), IRanges(start = start, end = end))
  tf  <- TabixFile(url, index = paste0(url, ".tbi"))
  cat(sprintf("  tabix slice: %s @ %s:%d-%d\n",
              source_tag, as.character(chrom), start, end))
  for (k in seq_len(retries)) {
    txt <- tryCatch(scanTabix(tf, param = rng)[[1]],
                    error = function(e) {
                      cat(sprintf("    [retry %d] %s\n", k, conditionMessage(e)))
                      NULL
                    })
    if (!is.null(txt)) break
    Sys.sleep(2 * k)
  }
  if (is.null(txt) || length(txt) == 0L) {
    cat("    [empty slice]\n")
    return(data.table())
  }
  # Write to cache + return
  con <- gzfile(cache, "wb")
  writeLines(txt, con); close(con)
  fread(cache, sep = "\t", header = FALSE, showProgress = FALSE)
}

# ---------------------------------------------------------------------------
# Per-anchor analysis
# ---------------------------------------------------------------------------
process_anchor <- function(a) {
  cat(sprintf("\n=== %s (%s axis) chr%s:%d-%d ===\n",
              a$symbol, a$axis, a$chrom, a$start, a$end))
  win_lo <- max(1L, a$start - WINDOW_BP)
  win_hi <- a$end + WINDOW_BP

  # ---- eQTL slice ----
  e <- slice_tabix(EQTL_URL, a$chrom, win_lo, win_hi,
                   sprintf("eqtl_%s", a$symbol))
  if (nrow(e) == 0L) {
    return(data.table(
      anchor = a$symbol, axis = a$axis,
      n_eqtl = 0L, lead_eqtl_p = NA_real_, gate_passed = FALSE,
      n_gwas = NA_integer_, n_overlap = NA_integer_,
      gwas_min_p = NA_real_,
      coloc_PP0 = NA_real_, coloc_PP1 = NA_real_, coloc_PP2 = NA_real_,
      coloc_PP3 = NA_real_, coloc_PP4 = NA_real_,
      verdict = "not supported",
      reason = "eQTL slice empty"))
  }
  setnames(e, EQTL_NAMES[1:ncol(e)])
  # Filter to anchor gene (gene_id may carry version suffix or not)
  e <- e[grepl(paste0("^", a$gene_id, "(\\.|$)"), gene_id) |
         grepl(paste0("^", a$gene_id, "(\\.|$)"), molecular_trait_id)]
  cat(sprintf("  eQTL rows for %s: %d\n", a$symbol, nrow(e)))
  if (nrow(e) == 0L) {
    return(data.table(
      anchor = a$symbol, axis = a$axis,
      n_eqtl = 0L, lead_eqtl_p = NA_real_, gate_passed = FALSE,
      n_gwas = NA_integer_, n_overlap = NA_integer_,
      gwas_min_p = NA_real_,
      coloc_PP0 = NA_real_, coloc_PP1 = NA_real_, coloc_PP2 = NA_real_,
      coloc_PP3 = NA_real_, coloc_PP4 = NA_real_,
      verdict = "not supported",
      reason = "no eQTL rows for gene_id in window"))
  }
  e[, position := as.integer(position)]
  e[, pvalue := as.numeric(pvalue)]
  e[, beta := as.numeric(beta)]; e[, se := as.numeric(se)]
  e[, maf := as.numeric(maf)]
  e[, an := as.integer(an)]
  e <- e[!is.na(pvalue) & !is.na(beta) & !is.na(se) & se > 0 & maf > 0 & maf < 0.5]
  lead_eqtl_p <- min(e$pvalue, na.rm = TRUE)
  gate_passed <- isTRUE(lead_eqtl_p < EQTL_GATE_P)
  cat(sprintf("  eQTL lead p: %.2e   gate (<5e-8): %s\n",
              lead_eqtl_p, gate_passed))

  # ---- GWAS slice ----
  g <- slice_tabix(GWAS_URL, a$chrom, win_lo, win_hi,
                   sprintf("gwas_%s", a$symbol))
  setnames(g, GWAS_NAMES[1:ncol(g)])
  g[, base_pair_location := as.integer(base_pair_location)]
  g[, beta := as.numeric(beta)]
  g[, standard_error := as.numeric(standard_error)]
  g[, p_value := as.numeric(p_value)]
  g[, n := as.integer(n)]
  g <- g[!is.na(beta) & !is.na(standard_error) & standard_error > 0
         & !is.na(p_value)]
  gwas_min_p <- if (nrow(g)) min(g$p_value, na.rm = TRUE) else NA_real_
  cat(sprintf("  GWAS rows: %d   min p: %.2e\n", nrow(g), gwas_min_p))

  # ---- Variant matching ----
  e[, key1 := paste(chromosome, position, ref, alt, sep = "_")]
  e[, key2 := paste(chromosome, position, alt, ref, sep = "_")]
  g[, key  := paste(chromosome, base_pair_location, effect_allele,
                    other_allele, sep = "_")]

  joined_a <- merge(e, g, by.x = "key1", by.y = "key", suffixes = c(".e", ".g"))
  joined_a[, allele_flip := FALSE]
  joined_b <- merge(e, g, by.x = "key2", by.y = "key", suffixes = c(".e", ".g"))
  if (nrow(joined_b)) {
    joined_b[, allele_flip := TRUE]
    joined_b[, beta.g := -beta.g]   # flip GWAS beta to align alleles with eQTL
  }
  setnames(joined_a, "key1", "key"); setnames(joined_b, "key2", "key")
  joined <- rbindlist(list(joined_a, joined_b),
                      use.names = TRUE, fill = TRUE)
  joined <- unique(joined, by = "key")
  cat(sprintf("  overlap variants: %d (allele-flip applied: %d)\n",
              nrow(joined), sum(joined$allele_flip)))

  # We need MAF and N estimates for coloc
  joined <- joined[!is.na(maf) & maf > 0 & maf < 0.5]
  joined[, n_eqtl := an %/% 2L]    # an = 2 * sample_size for diploid

  # ---- Coloc ----
  if (nrow(joined) < 50L) {
    return(data.table(
      anchor = a$symbol, axis = a$axis,
      n_eqtl = nrow(e), lead_eqtl_p = lead_eqtl_p, gate_passed = gate_passed,
      n_gwas = nrow(g), n_overlap = nrow(joined),
      gwas_min_p = gwas_min_p,
      coloc_PP0 = NA_real_, coloc_PP1 = NA_real_, coloc_PP2 = NA_real_,
      coloc_PP3 = NA_real_, coloc_PP4 = NA_real_,
      verdict = if (gate_passed) "inconclusive" else "not supported",
      reason = sprintf("overlap n=%d < 50", nrow(joined))))
  }

  # AF prevalence assumption for coloc: ~5% in EUR (binary case-control)
  # For multi-ancestry, n is per-variant; coloc only needs N.
  ds_eqtl <- list(snp = joined$key, beta = joined$beta.e,
                  varbeta = joined$se^2, MAF = joined$maf,
                  N = max(joined$n_eqtl, na.rm = TRUE),
                  type = "quant", sdY = 1)
  ds_gwas <- list(snp = joined$key, beta = joined$beta.g,
                  varbeta = joined$standard_error^2, MAF = joined$maf,
                  N = as.integer(median(joined$n, na.rm = TRUE)),
                  s = 0.05,    # case fraction proxy (locked)
                  type = "cc")

  res <- tryCatch(
    coloc.abf(dataset1 = ds_eqtl, dataset2 = ds_gwas,
              p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
    error = function(e) list(error = conditionMessage(e))
  )
  if (!is.null(res$error)) {
    cat(sprintf("  coloc.abf error: %s\n", res$error))
    return(data.table(
      anchor = a$symbol, axis = a$axis,
      n_eqtl = nrow(e), lead_eqtl_p = lead_eqtl_p, gate_passed = gate_passed,
      n_gwas = nrow(g), n_overlap = nrow(joined),
      gwas_min_p = gwas_min_p,
      coloc_PP0 = NA_real_, coloc_PP1 = NA_real_, coloc_PP2 = NA_real_,
      coloc_PP3 = NA_real_, coloc_PP4 = NA_real_,
      verdict = if (gate_passed) "inconclusive" else "not supported",
      reason = paste("coloc error:", res$error)))
  }
  pp <- as.numeric(res$summary[c("PP.H0.abf","PP.H1.abf","PP.H2.abf",
                                  "PP.H3.abf","PP.H4.abf")])
  cat(sprintf("  coloc PP4 = %.3f   (PP3 = %.3f)\n", pp[5], pp[4]))

  verdict <- if (!gate_passed) {
    "not supported"
  } else if (pp[5] >= COLOC_VERDICT_HI) {
    "supported"
  } else if (pp[5] >= COLOC_VERDICT_LO) {
    "inconclusive"
  } else {
    "not supported"
  }
  reason <- sprintf("gate=%s PP4=%.2f", gate_passed, pp[5])

  data.table(
    anchor = a$symbol, axis = a$axis,
    n_eqtl = nrow(e), lead_eqtl_p = lead_eqtl_p, gate_passed = gate_passed,
    n_gwas = nrow(g), n_overlap = nrow(joined),
    gwas_min_p = gwas_min_p,
    coloc_PP0 = pp[1], coloc_PP1 = pp[2], coloc_PP2 = pp[3],
    coloc_PP3 = pp[4], coloc_PP4 = pp[5],
    verdict = verdict, reason = reason
  )
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
do_subset <- if (length(args)) args else NULL

results <- list()
for (a in ANCHORS) {
  if (!is.null(do_subset) && !(a$symbol %in% do_subset)) next
  results[[a$symbol]] <- process_anchor(a)
}

out <- rbindlist(results, use.names = TRUE, fill = TRUE)

OUT_TSV <- file.path(OUT_DIR, "figure6_mr_anchors.tsv")
existing <- if (file.exists(OUT_TSV)) fread(OUT_TSV) else data.table()
if (nrow(existing)) existing <- existing[!anchor %in% out$anchor]
out_combined <- rbindlist(list(existing, out), use.names = TRUE, fill = TRUE)
fwrite(out_combined, OUT_TSV, sep = "\t")

cat("\n=== Final summary ===\n")
print(out)
cat("\nH4 DONE.\n")
