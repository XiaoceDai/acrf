#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# H4c — Positive-control cis-eQTL × AF-GWAS coloc on canonical AF anchor
#       genes.  Used to demonstrate that the colocalisation pipeline
#       used for the four pre-registered framework anchors (PRELP, OMD,
#       IDH2, NDUFA8; h4_mr_coloc_anchors.R + h4b_mr_coloc_blood.R) is
#       analytically capable of detecting a positive cis-eQTL × AF-GWAS
#       signal where one is plausibly expected to exist.
#
# Positive controls (locked here; not part of the original framework
# anchor set):
#   PITX2  (ENSG00000164093, chr4:110,610,001-110,633,149) — canonical
#          4q25 AF GWAS locus; reported to have weak / variable mRNA
#          cis-eQTL in atrial tissue at common-variant resolution
#          because the locus is enhancer-mediated.
#   KCNN3  (ENSG00000143603, chr1:154,837,802-154,961,963) — 1q21.3 AF
#          GWAS locus; reported to coloc through cis-eQTL in atrial
#          tissue in independent published analyses.
#
# Tissue panels (same as h4 / h4b):
#   atrial_appendage : eQTL Catalogue QTD000251 (GTEx v8, n = 372)
#   whole_blood      : eQTL Catalogue QTD000356 (GTEx v8, n = 670)
#
# Pipeline mirrors h4 / h4b exactly:
#   ±1 Mb window, lead-p gate at 5e-8, coloc.abf with default priors
#   (p1 = p2 = 1e-4, p12 = 1e-5), allele-flip-aware variant matching,
#   and the same {supported / inconclusive / not supported} verdict
#   bins (PP4 ≥ 0.7 / 0.3 ≤ PP4 < 0.7 / PP4 < 0.3 OR gate-fail).
#
# Output:
#   results/figure6_mr_positive_controls.tsv
#       columns: anchor, tissue, axis, n_eqtl, lead_eqtl_p, gate_passed,
#                n_gwas, n_overlap, gwas_min_p, coloc_PP0..4, verdict, reason
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Rsamtools); library(data.table); library(coloc)
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
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GWAS_URL <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000/GCST90624411/harmonised/GCST90624411.h.tsv.gz"

EQTL_URLS <- list(
  atrial_appendage = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000251/QTD000251.all.tsv.gz",
  whole_blood      = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000356/QTD000356.all.tsv.gz"
)

EQTL_NAMES <- c("molecular_trait_id","chromosome","position","ref","alt",
                "variant","ma_samples","maf","pvalue","beta","se","type",
                "ac","an","r2","molecular_trait_object_id","gene_id",
                "median_tpm","rsid","nlog10p")
GWAS_NAMES <- c("chromosome","base_pair_location","effect_allele",
                "other_allele","beta","standard_error",
                "effect_allele_frequency","p_value","rsid","rs_id","n",
                "hm_coordinate_conversion","hm_code","variant_id")

# ---------------------------------------------------------------------------
# Pre-registered positive-control anchors (locked here in script header).
# ---------------------------------------------------------------------------
POSITIVE_CONTROLS <- list(
  list(symbol = "PITX2", axis = "positive_control",
       gene_id = "ENSG00000164093",
       chrom = 4L,  start = 110610001L, end = 110633149L,
       expected_gate = "uncertain_enhancer_mediated"),
  list(symbol = "KCNN3", axis = "positive_control",
       gene_id = "ENSG00000143603",
       chrom = 1L,  start = 154837802L, end = 154961963L,
       expected_gate = "expected_to_pass")
)

WINDOW_BP    <- 1e6
EQTL_GATE_P  <- 5e-8
COLOC_VERDICT_HI <- 0.7
COLOC_VERDICT_LO <- 0.3

# ---------------------------------------------------------------------------
# Tabix slice with on-disk cache (identical to h4 / h4b convention).
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
  con <- gzfile(cache, "wb")
  writeLines(txt, con); close(con)
  fread(cache, sep = "\t", header = FALSE, showProgress = FALSE)
}

# ---------------------------------------------------------------------------
# Per-anchor × per-tissue analysis.
# ---------------------------------------------------------------------------
process_anchor <- function(a, tissue_label, eqtl_url) {
  cat(sprintf("\n=== %s [%s axis] in %s   chr%s:%d-%d ===\n",
              a$symbol, a$axis, tissue_label,
              a$chrom, a$start, a$end))
  win_lo <- max(1L, a$start - WINDOW_BP)
  win_hi <- a$end + WINDOW_BP

  eqtl_tag <- if (tissue_label == "whole_blood")
    sprintf("eqtl_blood_%s", a$symbol)
  else
    sprintf("eqtl_%s", a$symbol)

  e <- slice_tabix(eqtl_url, a$chrom, win_lo, win_hi, eqtl_tag)
  empty_row <- data.table(
    anchor = a$symbol, tissue = tissue_label, axis = a$axis,
    n_eqtl = 0L, lead_eqtl_p = NA_real_, gate_passed = FALSE,
    n_gwas = NA_integer_, n_overlap = NA_integer_,
    gwas_min_p = NA_real_,
    coloc_PP0 = NA_real_, coloc_PP1 = NA_real_, coloc_PP2 = NA_real_,
    coloc_PP3 = NA_real_, coloc_PP4 = NA_real_,
    verdict = "not supported",
    reason = "eQTL slice empty"
  )
  if (nrow(e) == 0L) return(empty_row)

  setnames(e, EQTL_NAMES[1:ncol(e)])
  e <- e[grepl(paste0("^", a$gene_id, "(\\.|$)"), gene_id) |
         grepl(paste0("^", a$gene_id, "(\\.|$)"), molecular_trait_id)]
  cat(sprintf("  eQTL rows for %s: %d\n", a$symbol, nrow(e)))
  if (nrow(e) == 0L) {
    empty_row$reason <- "no eQTL rows for gene_id in window"
    return(empty_row)
  }
  e[, position := as.integer(position)]
  e[, pvalue := as.numeric(pvalue)]
  e[, beta := as.numeric(beta)]; e[, se := as.numeric(se)]
  e[, maf := as.numeric(maf)]
  e[, an := as.integer(an)]
  e <- e[!is.na(pvalue) & !is.na(beta) & !is.na(se) & se > 0
         & maf > 0 & maf < 0.5]
  lead_eqtl_p <- min(e$pvalue, na.rm = TRUE)
  gate_passed <- isTRUE(lead_eqtl_p < EQTL_GATE_P)
  cat(sprintf("  eQTL lead p: %.2e   gate (<5e-8): %s\n",
              lead_eqtl_p, gate_passed))

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

  e[, key1 := paste(chromosome, position, ref, alt, sep = "_")]
  e[, key2 := paste(chromosome, position, alt, ref, sep = "_")]
  g[, key  := paste(chromosome, base_pair_location, effect_allele,
                    other_allele, sep = "_")]

  joined_a <- merge(e, g, by.x = "key1", by.y = "key", suffixes = c(".e", ".g"))
  joined_a[, allele_flip := FALSE]
  joined_b <- merge(e, g, by.x = "key2", by.y = "key", suffixes = c(".e", ".g"))
  if (nrow(joined_b)) {
    joined_b[, allele_flip := TRUE]
    joined_b[, beta.g := -beta.g]
  }
  setnames(joined_a, "key1", "key"); setnames(joined_b, "key2", "key")
  joined <- rbindlist(list(joined_a, joined_b),
                      use.names = TRUE, fill = TRUE)
  joined <- unique(joined, by = "key")
  cat(sprintf("  overlap variants: %d (allele-flip applied: %d)\n",
              nrow(joined), sum(joined$allele_flip)))

  joined <- joined[!is.na(maf) & maf > 0 & maf < 0.5]
  joined[, n_eqtl := an %/% 2L]

  base_row <- data.table(
    anchor = a$symbol, tissue = tissue_label, axis = a$axis,
    n_eqtl = nrow(e), lead_eqtl_p = lead_eqtl_p, gate_passed = gate_passed,
    n_gwas = nrow(g), n_overlap = nrow(joined),
    gwas_min_p = gwas_min_p,
    coloc_PP0 = NA_real_, coloc_PP1 = NA_real_, coloc_PP2 = NA_real_,
    coloc_PP3 = NA_real_, coloc_PP4 = NA_real_,
    verdict = "not supported",
    reason = sprintf("overlap n=%d", nrow(joined))
  )
  if (nrow(joined) < 50L) {
    base_row$verdict <- if (gate_passed) "inconclusive" else "not supported"
    base_row$reason  <- sprintf("overlap n=%d < 50", nrow(joined))
    return(base_row)
  }

  ds_eqtl <- list(snp = joined$key, beta = joined$beta.e,
                  varbeta = joined$se^2, MAF = joined$maf,
                  N = max(joined$n_eqtl, na.rm = TRUE),
                  type = "quant", sdY = 1)
  ds_gwas <- list(snp = joined$key, beta = joined$beta.g,
                  varbeta = joined$standard_error^2, MAF = joined$maf,
                  N = as.integer(median(joined$n, na.rm = TRUE)),
                  s = 0.05, type = "cc")

  res <- tryCatch(
    coloc.abf(dataset1 = ds_eqtl, dataset2 = ds_gwas,
              p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
    error = function(e) list(error = conditionMessage(e))
  )
  if (!is.null(res$error)) {
    cat(sprintf("  coloc.abf error: %s\n", res$error))
    base_row$verdict <- if (gate_passed) "inconclusive" else "not supported"
    base_row$reason  <- paste("coloc error:", res$error)
    return(base_row)
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

  base_row$coloc_PP0 <- pp[1]; base_row$coloc_PP1 <- pp[2]
  base_row$coloc_PP2 <- pp[3]; base_row$coloc_PP3 <- pp[4]
  base_row$coloc_PP4 <- pp[5]; base_row$verdict <- verdict
  base_row$reason   <- reason
  base_row
}

# ---------------------------------------------------------------------------
# Run all positive controls × all tissues.
# ---------------------------------------------------------------------------
results <- list()
for (a in POSITIVE_CONTROLS) {
  for (tissue_label in names(EQTL_URLS)) {
    key <- sprintf("%s_%s", a$symbol, tissue_label)
    results[[key]] <- process_anchor(a, tissue_label, EQTL_URLS[[tissue_label]])
  }
}

out <- rbindlist(results, use.names = TRUE, fill = TRUE)
OUT_TSV <- file.path(OUT_DIR, "figure6_mr_positive_controls.tsv")
fwrite(out, OUT_TSV, sep = "\t")

cat("\n=== Final positive-control summary ===\n")
print(out[, .(anchor, tissue, lead_eqtl_p, gate_passed,
              coloc_PP4, verdict)])
cat(sprintf("\nWrote %s\n", OUT_TSV))
cat("H4c DONE.\n")
