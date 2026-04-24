#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# H4d — v3 LE-extension cis-eQTL × AF-GWAS coloc on six new anchors.
#
# Spec lock: FIGURE6_WORKFLOW.md §10.2 / §10.3 / §10.4.
#
# Anchors (locked, §10.2 v3 LE extension):
#   F15 (NTRK2/TrkB, RA-biased): NTRK2 (receptor), NTF3 (ligand),
#                                FRS2 (direct adaptor)
#   F11 (BMP, RA-biased)       : BMP10 (atrial-enriched ligand,
#                                AF GWAS hit), SMAD9 (R-SMAD),
#                                SMAD6 (i-SMAD)
#
# Sources (locked, §10.3, identical to v2 §7):
#   AF GWAS  : GCST90624411 (Miyazawa 2023 multi-ancestry, harmonised tabix)
#   eQTL prim: eQTL Catalogue dataset QTD000251
#              (GTEx, heart_atrial_appendage, n=372, gene-level)
#   eQTL sens: eQTL Catalogue dataset QTD000356
#              (GTEx, whole_blood,         n=670, gene-level)
#
# Pipeline mirrors h4 / h4b / h4c exactly:
#   ±1 Mb cis window, lead-p eligibility gate at 5e-8,
#   coloc.abf with default priors (p1 = p2 = 1e-4, p12 = 1e-5),
#   allele-flip-aware variant matching.
#
# v3 NEW: also reports
#   prefilter_passed = TRUE iff cis window contains ≥1 SNP with
#                       nominal p_eQTL < 1e-4 (regardless of gate).
#
# Output:
#   results/figure6_mr_anchors_v3.tsv
#       columns: anchor, family, role, source, tissue, axis,
#                n_eqtl, lead_eqtl_p, gate_passed, prefilter_passed,
#                n_gwas, n_overlap, gwas_min_p,
#                coloc_PP0..4, reason
#       (verdict / shape_tag / family_verdict / scenario are computed
#        downstream by build_v3_consolidated.py per §10.5–10.8.)
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
dir.create(CACHE,  recursive = TRUE, showWarnings = FALSE)
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
# v3 extension anchors (locked here in script header per §10.2).
# Coordinates: GRCh38 (Ensembl).
# ---------------------------------------------------------------------------
ANCHORS <- list(
  # F15 NTRK2 / TrkB family (RA-biased)
  list(symbol = "NTRK2",  family = "F15", role = "receptor",
       gene_id = "ENSG00000148053", chrom = 9L,
       start = 84668522L,  end = 85029699L),
  list(symbol = "NTF3",   family = "F15", role = "ligand",
       gene_id = "ENSG00000185652", chrom = 12L,
       start = 5404888L,   end = 5495026L),
  list(symbol = "FRS2",   family = "F15", role = "adaptor",
       gene_id = "ENSG00000166225", chrom = 12L,
       start = 69464533L,  end = 69513650L),
  # F11 BMP family (RA-biased)
  list(symbol = "BMP10",  family = "F11", role = "ligand",
       gene_id = "ENSG00000163217", chrom = 2L,
       start = 68869054L,  end = 68886842L),
  list(symbol = "SMAD9",  family = "F11", role = "R-SMAD",
       gene_id = "ENSG00000120693", chrom = 13L,
       start = 36879326L,  end = 36925087L),
  list(symbol = "SMAD6",  family = "F11", role = "i-SMAD",
       gene_id = "ENSG00000137834", chrom = 15L,
       start = 66703041L,  end = 66782082L)
)

WINDOW_BP        <- 1e6   # ±1 Mb cis window (locked, §10.3)
EQTL_GATE_P      <- 5e-8  # eligibility gate (locked, §10.3, inherited from v2)
EQTL_PREFILTER_P <- 1e-4  # coloc input prefilter (locked, §10.3 v3 NEW)

# ---------------------------------------------------------------------------
# Tabix slice with on-disk cache (identical to h4 / h4b / h4c).
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
empty_row <- function(a, tissue_label, reason) {
  data.table(
    anchor = a$symbol, family = a$family, role = a$role,
    source = "v3_LE_extension", tissue = tissue_label,
    axis = paste0(a$family, "_", a$role),
    n_eqtl = 0L, lead_eqtl_p = NA_real_,
    gate_passed = FALSE, prefilter_passed = FALSE,
    n_gwas = NA_integer_, n_overlap = NA_integer_,
    gwas_min_p = NA_real_,
    coloc_PP0 = NA_real_, coloc_PP1 = NA_real_, coloc_PP2 = NA_real_,
    coloc_PP3 = NA_real_, coloc_PP4 = NA_real_,
    reason = reason
  )
}

process_anchor <- function(a, tissue_label, eqtl_url) {
  cat(sprintf("\n=== %s [%s · %s] in %s   chr%s:%d-%d ===\n",
              a$symbol, a$family, a$role, tissue_label,
              a$chrom, a$start, a$end))
  win_lo <- max(1L, a$start - WINDOW_BP)
  win_hi <- a$end + WINDOW_BP

  eqtl_tag <- if (tissue_label == "whole_blood")
    sprintf("eqtl_blood_%s", a$symbol)
  else
    sprintf("eqtl_%s", a$symbol)

  e <- slice_tabix(eqtl_url, a$chrom, win_lo, win_hi, eqtl_tag)
  if (nrow(e) == 0L) return(empty_row(a, tissue_label, "eQTL slice empty"))

  setnames(e, EQTL_NAMES[1:ncol(e)])
  e <- e[grepl(paste0("^", a$gene_id, "(\\.|$)"), gene_id) |
         grepl(paste0("^", a$gene_id, "(\\.|$)"), molecular_trait_id)]
  cat(sprintf("  eQTL rows for %s: %d\n", a$symbol, nrow(e)))
  if (nrow(e) == 0L)
    return(empty_row(a, tissue_label, "no eQTL rows for gene_id in window"))

  e[, position := as.integer(position)]
  e[, pvalue := as.numeric(pvalue)]
  e[, beta := as.numeric(beta)]; e[, se := as.numeric(se)]
  e[, maf := as.numeric(maf)]
  e[, an := as.integer(an)]
  e <- e[!is.na(pvalue) & !is.na(beta) & !is.na(se) & se > 0
         & maf > 0 & maf < 0.5]
  lead_eqtl_p <- min(e$pvalue, na.rm = TRUE)
  gate_passed      <- isTRUE(lead_eqtl_p < EQTL_GATE_P)
  prefilter_passed <- isTRUE(lead_eqtl_p < EQTL_PREFILTER_P)
  cat(sprintf("  eQTL lead p: %.2e   gate (<5e-8): %s   prefilter (<1e-4): %s\n",
              lead_eqtl_p, gate_passed, prefilter_passed))

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

  # ---- Variant matching (allele-flip aware) ----
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
    anchor = a$symbol, family = a$family, role = a$role,
    source = "v3_LE_extension", tissue = tissue_label,
    axis = paste0(a$family, "_", a$role),
    n_eqtl = nrow(e), lead_eqtl_p = lead_eqtl_p,
    gate_passed = gate_passed, prefilter_passed = prefilter_passed,
    n_gwas = nrow(g), n_overlap = nrow(joined),
    gwas_min_p = gwas_min_p,
    coloc_PP0 = NA_real_, coloc_PP1 = NA_real_, coloc_PP2 = NA_real_,
    coloc_PP3 = NA_real_, coloc_PP4 = NA_real_,
    reason = sprintf("overlap n=%d", nrow(joined))
  )
  if (nrow(joined) < 50L) {
    base_row$reason <- sprintf("overlap n=%d < 50", nrow(joined))
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
    base_row$reason <- paste("coloc error:", res$error)
    return(base_row)
  }
  pp <- as.numeric(res$summary[c("PP.H0.abf","PP.H1.abf","PP.H2.abf",
                                  "PP.H3.abf","PP.H4.abf")])
  cat(sprintf("  coloc PP4 = %.3f   (PP3 = %.3f)\n", pp[5], pp[4]))

  base_row$coloc_PP0 <- pp[1]; base_row$coloc_PP1 <- pp[2]
  base_row$coloc_PP2 <- pp[3]; base_row$coloc_PP3 <- pp[4]
  base_row$coloc_PP4 <- pp[5]
  base_row$reason   <- sprintf("gate=%s pre=%s PP4=%.2f",
                               gate_passed, prefilter_passed, pp[5])
  base_row
}

# ---------------------------------------------------------------------------
# Run all v3 anchors × all tissues
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
do_subset <- if (length(args)) args else NULL

results <- list()
for (a in ANCHORS) {
  if (!is.null(do_subset) && !(a$symbol %in% do_subset)) next
  for (tissue_label in names(EQTL_URLS)) {
    key <- sprintf("%s_%s", a$symbol, tissue_label)
    results[[key]] <- process_anchor(a, tissue_label,
                                     EQTL_URLS[[tissue_label]])
  }
}

out <- rbindlist(results, use.names = TRUE, fill = TRUE)
OUT_TSV <- file.path(OUT_DIR, "figure6_mr_anchors_v3.tsv")
fwrite(out, OUT_TSV, sep = "\t")

cat("\n=== Final v3-extension summary ===\n")
print(out[, .(anchor, family, role, tissue,
              lead_eqtl_p, gate_passed, prefilter_passed,
              coloc_PP4, coloc_PP3)])
cat(sprintf("\nWrote %s\n", OUT_TSV))
cat("H4d DONE.\n")
