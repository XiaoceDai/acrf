#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# H4d — Unified cis-eQTL × AF-GWAS coloc for all 12 pre-specified anchors.
#
# Mirrors h4_mr_coloc_anchors.R / h4c_mr_coloc_positive_controls.R pipeline.
# Writes a dated rerun TSV for diff against frozen pipeline_wide (does NOT
# overwrite frozen artifacts).
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
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
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

# Locked 12-gene anchor list (±1 Mb windows from build_figure6_mr_v3_manuscript_tables.py)
ANCHORS <- list(
  list(symbol = "PRELP",  axis = "matrix", gene_id = "ENSG00000188783",
       chrom = 1L,  start = 203475828L, end = 203491352L),
  list(symbol = "OMD",    axis = "matrix", gene_id = "ENSG00000127083",
       chrom = 9L,  start = 92414245L,  end = 92424461L),
  list(symbol = "IDH2",   axis = "metabolic", gene_id = "ENSG00000182054",
       chrom = 15L, start = 90083045L, end = 90102504L),
  list(symbol = "NDUFA8", axis = "metabolic", gene_id = "ENSG00000119421",
       chrom = 9L,  start = 122144058L, end = 122159819L),
  list(symbol = "NTRK2",  axis = "neurotrophin", gene_id = "ENSG00000148053",
       chrom = 9L,  start = 84668522L, end = 85029699L),
  list(symbol = "NTF3",   axis = "neurotrophin", gene_id = "ENSG00000185652",
       chrom = 12L, start = 5404888L, end = 5495026L),
  list(symbol = "FRS2",   axis = "neurotrophin", gene_id = "ENSG00000166225",
       chrom = 12L, start = 69464533L, end = 69513650L),
  list(symbol = "BMP10",  axis = "BMP-SMAD", gene_id = "ENSG00000163217",
       chrom = 2L,  start = 68869054L, end = 68886842L),
  list(symbol = "SMAD9",  axis = "BMP-SMAD", gene_id = "ENSG00000120693",
       chrom = 13L, start = 36879326L, end = 36925087L),
  list(symbol = "SMAD6",  axis = "BMP-SMAD", gene_id = "ENSG00000137834",
       chrom = 15L, start = 66703041L, end = 66782082L),
  list(symbol = "PITX2",  axis = "calibration", gene_id = "ENSG00000164093",
       chrom = 4L,  start = 110610001L, end = 110633149L),
  list(symbol = "KCNN3",  axis = "calibration", gene_id = "ENSG00000143603",
       chrom = 1L,  start = 154837802L, end = 154961963L)
)

WINDOW_BP    <- 1e6
EQTL_GATE_P  <- 5e-8
COLOC_VERDICT_HI <- 0.7
COLOC_VERDICT_LO <- 0.3

slice_tabix <- function(url, chrom, start, end, source_tag, retries = 3) {
  key <- sprintf("%s_chr%s_%s_%s.tsv.gz", source_tag, chrom, start, end)
  cache <- file.path(CACHE, key)
  if (file.exists(cache) && file.size(cache) > 0) {
    cat(sprintf("  cached: %s\n", basename(cache)))
    return(fread(cache, sep = "\t", header = FALSE, showProgress = FALSE))
  }
  rng <- GRanges(as.character(chrom), IRanges(start = start, end = end))
  tf  <- TabixFile(url, index = paste0(url, ".tbi"))
  cat(sprintf("  tabix slice: %s @ %s:%d-%d\n", source_tag, chrom, start, end))
  for (k in seq_len(retries)) {
    txt <- tryCatch(scanTabix(tf, param = rng)[[1]],
                    error = function(e) { cat(sprintf("    [retry %d] %s\n", k, conditionMessage(e))); NULL })
    if (!is.null(txt)) break
    Sys.sleep(2 * k)
  }
  if (is.null(txt) || length(txt) == 0L) return(data.table())
  con <- gzfile(cache, "wb"); writeLines(txt, con); close(con)
  fread(cache, sep = "\t", header = FALSE, showProgress = FALSE)
}

process_anchor <- function(a) {
  cat(sprintf("\n=== %s (%s) chr%s:%d-%d ===\n", a$symbol, a$axis, a$chrom, a$start, a$end))
  win_lo <- max(1L, a$start - WINDOW_BP)
  win_hi <- a$end + WINDOW_BP

  e <- slice_tabix(EQTL_URL, a$chrom, win_lo, win_hi, sprintf("eqtl_%s", a$symbol))
  empty <- data.table(anchor = a$symbol, axis = a$axis, tissue = "atrial_appendage",
                      n_eqtl = 0L, lead_eqtl_p = NA_real_, gate_passed = FALSE,
                      n_gwas = NA_integer_, n_overlap = NA_integer_,
                      gwas_min_p = NA_real_,
                      coloc_PP0 = NA_real_, coloc_PP1 = NA_real_, coloc_PP2 = NA_real_,
                      coloc_PP3 = NA_real_, coloc_PP4 = NA_real_,
                      verdict = "not supported", reason = "eQTL slice empty")
  if (nrow(e) == 0L) return(empty)

  setnames(e, EQTL_NAMES[1:ncol(e)])
  e <- e[grepl(paste0("^", a$gene_id, "(\\.|$)"), gene_id) |
         grepl(paste0("^", a$gene_id, "(\\.|$)"), molecular_trait_id)]
  if (nrow(e) == 0L) { empty$reason <- "no eQTL rows for gene_id"; return(empty) }

  e[, position := as.integer(position)]
  e[, pvalue := as.numeric(pvalue)]; e[, beta := as.numeric(beta)]; e[, se := as.numeric(se)]
  e[, maf := as.numeric(maf)]; e[, an := as.integer(an)]
  e <- e[!is.na(pvalue) & !is.na(beta) & !is.na(se) & se > 0 & maf > 0 & maf < 0.5]
  lead_eqtl_p <- min(e$pvalue, na.rm = TRUE)
  gate_passed <- isTRUE(lead_eqtl_p < EQTL_GATE_P)

  g <- slice_tabix(GWAS_URL, a$chrom, win_lo, win_hi, sprintf("gwas_%s", a$symbol))
  setnames(g, GWAS_NAMES[1:ncol(g)])
  g[, base_pair_location := as.integer(base_pair_location)]
  g[, beta := as.numeric(beta)]; g[, standard_error := as.numeric(standard_error)]
  g[, p_value := as.numeric(p_value)]; g[, n := as.integer(n)]
  g <- g[!is.na(beta) & !is.na(standard_error) & standard_error > 0 & !is.na(p_value)]
  gwas_min_p <- if (nrow(g)) min(g$p_value, na.rm = TRUE) else NA_real_

  e[, key1 := paste(chromosome, position, ref, alt, sep = "_")]
  e[, key2 := paste(chromosome, position, alt, ref, sep = "_")]
  g[, key  := paste(chromosome, base_pair_location, effect_allele, other_allele, sep = "_")]
  joined_a <- merge(e, g, by.x = "key1", by.y = "key", suffixes = c(".e", ".g"))
  joined_a[, allele_flip := FALSE]
  joined_b <- merge(e, g, by.x = "key2", by.y = "key", suffixes = c(".e", ".g"))
  if (nrow(joined_b)) { joined_b[, allele_flip := TRUE]; joined_b[, beta.g := -beta.g] }
  setnames(joined_a, "key1", "key"); setnames(joined_b, "key2", "key")
  joined <- unique(rbindlist(list(joined_a, joined_b), use.names = TRUE, fill = TRUE), by = "key")
  joined <- joined[!is.na(maf) & maf > 0 & maf < 0.5]
  joined[, n_eqtl := an %/% 2L]

  if (nrow(joined) < 50L) {
    return(data.table(anchor = a$symbol, axis = a$axis, tissue = "atrial_appendage",
                        n_eqtl = nrow(e), lead_eqtl_p = lead_eqtl_p, gate_passed = gate_passed,
                        n_gwas = nrow(g), n_overlap = nrow(joined), gwas_min_p = gwas_min_p,
                        coloc_PP0 = NA_real_, coloc_PP1 = NA_real_, coloc_PP2 = NA_real_,
                        coloc_PP3 = NA_real_, coloc_PP4 = NA_real_,
                        verdict = if (gate_passed) "inconclusive" else "not supported",
                        reason = sprintf("overlap n=%d < 50", nrow(joined))))
  }

  ds_eqtl <- list(snp = joined$key, beta = joined$beta.e, varbeta = joined$se^2,
                  MAF = joined$maf, N = max(joined$n_eqtl, na.rm = TRUE),
                  type = "quant", sdY = 1)
  ds_gwas <- list(snp = joined$key, beta = joined$beta.g,
                  varbeta = joined$standard_error^2, MAF = joined$maf,
                  N = as.integer(median(joined$n, na.rm = TRUE)),
                  s = 0.05, type = "cc")
  res <- tryCatch(coloc.abf(dataset1 = ds_eqtl, dataset2 = ds_gwas,
                              p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
                  error = function(e) list(error = conditionMessage(e)))
  if (!is.null(res$error)) {
    return(data.table(anchor = a$symbol, axis = a$axis, tissue = "atrial_appendage",
                      n_eqtl = nrow(e), lead_eqtl_p = lead_eqtl_p, gate_passed = gate_passed,
                      n_gwas = nrow(g), n_overlap = nrow(joined), gwas_min_p = gwas_min_p,
                      coloc_PP0 = NA_real_, coloc_PP1 = NA_real_, coloc_PP2 = NA_real_,
                      coloc_PP3 = NA_real_, coloc_PP4 = NA_real_,
                      verdict = "not supported", reason = paste("coloc error:", res$error)))
  }
  pp <- as.numeric(res$summary[c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")])
  verdict <- if (!gate_passed) "not supported"
             else if (pp[5] >= COLOC_VERDICT_HI) "supported"
             else if (pp[5] >= COLOC_VERDICT_LO) "inconclusive"
             else "not supported"
  data.table(anchor = a$symbol, axis = a$axis, tissue = "atrial_appendage",
             n_eqtl = nrow(e), lead_eqtl_p = lead_eqtl_p, gate_passed = gate_passed,
             n_gwas = nrow(g), n_overlap = nrow(joined), gwas_min_p = gwas_min_p,
             coloc_PP0 = pp[1], coloc_PP1 = pp[2], coloc_PP2 = pp[3],
             coloc_PP3 = pp[4], coloc_PP4 = pp[5],
             verdict = verdict, reason = sprintf("gate=%s PP4=%.3f", gate_passed, pp[5]))
}

results <- lapply(ANCHORS, process_anchor)
out <- rbindlist(results, use.names = TRUE, fill = TRUE)
stamp <- format(Sys.Date(), "%Y%m%d")
OUT_TSV <- file.path(OUT_DIR, sprintf("figure6_mr_anchors_v3_rerun_%s.tsv", stamp))
fwrite(out, OUT_TSV, sep = "\t")
cat(sprintf("\nWrote %s (%d rows)\n", OUT_TSV, nrow(out)))
cat("H4d DONE — frozen pipeline_wide NOT modified.\n")
