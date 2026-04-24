#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# H4e — Per-anchor GTEx eGene status (Plan C eligibility primary gate).
#
# Spec lock: FIGURE6_WORKFLOW.md §10.3 (v3 follow-up, 2026-04-23).
#
# Purpose: replace the methodologically-misapplied lead-p < 5e-8 GWAS
#          significance gate with the canonical GTEx v8 eGene definition,
#          i.e. BH-FDR q < 0.05 on permutation-adjusted p_beta from the
#          eQTL Catalogue per-tissue .permuted.tsv.gz files.  Reported
#          per anchor for both atrial appendage (PRIMARY) and whole blood
#          (sensitivity informational only).
#
# Sources (locked, identical to h4 / h4b / h4c):
#   atrial_appendage : QTD000251.permuted.tsv.gz  (GTEx v8, n = 372)
#   whole_blood      : QTD000356.permuted.tsv.gz  (GTEx v8, n = 670)
#
# Output:
#   results/figure6_egene_status.tsv
#       columns: anchor, gene_id, tissue, n_genes_in_file,
#                p_beta, q_BH, is_eGene_q05
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
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

# Per-tissue permuted FDR files — single download, then cached locally.
PERM_URLS <- list(
  atrial_appendage = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000251/QTD000251.permuted.tsv.gz",
  whole_blood      = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000356/QTD000356.permuted.tsv.gz"
)

# All 11 anchor + calibrator gene_ids (locked, FIGURE6_WORKFLOW.md §10.2).
ANCHOR_GENES <- data.table(
  anchor = c("PRELP", "OMD", "IDH2", "NDUFA8",
             "NTRK2", "NTF3", "FRS2",
             "BMP10", "SMAD9", "SMAD6",
             "PITX2"),
  gene_id = c("ENSG00000188783", "ENSG00000127083",
              "ENSG00000182054", "ENSG00000119421",
              "ENSG00000148053", "ENSG00000185652", "ENSG00000166225",
              "ENSG00000163217", "ENSG00000120693", "ENSG00000137834",
              "ENSG00000164093")
)

cache_path <- function(tissue) {
  file.path(CACHE, sprintf("%s.permuted.tsv.gz",
                           c(atrial_appendage = "QTD000251",
                             whole_blood      = "QTD000356")[tissue]))
}

download_permuted <- function(tissue, url) {
  dest <- cache_path(tissue)
  if (file.exists(dest) && file.size(dest) > 1000L) {
    cat(sprintf("  cached: %s (%.1f MB)\n",
                basename(dest), file.size(dest) / 1024^2))
    return(dest)
  }
  cat(sprintf("  downloading %s\n  → %s\n", url, dest))
  utils::download.file(url, destfile = dest, mode = "wb",
                       quiet = FALSE, method = "libcurl")
  if (!file.exists(dest) || file.size(dest) < 1000L)
    stop("Download failed for ", tissue)
  dest
}

process_tissue <- function(tissue, url) {
  cat(sprintf("\n=== %s ===\n", tissue))
  f <- download_permuted(tissue, url)
  perm <- fread(f, sep = "\t", header = TRUE,
                showProgress = FALSE)
  cat(sprintf("  rows = %d   cols = %d\n", nrow(perm), ncol(perm)))
  cat(sprintf("  cols: %s\n",
              paste(names(perm), collapse = ", ")))

  # eQTL Catalogue permuted schema column for permutation-adjusted lead-p
  # is `p_beta` (canonical FastQTL --permute output).  Defensive lookup
  # in case of column-name drift across releases.
  p_col <- if ("p_beta" %in% names(perm)) "p_beta" else
           if ("pvalue_beta" %in% names(perm)) "pvalue_beta" else
           if ("p.beta" %in% names(perm)) "p.beta" else
           stop("Cannot find p_beta column in ", tissue)

  perm[, p_beta_num := as.numeric(get(p_col))]
  perm <- perm[!is.na(p_beta_num)]
  cat(sprintf("  with valid p_beta: %d\n", nrow(perm)))

  # BH-FDR across ALL gene-level permuted lead p-values
  perm[, q_BH := stats::p.adjust(p_beta_num, method = "BH")]
  n_eGene_all <- sum(perm$q_BH < 0.05)
  cat(sprintf("  all-tissue eGene count (q<0.05): %d / %d\n",
              n_eGene_all, nrow(perm)))

  # Filter to our 11 anchor gene_ids — eQTL Catalogue stores
  # molecular_trait_id with the bare ENSG (no version) for gene-level.
  gene_id_col <- if ("molecular_trait_id" %in% names(perm)) "molecular_trait_id"
                 else if ("gene_id" %in% names(perm)) "gene_id"
                 else stop("Cannot find gene_id column in ", tissue)

  perm[, gene_id_bare := sub("\\.[0-9]+$", "", get(gene_id_col))]
  out <- ANCHOR_GENES[
    perm,
    on = c(gene_id = "gene_id_bare"),
    nomatch = NULL
  ][, .(anchor, gene_id,
         p_beta       = p_beta_num,
         q_BH         = q_BH,
         is_eGene_q05 = q_BH < 0.05)]
  out[, tissue := tissue]
  out[, n_genes_in_file := nrow(perm)]

  # Make sure every anchor appears, even if absent from the permuted file.
  missing <- setdiff(ANCHOR_GENES$anchor, out$anchor)
  if (length(missing)) {
    miss <- ANCHOR_GENES[anchor %in% missing][, .(
      anchor, gene_id,
      p_beta = NA_real_, q_BH = NA_real_,
      is_eGene_q05 = FALSE,
      tissue = tissue,
      n_genes_in_file = nrow(perm)
    )]
    out <- rbind(out, miss, use.names = TRUE)
    cat(sprintf("  anchors absent from permuted file: %s\n",
                paste(missing, collapse = ", ")))
  }

  setcolorder(out, c("anchor", "gene_id", "tissue",
                      "n_genes_in_file",
                      "p_beta", "q_BH", "is_eGene_q05"))
  out
}

# Run both tissues
results <- lapply(names(PERM_URLS),
                  function(t) process_tissue(t, PERM_URLS[[t]]))
all_out <- rbindlist(results, use.names = TRUE)

# Sort: anchor canonical order, then tissue
canonical_order <- ANCHOR_GENES$anchor
all_out[, anchor := factor(anchor, levels = canonical_order)]
all_out[, tissue := factor(tissue,
                            levels = c("atrial_appendage", "whole_blood"))]
setorder(all_out, anchor, tissue)
all_out[, anchor := as.character(anchor)]
all_out[, tissue := as.character(tissue)]

OUT_TSV <- file.path(OUT_DIR, "figure6_egene_status.tsv")
fwrite(all_out, OUT_TSV, sep = "\t")

cat("\n=== eGene status summary ===\n")
print(all_out[, .(anchor, tissue, p_beta, q_BH, is_eGene_q05)])
cat(sprintf("\nWrote %s\n", OUT_TSV))
cat("H4e DONE.\n")
