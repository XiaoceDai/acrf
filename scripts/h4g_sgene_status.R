#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# H4g — Per-anchor GTEx atrial sGene status (Plan C eligibility for sQTL).
#
# Spec lock: FIGURE6_WORKFLOW.md §10.14 (sQTL RNA-layer extension).
#
# Mirrors h4e_egene_status.R but for splice QTL.  An anchor is called an
# "sGene" (q < 0.05) if ANY leafcutter cluster mapped to that gene's
# Ensembl ID has BH-FDR-corrected p_beta < 0.05 in the per-cluster
# permutation file.
#
# Source: eQTL Catalogue QTD000255.permuted.tsv.gz
#         (GTEx v8 atrial appendage, leafcutter, n=372)
#
# Output:
#   results/figure6_sgene_status.tsv
#       anchor, gene_id, n_clusters_in_file, p_beta_min, q_BH_min, is_sGene_q05
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
})

# Allow long downloads (eQTL Catalogue mirror is slow at times)
options(timeout = 1800)

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

PERM_URL <- "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000255/QTD000255.permuted.tsv.gz"
PERM_CACHE <- file.path(CACHE, "QTD000255.permuted.tsv.gz")

# h4g depends on h4f having cached the per-anchor sQTL slices, because the
# leafcutter permuted file does not store a gene_id column — we therefore
# build the cluster_id -> gene_id mapping from the per-anchor cc.tsv.gz
# slices that h4f writes to Project_P6/cache/.
SLICE_GLOB <- file.path(CACHE, "sqtl_*.tsv.gz")
SLICE_NAMES <- c("molecular_trait_id","chromosome","position","ref","alt",
                  "variant","ma_samples","maf","pvalue","beta","se","type",
                  "ac","an","r2","molecular_trait_object_id","gene_id",
                  "median_tpm","rsid")

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

# ---- Download once ----
if (!file.exists(PERM_CACHE) || file.size(PERM_CACHE) < 1000L) {
  cat(sprintf("downloading %s\n  → %s\n", PERM_URL, PERM_CACHE))
  utils::download.file(PERM_URL, destfile = PERM_CACHE, mode = "wb",
                       quiet = FALSE, method = "libcurl")
} else {
  cat(sprintf("cached: %s (%.1f MB)\n",
              basename(PERM_CACHE), file.size(PERM_CACHE) / 1024^2))
}

perm <- fread(PERM_CACHE, sep = "\t", header = TRUE, showProgress = FALSE)
cat(sprintf("rows = %d   cols = %d\n", nrow(perm), ncol(perm)))
cat(sprintf("cols: %s\n", paste(names(perm), collapse = ", ")))

p_col <- if ("p_beta" %in% names(perm)) "p_beta" else
         if ("pvalue_beta" %in% names(perm)) "pvalue_beta" else
         stop("Cannot find p_beta column")
perm[, p_beta_num := as.numeric(get(p_col))]
perm <- perm[!is.na(p_beta_num)]
cat(sprintf("with valid p_beta: %d\n", nrow(perm)))

# BH-FDR across ALL clusters (cluster-level FDR, not gene-level)
perm[, q_BH := stats::p.adjust(p_beta_num, method = "BH")]

# Build cluster_id -> gene_id mapping from per-anchor sQTL slice cache.
# The eQTL Catalogue leafcutter permuted file does NOT contain gene_id;
# the mapping lives in the .cc.tsv.gz per-variant rows that h4f cached.
slice_files <- Sys.glob(SLICE_GLOB)
if (length(slice_files) == 0L)
  stop("No per-anchor sqtl_*.tsv.gz slices found under ", CACHE,
       "; run scripts/h4f_sqtl_coloc.R first.")

cat(sprintf("\nBuilding cluster_id -> gene_id map from %d cached slice(s):\n",
            length(slice_files)))
cluster_gene_map <- rbindlist(lapply(slice_files, function(f) {
  cat(sprintf("  %s\n", basename(f)))
  d <- fread(f, sep = "\t", header = FALSE, showProgress = FALSE)
  if (ncol(d) < 17L) return(data.table())
  setnames(d, SLICE_NAMES[1:ncol(d)])
  unique(d[, .(molecular_trait_object_id, gene_id)])
}), use.names = TRUE, fill = TRUE)
cluster_gene_map <- unique(cluster_gene_map)
cat(sprintf("  unique cluster_id <-> gene_id pairs: %d\n",
            nrow(cluster_gene_map)))
cluster_gene_map[, gene_id_bare := sub("\\.[0-9]+$", "", gene_id)]

# Join cluster->gene mapping into the permuted table.
# molecular_trait_object_id is the cluster id in the permuted file.
perm <- merge(perm,
              cluster_gene_map[, .(molecular_trait_object_id, gene_id_bare)],
              by = "molecular_trait_object_id",
              all.x = TRUE)
n_with_gene <- sum(!is.na(perm$gene_id_bare))
cat(sprintf("permuted clusters with gene_id mapping: %d / %d\n",
            n_with_gene, nrow(perm)))

# Per anchor: is_sGene = ANY cluster of that gene has q_BH < 0.05
out <- ANCHOR_GENES[, {
  sub <- perm[gene_id_bare == gene_id]
  list(
    n_clusters_in_file = nrow(sub),
    p_beta_min   = if (nrow(sub)) min(sub$p_beta_num, na.rm = TRUE) else NA_real_,
    q_BH_min     = if (nrow(sub)) min(sub$q_BH,         na.rm = TRUE) else NA_real_,
    is_sGene_q05 = if (nrow(sub)) any(sub$q_BH < 0.05, na.rm = TRUE) else FALSE
  )
}, by = .(anchor, gene_id)]

n_sGene_all <- sum(perm$q_BH < 0.05)
cat(sprintf("\nAll-tissue cluster-level sQTL count (q<0.05): %d / %d\n",
            n_sGene_all, nrow(perm)))

OUT_TSV <- file.path(OUT_DIR, "figure6_sgene_status.tsv")
fwrite(out, OUT_TSV, sep = "\t")

cat("\n=== Per-anchor sGene status ===\n")
print(out)
cat(sprintf("\nWrote %s\n", OUT_TSV))
cat("H4g DONE.\n")
