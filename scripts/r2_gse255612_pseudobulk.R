#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# R2 - Donor x compartment pseudobulk for GSE255612.
#
# Input matrix is log-normalized per-cell expression (author-processed in
# GEO, not raw counts).  For Figure 3 we only feed these pseudobulks into
# singscore (rank-based), so the sample-level aggregator must be
# rank-preserving within a sample.  Choice (pseudobulk normalisation note in header):
#
#   pseudobulk_{donor,ct}(gene) = mean over cells in (donor, ct) of the
#                                 log-normalized expression
#
# Duplicate gene symbols (multiple ENSGs per HGNC symbol) are collapsed
# AFTER aggregation by taking the mean of their log-normalized values
# (also rank-preserving for singscore).
#
# Gate (spec §3.6): only donor x ct with n_cells >= 30 are emitted; gate
# violations are logged to qc_pseudobulk.tsv.  r1 already confirmed no
# violation in GSE255612.
#
# Output layout:
#   Project_P3/GSE255612/pseudobulk/
#   ├── <compartment4>/                    (atrialCM, fibroblast, endothelial_or_vascular, immune)
#   │   ├── expr.tsv                       gene x donor mean log-norm matrix
#   │   └── samples.tsv                    donor x (rhythm, sex, n_cells)
#   ├── supp_pericyte/                     same schema, for Supp
#   ├── supp_vascular_smooth_muscle/       same schema, for Supp
#   ├── summary.tsv                        per-compartment donor counts
#   └── qc_pseudobulk.tsv                  any (donor, ct) dropped by gate
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
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
  cwd <- getwd()
  if (basename(cwd) == "scripts") return(dirname(cwd))
  cwd
}
ROOT    <- find_root()
IN_DIR  <- file.path(ROOT, "Project_P3", "GSE255612")
OUT_DIR <- file.path(IN_DIR, "pseudobulk")
QC_DIR  <- IN_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("ROOT    = %s\n", ROOT))
cat(sprintf("IN_DIR  = %s\n", IN_DIR))
cat(sprintf("OUT_DIR = %s\n", OUT_DIR))

MIN_CELLS <- 30L      # spec §3.6 gate


# --------------------------------------------------------------------
# 1. Load barcodes, genes, obs
# --------------------------------------------------------------------

cat("reading barcodes / genes / obs ...\n")
barcodes <- fread(file.path(IN_DIR, "GSE255612_AF_snRNA_Processed_Expression_Matrix_barcodes_V1.tsv.gz"),
                  header = FALSE, col.names = "barcode")$barcode
genes <- fread(file.path(IN_DIR, "GSE255612_AF_snRNA_Processed_Expression_Matrix_genes_final.tsv.gz"),
               header = FALSE, col.names = c("ensg", "symbol"))
obs <- fread(file.path(IN_DIR, "derived", "gse255612_obs_4compartment.tsv"))

cat(sprintf("  barcodes: %d\n", length(barcodes)))
cat(sprintf("  genes:    %d  (symbols unique: %d)\n",
            nrow(genes), length(unique(genes$symbol))))
cat(sprintf("  obs kept: %d  (main=%d, supp=%d)\n",
            nrow(obs), sum(obs$keep_for_main), sum(obs$flagged_for_supp)))


# --------------------------------------------------------------------
# 2. Load the sparse matrix (genes x cells)
# --------------------------------------------------------------------
# readMM handles gzfile connections directly in recent Matrix versions.
cat("reading mtx (this takes ~1 min) ...\n")
t0 <- Sys.time()
M <- readMM(gzfile(file.path(IN_DIR, "GSE255612_AF_snRNA_Matrix_V1.mtx.gz")))
t1 <- Sys.time()
cat(sprintf("  mtx: %dx%d  %.0f Mnnz  (%s)\n",
            nrow(M), ncol(M), nnzero(M) / 1e6,
            format(difftime(t1, t0, units = "secs"))))
stopifnot(nrow(M) == nrow(genes), ncol(M) == length(barcodes))
# Convert to column-sparse for fast column extraction.
M <- as(M, "CsparseMatrix")
rownames(M) <- genes$ensg
colnames(M) <- barcodes


# --------------------------------------------------------------------
# 3. Restrict to obs-kept cells, attach labels
# --------------------------------------------------------------------

keep_idx <- match(obs$cell_id, colnames(M))
stopifnot(!anyNA(keep_idx))
M <- M[, keep_idx, drop = FALSE]
cat(sprintf("  after obs filter: %dx%d\n", nrow(M), ncol(M)))

# Build group key: either compartment4 (main) or "supp_<original>" for supp.
obs[, group := ifelse(
  keep_for_main, compartment4,
  paste0("supp_", tolower(gsub("[^A-Za-z0-9]+", "_", cell_type_original)))
)]

groups <- sort(unique(obs$group))
cat(sprintf("  groups to emit: %s\n", paste(groups, collapse = ", ")))


# --------------------------------------------------------------------
# 4. Aggregate: mean per (donor, group), at ENSG level
# --------------------------------------------------------------------
# Build a one-hot (cells x samples) sparse matrix W where
#   W[i, s] = 1 / n_cells(s)   if cell i belongs to sample s
# then   pseudobulk = M %*% W    is  genes x samples  dense-enough.

obs[, sample_id := paste(group, donor_id, sep = "__")]
samples <- obs[, .(
  group     = group[1],
  donor_id  = donor_id[1],
  biosample_id = biosample_id[1],
  rhythm    = rhythm[1],
  sex       = sex[1],
  n_cells   = .N
), by = sample_id]
setorder(samples, group, rhythm, donor_id)

cat(sprintf("  n samples (donor x group): %d\n", nrow(samples)))

# Sanity: gate check (redundant with r1 but cheap).
gate_fail <- samples[n_cells < MIN_CELLS]
if (nrow(gate_fail) > 0) {
  cat(sprintf("  WARNING %d (donor, ct) combos fail MIN_CELLS=%d; dropped:\n",
              nrow(gate_fail), MIN_CELLS))
  print(as.data.frame(gate_fail))
  samples <- samples[n_cells >= MIN_CELLS]
}

# Build W: 165740 cells x 204 samples, sparse.
s_idx <- match(obs$sample_id, samples$sample_id)
keep_cell <- !is.na(s_idx)
cat(sprintf("  cells with valid sample: %d / %d\n", sum(keep_cell), length(keep_cell)))

# compute 1 / n_cells per sample
inv_n <- 1 / samples$n_cells

W <- sparseMatrix(
  i = which(keep_cell),
  j = s_idx[keep_cell],
  x = inv_n[s_idx[keep_cell]],
  dims = c(ncol(M), nrow(samples))
)
colnames(W) <- samples$sample_id

cat("  aggregating mean per (donor, group) ... ")
t0 <- Sys.time()
PB <- M %*% W                        # genes x samples
PB <- as.matrix(PB)
t1 <- Sys.time()
cat(sprintf("done (%s)\n", format(difftime(t1, t0, units = "secs"))))
rownames(PB) <- genes$ensg
colnames(PB) <- samples$sample_id


# --------------------------------------------------------------------
# 5. Collapse ENSG -> SYMBOL (mean over ENSGs, rank-preserving)
# --------------------------------------------------------------------

collapse_to_symbol <- function(mat, ensg2sym) {
  # ensg2sym: character vector in the same order as rownames(mat)
  valid <- !is.na(ensg2sym) & nzchar(ensg2sym)
  mat <- mat[valid, , drop = FALSE]
  sym <- ensg2sym[valid]
  if (!anyDuplicated(sym)) {
    rownames(mat) <- sym
    return(mat)
  }
  # Group rows by symbol; take column-mean across duplicates.
  ord <- order(sym)
  sym <- sym[ord]
  mat <- mat[ord, , drop = FALSE]
  rle_s <- rle(sym)
  starts <- c(1, cumsum(rle_s$lengths) + 1)[seq_along(rle_s$values)]
  ends   <- cumsum(rle_s$lengths)
  out <- matrix(NA_real_, nrow = length(rle_s$values), ncol = ncol(mat),
                dimnames = list(rle_s$values, colnames(mat)))
  for (k in seq_along(rle_s$values)) {
    if (starts[k] == ends[k]) {
      out[k, ] <- mat[starts[k], ]
    } else {
      out[k, ] <- colMeans(mat[starts[k]:ends[k], , drop = FALSE])
    }
  }
  out
}

cat("collapsing ENSG -> SYMBOL ... ")
t0 <- Sys.time()
PB_sym <- collapse_to_symbol(PB, genes$symbol)
t1 <- Sys.time()
cat(sprintf("done (%d symbols, %s)\n", nrow(PB_sym),
            format(difftime(t1, t0, units = "secs"))))


# --------------------------------------------------------------------
# 6. Write per-group output
# --------------------------------------------------------------------

for (g in groups) {
  samps_g <- samples[group == g]
  if (nrow(samps_g) == 0) next
  outdir_g <- file.path(OUT_DIR, g)
  dir.create(outdir_g, recursive = TRUE, showWarnings = FALSE)

  expr <- PB_sym[, samps_g$sample_id, drop = FALSE]
  # header: gene then donor ids (sample_id is <group>__<donor>, but for the
  # final expression table we use the donor_id directly to keep columns tight)
  colnames(expr) <- samps_g$donor_id

  write.table(
    cbind(gene = rownames(expr), as.data.frame(expr)),
    file.path(outdir_g, "expr.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  fwrite(samps_g[, .(sample_id, donor_id, biosample_id, rhythm, sex, n_cells)],
         file.path(outdir_g, "samples.tsv"), sep = "\t")

  cat(sprintf("  wrote %s  (%d donors, %d genes)\n",
              sub(paste0("^", ROOT, "/"), "", outdir_g),
              nrow(samps_g), nrow(expr)))
}


# --------------------------------------------------------------------
# 7. Summary + QC outputs
# --------------------------------------------------------------------

summary_tbl <- samples[, .(
  n_samples   = .N,
  n_donors    = uniqueN(donor_id),
  n_AF        = uniqueN(donor_id[rhythm == "AF"]),
  n_Ctrl      = uniqueN(donor_id[rhythm == "Ctrl"]),
  median_cells = as.integer(median(n_cells)),
  min_cells   = min(n_cells),
  max_cells   = max(n_cells)
), by = group]
setorder(summary_tbl, group)
fwrite(summary_tbl, file.path(OUT_DIR, "summary.tsv"), sep = "\t")
cat("\n=== pseudobulk summary ===\n")
print(as.data.frame(summary_tbl), row.names = FALSE)

# qc_pseudobulk.tsv: records of gate failures (or confirmation that none).
qc <- if (nrow(gate_fail) > 0) gate_fail else
  data.table(sample_id = character(0), group = character(0),
             donor_id = character(0), biosample_id = character(0),
             rhythm = character(0), sex = character(0), n_cells = integer(0))
fwrite(qc, file.path(QC_DIR, "qc_pseudobulk.tsv"), sep = "\t")

cat("\nR2 DONE.\n")
