#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# P2 — Build gene × sample expression matrices for all three cohorts.
#
# Outputs (all in log-space, symbol-indexed, columns in the order the sample
# metadata TSV prescribes):
#
#   Project_P1/GSE128188/expr_symbol.tsv    voom logCPM, 20 samples
#   Project_P1/GSE79768/expr_symbol.tsv     RMA log2,    26 samples
#   Project_P1/GSE115574/expr_symbol.tsv    RMA log2,    59 samples
#
# Plus per-cohort QC PDFs with sample PCA (colour = chamber, shape = rhythm):
#
#   Project_P1/<cohort>/qc_pca.pdf
#
# Rules locked for paired-chamber aCAI / aCAI_rel construction:
#   * GSE128188  edgeR::filterByExpr + voom -> logCPM
#   * GSE79768   affy::ReadAffy + rma       -> log2
#   * GSE115574  affy::ReadAffy + rma       -> log2
#   * probe -> gene symbol collapse rule: max-IQR probeset per symbol
#   * NO cross-cohort merging, NO batch correction
#   * column order follows sample_metadata.tsv so downstream code never
#     has to re-match anything.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(affy)
  library(hgu133plus2.db)
  library(AnnotationDbi)
})

# Robustly locate the project root. Works when run:
#   * from project root:     `Rscript scripts/p2_build_matrices.R`
#   * from scripts/ dir:     `Rscript p2_build_matrices.R`
#   * via an absolute path:  `Rscript /abs/path/scripts/p2_build_matrices.R`
find_root <- function() {
  # 1. command-line invocation
  ca <- commandArgs(trailingOnly = FALSE)
  flag <- "--file="
  f <- sub(flag, "", grep(flag, ca, value = TRUE)[1])
  if (length(f) && !is.na(f) && nzchar(f)) {
    d <- normalizePath(dirname(f), mustWork = FALSE)
    if (basename(d) == "scripts") return(dirname(d))
    return(d)
  }
  # 2. fall back to cwd (works when pasted interactively)
  cwd <- getwd()
  if (basename(cwd) == "scripts") return(dirname(cwd))
  cwd
}
ROOT <- find_root()

P1 <- file.path(ROOT, "Project_P1")
cat(sprintf("ROOT = %s\n", ROOT))
stopifnot(dir.exists(P1))


# --------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------

load_meta <- function(cohort) {
  p <- file.path(P1, cohort, "sample_metadata.tsv")
  read.delim(p, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}

# Collapse a probe-level matrix to symbol-level by taking, for each symbol,
# the probeset with the largest IQR across samples.
collapse_probes_maxIQR <- function(mat, probe2symbol) {
  stopifnot(is.matrix(mat))
  common <- intersect(rownames(mat), names(probe2symbol))
  mat <- mat[common, , drop = FALSE]
  sym <- probe2symbol[common]
  iqr <- apply(mat, 1, IQR, na.rm = TRUE)
  keep_row <- tapply(seq_along(sym), sym, function(idx) idx[which.max(iqr[idx])])
  keep_row <- unlist(keep_row, use.names = FALSE)
  out <- mat[keep_row, , drop = FALSE]
  rownames(out) <- sym[keep_row]
  out[order(rownames(out)), , drop = FALSE]
}

# Extract the per-cohort tar into a CEL/ subdir only if it hasn't been
# extracted yet, and return the list of local CEL paths in the order given
# by `expected_basenames`.
extract_cels <- function(cohort, tar_name, expected_basenames) {
  tar_path <- file.path(P1, cohort, tar_name)
  cel_dir  <- file.path(P1, cohort, "CEL")
  dir.create(cel_dir, recursive = TRUE, showWarnings = FALSE)
  present <- list.files(cel_dir, pattern = "\\.CEL(\\.gz)?$", full.names = FALSE)
  if (!all(expected_basenames %in% present)) {
    cat(sprintf("  extracting %s -> %s\n", tar_name, cel_dir))
    untar(tar_path, exdir = cel_dir)
  } else {
    cat(sprintf("  CEL files already extracted in %s\n", cel_dir))
  }
  file.path(cel_dir, expected_basenames)
}

# Run RMA on a fixed-order vector of CEL paths and return a probeset × sample
# log2 intensity matrix with column names == GSMs.
rma_matrix <- function(cel_paths, gsm_names) {
  stopifnot(length(cel_paths) == length(gsm_names))
  stopifnot(all(file.exists(cel_paths)))
  ab <- ReadAffy(filenames = cel_paths)
  eset <- rma(ab, verbose = FALSE)
  mat <- exprs(eset)
  colnames(mat) <- gsm_names
  mat
}

# Build probeset -> SYMBOL mapping for HG-U133 Plus 2.0.  AnnotationDbi
# returns NA where no symbol is available; those probesets are dropped.
hgu133plus2_probe2symbol <- function() {
  probesets <- keys(hgu133plus2.db, keytype = "PROBEID")
  map <- AnnotationDbi::select(
    hgu133plus2.db,
    keys = probesets,
    columns = c("PROBEID", "SYMBOL"),
    keytype = "PROBEID"
  )
  # One row per probeset; drop NAs and duplicate mappings.
  map <- map[!is.na(map$SYMBOL) & !duplicated(map$PROBEID), ]
  setNames(map$SYMBOL, map$PROBEID)
}

# Minimal PCA QC plot, written as a vector PDF.
write_pca_pdf <- function(mat, meta, out_pdf) {
  # top 2000 most variable genes
  v <- apply(mat, 1, var, na.rm = TRUE)
  v[is.na(v)] <- 0
  keep <- order(v, decreasing = TRUE)[seq_len(min(2000, length(v)))]
  pr <- prcomp(t(mat[keep, , drop = FALSE]), center = TRUE, scale. = FALSE)
  pvar <- pr$sdev^2 / sum(pr$sdev^2)

  m <- meta[match(rownames(pr$x), meta$gsm), , drop = FALSE]
  ch_col <- ifelse(m$chamber == "LA", "#2166AC", "#B2182B")  # LA=blue, RA=red
  rh_pch <- ifelse(m$rhythm == "SR", 16, 17)                 # SR=circle, AF=triangle

  pdf(out_pdf, width = 5, height = 4.2, useDingbats = FALSE)
  op <- par(mar = c(4, 4, 3, 8), xpd = NA)
  plot(pr$x[, 1], pr$x[, 2],
       col = ch_col, pch = rh_pch, cex = 1.1,
       xlab = sprintf("PC1 (%.1f%%)", 100 * pvar[1]),
       ylab = sprintf("PC2 (%.1f%%)", 100 * pvar[2]),
       main = sprintf("%s sample PCA (top 2000 var)", m$cohort[1]),
       las = 1)
  # label outliers (points >3 SD from centroid on PC1 or PC2)
  z <- scale(pr$x[, 1:2])
  outl <- which(abs(z[, 1]) > 3 | abs(z[, 2]) > 3)
  if (length(outl)) {
    text(pr$x[outl, 1], pr$x[outl, 2], labels = m$donor_id[outl],
         cex = 0.6, pos = 3, offset = 0.3)
  }
  # legend outside
  usr <- par("usr")
  legend(x = usr[2] + 0.02 * diff(usr[1:2]), y = usr[4],
         legend = c("LA", "RA", "SR", "AF"),
         col    = c("#2166AC", "#B2182B", "black", "black"),
         pch    = c(15, 15, 16, 17),
         bty    = "n", cex = 0.75)
  par(op)
  dev.off()
}

`%||%` <- function(a, b) if (is.null(a)) b else a


# --------------------------------------------------------------------
# GSE128188 (bulk RNA-seq, counts already symbol-indexed)
# --------------------------------------------------------------------

build_gse128188 <- function() {
  cohort <- "GSE128188"
  cat(sprintf("\n=== %s (bulk RNA-seq) ===\n", cohort))
  meta <- load_meta(cohort)
  stopifnot(nrow(meta) == 20)

  counts_path <- file.path(P1, cohort, "GSE128188_CountsV2_EdgeR.txt.gz")
  raw <- read.delim(gzfile(counts_path), sep = "\t", stringsAsFactors = FALSE,
                    check.names = FALSE, row.names = 1)
  cat(sprintf("  raw counts: %d genes x %d samples\n", nrow(raw), ncol(raw)))

  # Column order = meta$count_column, then rename to meta$gsm so downstream is GSM-keyed.
  stopifnot(all(meta$count_column %in% colnames(raw)))
  raw <- raw[, meta$count_column, drop = FALSE]
  colnames(raw) <- meta$gsm

  # edgeR filterByExpr + voom
  dge <- DGEList(counts = as.matrix(raw))
  group <- paste(meta$chamber, meta$rhythm, sep = "_")
  keep <- filterByExpr(dge, group = group)
  cat(sprintf("  filterByExpr kept %d / %d genes\n", sum(keep), length(keep)))
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")

  design <- model.matrix(~0 + group)
  colnames(design) <- levels(factor(group))
  v <- voom(dge, design, plot = FALSE)
  mat <- v$E
  cat(sprintf("  voom logCPM: %d genes x %d samples\n", nrow(mat), ncol(mat)))

  out_tsv <- file.path(P1, cohort, "expr_symbol.tsv")
  write.table(cbind(gene = rownames(mat), mat), out_tsv,
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  wrote %s\n", out_tsv))

  write_pca_pdf(mat, meta, file.path(P1, cohort, "qc_pca.pdf"))
  cat(sprintf("  wrote %s\n", file.path(P1, cohort, "qc_pca.pdf")))
}


# --------------------------------------------------------------------
# GSE79768 / GSE115574 (Affymetrix arrays)
# --------------------------------------------------------------------

build_array_cohort <- function(cohort, tar_name) {
  cat(sprintf("\n=== %s (Affymetrix HG-U133 Plus 2.0) ===\n", cohort))
  meta <- load_meta(cohort)

  cel_paths <- extract_cels(cohort, tar_name, meta$cel_file)

  cat(sprintf("  running RMA on %d CEL files...\n", length(cel_paths)))
  probe_mat <- rma_matrix(cel_paths, meta$gsm)
  cat(sprintf("  RMA: %d probesets x %d samples\n", nrow(probe_mat), ncol(probe_mat)))

  p2s <- hgu133plus2_probe2symbol()
  sym_mat <- collapse_probes_maxIQR(probe_mat, p2s)
  cat(sprintf("  collapsed (max-IQR) to %d gene symbols\n", nrow(sym_mat)))

  out_tsv <- file.path(P1, cohort, "expr_symbol.tsv")
  write.table(cbind(gene = rownames(sym_mat), sym_mat), out_tsv,
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  wrote %s\n", out_tsv))

  write_pca_pdf(sym_mat, meta, file.path(P1, cohort, "qc_pca.pdf"))
  cat(sprintf("  wrote %s\n", file.path(P1, cohort, "qc_pca.pdf")))
}


# --------------------------------------------------------------------
# main
# --------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args)) args[1] else "all"

if (target %in% c("all", "GSE128188"))  build_gse128188()
if (target %in% c("all", "GSE79768"))   build_array_cohort("GSE79768",  "GSE79768_RAW.tar")
if (target %in% c("all", "GSE115574"))  build_array_cohort("GSE115574", "GSE115574_RAW.tar")

cat("\nDONE.\n")
