#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# D2 - RMA + symbol-collapse for GSE41177 (Affymetrix HG-U133 Plus 2.0).
#
# Same pipeline as Fig 2's scripts/p2_build_matrices.R:
#   - extract CELs from GSE41177_RAW.tar
#   - affy::ReadAffy + rma  (requires preprocessCore built with
#     --disable-threading; already installed for Fig 2)
#   - collapse probesets -> gene symbols via max-IQR-per-symbol rule
#
# Outputs (duration-cohort expression matrices):
#   Project_P4/GSE41177/expr_symbol_LAA.tsv
#   Project_P4/GSE41177/expr_symbol_LAPV.tsv
#   Project_P4/GSE41177/qc_pca.pdf           (coloured by rhythm, shaped by tissue)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(affy)
  library(hgu133plus2.db)
  library(AnnotationDbi)
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
ROOT <- find_root()
IN_DIR  <- file.path(ROOT, "Project_P4", "GSE41177")
TAR_PATH <- file.path(IN_DIR, "GSE41177_RAW.tar")
CEL_DIR  <- file.path(IN_DIR, "CEL")
META_PATH <- file.path(IN_DIR, "sample_metadata.tsv")
cat(sprintf("IN_DIR = %s\n", IN_DIR))

# --------------------------------------------------------------------
# 1. Extract CEL files
# --------------------------------------------------------------------

dir.create(CEL_DIR, recursive = TRUE, showWarnings = FALSE)
meta <- fread(META_PATH)
stopifnot(all(!is.na(meta$cel_filename)))

already <- list.files(CEL_DIR, pattern = "\\.CEL\\.gz$", full.names = FALSE)
if (!all(meta$cel_filename %in% already)) {
  cat("extracting CEL files from tar ...\n")
  # only extract .CEL.gz to avoid the .rma.chp.gz clutter
  system(paste0("tar -xf ", shQuote(TAR_PATH),
                " -C ", shQuote(CEL_DIR),
                " --wildcards '*.CEL.gz'"))
}

cel_paths <- file.path(CEL_DIR, meta$cel_filename)
stopifnot(all(file.exists(cel_paths)))
cat(sprintf("CEL files ready: %d\n", length(cel_paths)))


# --------------------------------------------------------------------
# 2. Affy RMA
# --------------------------------------------------------------------

cat(sprintf("running RMA on %d CELs ...\n", length(cel_paths)))
ab <- ReadAffy(filenames = cel_paths)
eset <- rma(ab, verbose = FALSE)
probe_mat <- exprs(eset)
colnames(probe_mat) <- meta$gsm
cat(sprintf("  RMA: %d probesets x %d samples\n",
            nrow(probe_mat), ncol(probe_mat)))


# --------------------------------------------------------------------
# 3. Probe -> symbol collapse (max-IQR rule, same as Fig 2)
# --------------------------------------------------------------------

hgu133plus2_probe2symbol <- function() {
  probesets <- keys(hgu133plus2.db, keytype = "PROBEID")
  m <- AnnotationDbi::select(hgu133plus2.db,
                             keys = probesets,
                             columns = c("PROBEID", "SYMBOL"),
                             keytype = "PROBEID")
  m <- m[!is.na(m$SYMBOL) & !duplicated(m$PROBEID), ]
  setNames(m$SYMBOL, m$PROBEID)
}
collapse_probes_maxIQR <- function(mat, probe2symbol) {
  common <- intersect(rownames(mat), names(probe2symbol))
  mat <- mat[common, , drop = FALSE]
  sym <- probe2symbol[common]
  iqr <- apply(mat, 1, IQR, na.rm = TRUE)
  keep_row <- tapply(seq_along(sym), sym,
                     function(idx) idx[which.max(iqr[idx])])
  keep_row <- unlist(keep_row, use.names = FALSE)
  out <- mat[keep_row, , drop = FALSE]
  rownames(out) <- sym[keep_row]
  out[order(rownames(out)), , drop = FALSE]
}

p2s <- hgu133plus2_probe2symbol()
sym_mat <- collapse_probes_maxIQR(probe_mat, p2s)
cat(sprintf("  collapsed to %d gene symbols\n", nrow(sym_mat)))


# --------------------------------------------------------------------
# 4. Split by tissue and write
# --------------------------------------------------------------------

for (tiss in c("LAA", "LAPV")) {
  keep <- meta$tissue == tiss
  stopifnot(sum(keep) == 19)
  m_t <- sym_mat[, meta$gsm[keep], drop = FALSE]
  out <- data.frame(gene = rownames(m_t), m_t, check.names = FALSE)
  out_path <- file.path(IN_DIR, sprintf("expr_symbol_%s.tsv", tiss))
  write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("wrote %s  (%d genes x %d samples)\n",
              sub(paste0("^", ROOT, "/"), "", out_path),
              nrow(m_t), ncol(m_t)))
}


# --------------------------------------------------------------------
# 5. QC PCA (both tissues in one PDF, shaped by tissue, coloured by rhythm)
# --------------------------------------------------------------------

qc_pdf <- file.path(IN_DIR, "qc_pca.pdf")
v <- apply(sym_mat, 1, var, na.rm = TRUE)
keep_g <- order(v, decreasing = TRUE)[seq_len(min(2000, length(v)))]
pr <- prcomp(t(sym_mat[keep_g, , drop = FALSE]), center = TRUE, scale. = FALSE)
pvar <- pr$sdev^2 / sum(pr$sdev^2)

mc <- meta[match(rownames(pr$x), meta$gsm), ]
col_rh <- ifelse(mc$rhythm == "SR", "#666666", "#8C510A")   # BrBG AF brown
pch_ts <- ifelse(mc$tissue == "LAA", 16, 17)                 # LAA circle, LAPV triangle

pdf(qc_pdf, width = 5.2, height = 4.4, useDingbats = FALSE)
op <- par(mar = c(4, 4, 3, 8), xpd = NA)
plot(pr$x[, 1], pr$x[, 2], col = col_rh, pch = pch_ts, cex = 1.1,
     xlab = sprintf("PC1 (%.1f%%)", 100 * pvar[1]),
     ylab = sprintf("PC2 (%.1f%%)", 100 * pvar[2]),
     main = "GSE41177 sample PCA", las = 1)
z <- scale(pr$x[, 1:2])
outl <- which(abs(z[, 1]) > 3 | abs(z[, 2]) > 3)
if (length(outl))
  text(pr$x[outl, 1], pr$x[outl, 2], labels = mc$donor_id[outl],
       cex = 0.6, pos = 3, offset = 0.3)
usr <- par("usr")
legend(x = usr[2] + 0.02 * diff(usr[1:2]), y = usr[4],
       legend = c("SR", "AF", "LAA", "LAPV"),
       col    = c("#666666", "#8C510A", "black", "black"),
       pch    = c(15, 15, 16, 17),
       bty = "n", cex = 0.75)
par(op)
dev.off()
cat(sprintf("wrote %s\n", sub(paste0("^", ROOT, "/"), "", qc_pdf)))

cat("\nD2 DONE.\n")
