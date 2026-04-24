#!/usr/bin/env Rscript
# F4b — Minimal sanity check on F4 pseudobulks.
# Two QC plots only (NOT for main figure):
#   qc_pseudobulk_pca.pdf          PCA, color = compartment, shape = rhythm, label = donor
#   qc_pseudobulk_correlation.pdf  Spearman correlation heatmap, annotated by compartment / rhythm
#
# Inputs: gse238242_pseudobulk_counts_main.rds + meta_main.tsv (from F4).
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(edgeR)
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
ROOT <- find_root()
PB   <- file.path(ROOT, "Project_P5", "GSE238242", "pseudobulk")
QC   <- file.path(ROOT, "Project_P5", "GSE238242", "qc")
dir.create(QC, recursive = TRUE, showWarnings = FALSE)

M <- as.matrix(readRDS(file.path(PB, "gse238242_pseudobulk_counts_main.rds")))
meta <- fread(file.path(PB, "gse238242_pseudobulk_meta_main.tsv"))
stopifnot(all(meta$pb_id %in% colnames(M)))
M <- M[, meta$pb_id]
cat(sprintf("loaded: %d genes x %d pb\n", nrow(M), ncol(M)))

# TMM + logCPM (per-pb library size)
y  <- DGEList(counts = M)
keep <- filterByExpr(y, group = meta$compartment, min.count = 5)
cat(sprintf("filterByExpr keeps %d / %d genes\n", sum(keep), nrow(y)))
y  <- y[keep, , keep.lib.sizes = FALSE]
y  <- calcNormFactors(y, method = "TMM")
lcpm <- cpm(y, log = TRUE, prior.count = 1)

# ---- PCA ----
top_v  <- order(matrixStats::rowVars(lcpm), decreasing = TRUE)[1:min(2000L, nrow(lcpm))]
pc     <- prcomp(t(lcpm[top_v, ]), scale. = TRUE)
var_pc <- round(100 * (pc$sdev^2 / sum(pc$sdev^2))[1:2], 1)

ct_col <- c(atrialCM = "#E6AB02", fibroblast = "#66A61E",
            endothelial_or_vascular = "#7570B3", immune = "#D95F02")
sh     <- ifelse(meta$rhythm == "AF", 17L, 16L)   # 17=triangle, 16=circle

pdf(file.path(QC, "qc_pseudobulk_pca.pdf"), width = 5.4, height = 4.4)
op <- par(mar = c(4, 4, 2.4, 1), cex = 0.75, las = 1)
plot(pc$x[, 1], pc$x[, 2],
     col = ct_col[meta$compartment], pch = sh, cex = 1.4,
     xlab = sprintf("PC1 (%.1f%%)", var_pc[1]),
     ylab = sprintf("PC2 (%.1f%%)", var_pc[2]),
     main = "GSE238242 pseudobulk PCA (TMM logCPM, top 2k var)")
text(pc$x[, 1], pc$x[, 2], labels = meta$donor, pos = 3, cex = 0.55, col = "#333333")
legend("topleft", names(ct_col), col = ct_col, pch = 16, bty = "n", cex = 0.7,
       title = "compartment", title.adj = 0)
legend("topright", c("AF", "SR"), pch = c(17, 16), col = "#444444",
       bty = "n", cex = 0.7, title = "rhythm", title.adj = 0)
par(op); dev.off()
cat("wrote ", file.path(QC, "qc_pseudobulk_pca.pdf"), "\n")

# ---- Sample correlation (Spearman) ----
S   <- cor(lcpm, method = "spearman")
ord <- order(meta$compartment, meta$rhythm, meta$donor)
S   <- S[ord, ord]; mo <- meta[ord]

pdf(file.path(QC, "qc_pseudobulk_correlation.pdf"), width = 6.6, height = 6.0)
op <- par(mar = c(0.5, 0.5, 2.4, 0.5))
n  <- ncol(S)
pal <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
zr  <- range(S, na.rm = TRUE)
image(1:n, 1:n, S, col = pal, xaxt = "n", yaxt = "n",
      xlab = "", ylab = "",
      main = sprintf("GSE238242 pseudobulk Spearman correlation (range %.2f-%.2f)", zr[1], zr[2]))

# Compartment + rhythm colored bars
top_y <- n + 0.6
left_x <- 0.4
rect(seq_len(n) - 0.5, top_y, seq_len(n) + 0.5, top_y + 0.6,
     col = ct_col[mo$compartment], border = NA, xpd = TRUE)
rect(left_x - 1.0, seq_len(n) - 0.5, left_x, seq_len(n) + 0.5,
     col = ct_col[mo$compartment], border = NA, xpd = TRUE)
rect(seq_len(n) - 0.5, top_y + 0.6, seq_len(n) + 0.5, top_y + 1.2,
     col = ifelse(mo$rhythm == "AF", "#666666", "#cccccc"), border = NA, xpd = TRUE)

# Donor labels (small, rotated)
text(seq_len(n), 0, labels = mo$donor, srt = 90, pos = 2, cex = 0.45,
     offset = -1.0, xpd = TRUE)
text(0, seq_len(n), labels = mo$donor, pos = 2, cex = 0.45,
     offset = 1.5, xpd = TRUE)
par(op); dev.off()
cat("wrote ", file.path(QC, "qc_pseudobulk_correlation.pdf"), "\n")

cat("F4b DONE.\n")
