#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# R2b — Donor x compartment pseudobulk for GSE255612 with the EXTENDED
#       6-compartment scheme that includes adipocytes and epicardial cells
#       and a CLEAN vascular-endothelial definition (Endothelial_1 +
#       Endothelial_2 only; lymphatic split off as a separate group).
#
# Rationale (post hoc Fig 1 framework upgrade):
#   - The previous 4-compartment scheme excluded Adipocytes and
#     Epicardial_Cells, both of which are biologically central to the
#     N1 data-driven modules M7 (adipose-program, ADIPOQ/LEP/PLIN1)
#     and M11 (PITX2-anchored mesothelial / epicardial; ITLN1, MSLN,
#     KRT19, PRG4, SLPI).  Attribution and carrier validation cannot
#     work without these compartments present.
#   - Endothelial_or_vascular previously mixed lymphatic endothelial in
#     with vascular endothelial; "endothelial" is now defined cleanly
#     as Endothelial_1 + Endothelial_2 (vascular only).
#
# Compartment scheme (NEW):
#     atrialCM        Cardiomyocytes                       n = 32,438
#     fibroblast      Fibroblasts                          n = 56,258
#     endothelial     Endothelial_1 + Endothelial_2        n = 25,763
#     immune          Macrophages + Lymphocytes +
#                     Mast_Cells + B_Cells                 n = 26,674
#     adipocyte       Adipocytes                           n = 11,243   <- NEW
#     epicardial      Epicardial_Cells                     n =  1,051   <- NEW
#     [supp_lymphatic] Lymphatic_Endothelial               n =  1,735   (supp)
#     [supp_pericyte / supp_VSMC]                                       (supp)
#
# Output layout (all under Project_P3/GSE255612/pseudobulk_6comp/):
#     atrialCM/expr.tsv + samples.tsv
#     fibroblast/expr.tsv + samples.tsv
#     endothelial/expr.tsv + samples.tsv      <- recomputed (clean)
#     immune/expr.tsv + samples.tsv
#     adipocyte/expr.tsv + samples.tsv        <- NEW
#     epicardial/expr.tsv + samples.tsv       <- NEW
#     supp_lymphatic_endothelial/...          <- separated out
#     summary.tsv
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
  cwd <- getwd(); if (basename(cwd) == "scripts") return(dirname(cwd)); cwd
}
ROOT    <- find_root()
IN_DIR  <- file.path(ROOT, "Project_P3", "GSE255612")
OUT_DIR <- file.path(IN_DIR, "pseudobulk_full")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("ROOT    = %s\n", ROOT))
cat(sprintf("OUT_DIR = %s\n", OUT_DIR))

MIN_CELLS <- 30L

# ---------------------------------------------------------------------------
# 1. Extended compartment mapping (cell_type_original -> compartment6_clean)
# ---------------------------------------------------------------------------
# Each cell type kept SEPARATE, with two exceptions:
#   - Endothelial_1 + Endothelial_2 = "endothelial" (clean vascular endothelial)
#   - Megakaryocytes excluded (n = 58, too sparse for any donor to pass gate)
COMPARTMENT_MAP <- c(
  "Cardiomyocytes"          = "atrialCM",
  "Fibroblasts"             = "fibroblast",
  "Endothelial_1"           = "endothelial",
  "Endothelial_2"           = "endothelial",
  "Lymphatic_Endothelial"   = "lymphatic_endothelial",
  "Macrophages"             = "macrophage",
  "Lymphocytes"             = "lymphocyte",
  "Mast_Cells"              = "mast_cell",
  "B_Cells"                 = "B_cell",
  "Pericytes"               = "pericyte",
  "Vascular_Smooth_Muscle"  = "vascular_smooth_muscle",
  "Adipocytes"              = "adipocyte",
  "Epicardial_Cells"        = "epicardial",
  "Neuronal"                = "neuronal",
  "Megakaryocytes"          = NA_character_   # too sparse (n = 58)
)

# ---------------------------------------------------------------------------
# 2. Load barcodes / genes; build obs from RAW metadata (project's pre-derived
#    obs file had already filtered out Adipocytes and Epicardial_Cells, so
#    we go back to the source of truth here).
# ---------------------------------------------------------------------------
cat("reading barcodes / genes / raw MetaData ...\n")
barcodes <- fread(file.path(IN_DIR, "GSE255612_AF_snRNA_Processed_Expression_Matrix_barcodes_V1.tsv.gz"),
                  header = FALSE, col.names = "barcode")$barcode
genes <- fread(file.path(IN_DIR, "GSE255612_AF_snRNA_Processed_Expression_Matrix_genes_final.tsv.gz"),
               header = FALSE, col.names = c("ensg", "symbol"))
# Raw metadata (header row 1, type row 2 starting from 'TYPE'); skip type row.
md_raw <- fread(file.path(IN_DIR, "GSE255612_AF_snRNA_MetaData.txt.gz"),
                skip = 2, header = FALSE)
md_header <- fread(file.path(IN_DIR, "GSE255612_AF_snRNA_MetaData.txt.gz"),
                   nrows = 1, header = FALSE)
setnames(md_raw, as.character(md_header))
setnames(md_raw, "NAME", "cell_id")
setnames(md_raw, "cell_type", "cell_type_original")
setnames(md_raw, "af", "rhythm_raw")
md_raw[, rhythm := ifelse(rhythm_raw == "case", "AF",
                  ifelse(rhythm_raw == "control", "Ctrl", rhythm_raw))]
obs <- md_raw[, .(cell_id, biosample_id, donor_id, cell_type_original,
                  sex, rhythm)]
cat(sprintf("  raw metadata: %d cells\n", nrow(obs)))

obs[, compartment6 := COMPARTMENT_MAP[cell_type_original]]
obs <- obs[!is.na(compartment6)]
# restrict to cells present in the matrix barcode order
obs <- obs[cell_id %in% barcodes]
cat(sprintf("  cells with compartment6 assignment present in matrix: %d\n",
            nrow(obs)))

cat("\n=== compartment6 cell counts ===\n")
print(obs[, .N, by = compartment6][order(-N)])

# ---------------------------------------------------------------------------
# 3. Load mtx and restrict to retained cells
# ---------------------------------------------------------------------------
cat("\nreading mtx ...\n")
t0 <- Sys.time()
M <- readMM(gzfile(file.path(IN_DIR, "GSE255612_AF_snRNA_Matrix_V1.mtx.gz")))
cat(sprintf("  mtx: %dx%d  (%s)\n", nrow(M), ncol(M),
            format(difftime(Sys.time(), t0, units = "secs"))))
stopifnot(nrow(M) == nrow(genes), ncol(M) == length(barcodes))
M <- as(M, "CsparseMatrix")
rownames(M) <- genes$ensg
colnames(M) <- barcodes

keep_idx <- match(obs$cell_id, colnames(M))
stopifnot(!anyNA(keep_idx))
M <- M[, keep_idx, drop = FALSE]
cat(sprintf("  after obs filter: %dx%d\n", nrow(M), ncol(M)))

# ---------------------------------------------------------------------------
# 4. Build (donor, compartment6) sample table + gate
# ---------------------------------------------------------------------------
obs[, sample_id := paste(compartment6, donor_id, sep = "__")]
samples <- obs[, .(
  group        = compartment6[1],
  donor_id     = donor_id[1],
  biosample_id = biosample_id[1],
  rhythm       = rhythm[1],
  sex          = sex[1],
  n_cells      = .N
), by = sample_id]
setorder(samples, group, rhythm, donor_id)

gate_fail <- samples[n_cells < MIN_CELLS]
if (nrow(gate_fail) > 0) {
  cat(sprintf("\n  WARNING: %d (donor x compartment) drop n_cells < %d:\n",
              nrow(gate_fail), MIN_CELLS))
  print(gate_fail[, .N, by = group])
  samples <- samples[n_cells >= MIN_CELLS]
}
cat(sprintf("\n  retained (donor x compartment) samples: %d\n", nrow(samples)))

# ---------------------------------------------------------------------------
# 5. Aggregate via sparse one-hot weight matrix
# ---------------------------------------------------------------------------
s_idx <- match(obs$sample_id, samples$sample_id)
keep_cell <- !is.na(s_idx)
inv_n <- 1 / samples$n_cells
W <- sparseMatrix(
  i = which(keep_cell),
  j = s_idx[keep_cell],
  x = inv_n[s_idx[keep_cell]],
  dims = c(ncol(M), nrow(samples))
)
colnames(W) <- samples$sample_id

cat("aggregating mean per (donor, compartment6) ... ")
t0 <- Sys.time()
PB <- as.matrix(M %*% W)
cat(sprintf("done (%s)\n", format(difftime(Sys.time(), t0, units = "secs"))))
rownames(PB) <- genes$ensg
colnames(PB) <- samples$sample_id

# ---------------------------------------------------------------------------
# 6. Collapse ENSG -> SYMBOL (mean across duplicates; rank-preserving)
# ---------------------------------------------------------------------------
collapse_to_symbol <- function(mat, ensg2sym) {
  valid <- !is.na(ensg2sym) & nzchar(ensg2sym)
  mat <- mat[valid, , drop = FALSE]
  sym <- ensg2sym[valid]
  if (!anyDuplicated(sym)) {
    rownames(mat) <- sym; return(mat)
  }
  ord <- order(sym); sym <- sym[ord]; mat <- mat[ord, , drop = FALSE]
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
cat(sprintf("done (%d symbols, %s)\n", nrow(PB_sym),
            format(difftime(Sys.time(), t0, units = "secs"))))

# ---------------------------------------------------------------------------
# 7. Write per-group output
# ---------------------------------------------------------------------------
groups <- sort(unique(samples$group))
for (g in groups) {
  samps_g <- samples[group == g]
  if (nrow(samps_g) == 0) next
  outdir_g <- file.path(OUT_DIR, g)
  dir.create(outdir_g, recursive = TRUE, showWarnings = FALSE)
  expr <- PB_sym[, samps_g$sample_id, drop = FALSE]
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

summary_tbl <- samples[, .(
  n_samples    = .N,
  n_donors     = uniqueN(donor_id),
  n_AF         = uniqueN(donor_id[rhythm == "AF"]),
  n_Ctrl       = uniqueN(donor_id[rhythm == "Ctrl"]),
  median_cells = as.integer(median(n_cells)),
  min_cells    = min(n_cells),
  max_cells    = max(n_cells)
), by = group]
setorder(summary_tbl, group)
fwrite(summary_tbl, file.path(OUT_DIR, "summary.tsv"), sep = "\t")
cat("\n=== pseudobulk summary (6-compartment) ===\n")
print(as.data.frame(summary_tbl), row.names = FALSE)

cat("\nR2b DONE.\n")
