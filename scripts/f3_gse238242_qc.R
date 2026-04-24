#!/usr/bin/env Rscript
# F3 — Minimal QC for GSE238242: metadata summary + paired RNA/ATAC file check.
#     One RNA file peek: matrix orientation and sparsity (first 500 cells if large).
#
# Inputs:
#   Project_P5/GSE238242/GSE238242_snAF.metadata.tsv.gz
#   Project_P5/GSE238242/extracted/*.counts.tsv.gz
#
# Output:
#   Project_P5/GSE238242/qc/qc_summary.tsv
#   Project_P5/GSE238242/qc/qc_report.txt
suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    message("Installing data.table ...")
    install.packages("data.table", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
ROOT <- normalizePath(file.path(dirname(script), ".."))
GDIR <- file.path(ROOT, "Project_P5", "GSE238242")
META <- file.path(GDIR, "GSE238242_snAF.metadata.tsv.gz")
EXT <- file.path(GDIR, "extracted")
QC <- file.path(GDIR, "qc")
dir.create(QC, showWarnings = FALSE, recursive = TRUE)

report <- function(...) cat(..., file = tf, append = TRUE, sep = "")

tf <- file.path(QC, "qc_report.txt")
cat("", file = tf)

rpt <- function(...) {
  cat(..., "\n")
  cat(..., file = tf, append = TRUE, "\n")
}

rpt("=== GSE238242 QC ", format(Sys.time(), "%Y-%m-%d %H:%M"), " ===\n")

if (!file.exists(META)) {
  stop("Missing metadata: ", META, "\nRun: bash scripts/f1_gse238242_download.sh")
}
md <- fread(META, showProgress = FALSE, header = FALSE, skip = 1L)
setnames(md, c("barcode", "cell_type", "sample", "sex", "Rhythm"))

rpt("Metadata rows (nuclei/barcode): ", nrow(md))
rpt("Columns: ", paste(names(md), collapse = ", "))

if (!"sample" %in% names(md)) stop("Expected column: sample")
if (!"Rhythm" %in% names(md) && !"rhythm" %in% names(md)) {
  warning("Rhythm column not found")
} else {
  rhythm_col <- if ("Rhythm" %in% names(md)) "Rhythm" else "rhythm"
  rpt("\nRhythm counts:")
  print(table(md[[rhythm_col]]))
  capture.output(print(table(md[[rhythm_col]])), file = tf, append = TRUE)
}

if ("cell_type" %in% names(md)) {
  rpt("\nCell type counts (top 15):")
  ct <- sort(table(md$cell_type), decreasing = TRUE)
  print(head(ct, 15))
  capture.output(print(head(ct, 15)), file = tf, append = TRUE)
}

## Per-sample nucleus counts
sm <- md[, .N, by = sample][order(sample)]
fwrite(sm, file.path(QC, "qc_nuclei_per_sample.tsv"), sep = "\t")
rpt("\nWrote qc_nuclei_per_sample.tsv")

## File pairing
if (!dir.exists(EXT)) {
  rpt("\nNo extracted/ directory — run: bash scripts/f2_gse238242_extract.sh")
  quit(save = "no", status = 0)
}

all_f <- list.files(EXT, pattern = "counts\\.tsv\\.gz$", full.names = TRUE)
rna <- grep("RNA\\.counts", all_f, value = TRUE)
atac <- grep("ATAC\\.counts", all_f, value = TRUE)
rpt("\nRNA files: ", length(rna))
rpt("ATAC files: ", length(atac))

donor_id <- function(path) {
  m <- regexpr("CF[0-9]+", basename(path), perl = TRUE)
  regmatches(basename(path), m)
}
pair <- rbindlist(lapply(rna, function(rp) {
  id <- donor_id(rp)
  m_at <- atac[grepl(paste0(id, "ATAC"), basename(atac))]
  data.table(
    donor = id,
    rna_file = basename(rp),
    atac_file = if (length(m_at)) basename(m_at[1]) else NA_character_,
    paired = length(m_at) > 0
  )
}))
fwrite(pair[order(donor)], file.path(QC, "qc_rna_atac_file_pairs.tsv"), sep = "\t")
rpt("Wrote qc_rna_atac_file_pairs.tsv")

## Peek first RNA matrix (smallest file if possible)
if (length(rna)) {
  sz <- file.info(rna)$size
  pick <- rna[which.min(sz)]
  rpt("\nPeek matrix (smallest RNA file): ", basename(pick))
  a <- fread(
    pick,
    nrows = 5L,
    header = TRUE,
    showProgress = FALSE
  )
  nc <- ncol(a)
  nr <- nrow(a)
  rpt("First-pass fread nrows=5: nrow=", nr, " ncol=", nc)
  rpt("Column 1 name: ", names(a)[1])
  rpt("Sample colnames (cells?): ", paste(head(names(a), 4), collapse = ", "), " ...")

  peek_info <- data.table(
    file = basename(pick),
    nrow_peek = nr,
    ncol = nc,
    first_col = names(a)[1],
    likely_genes_in_rows = isTRUE(!grepl("^[0-9]+-[A-Z]", names(a)[2]))
  )
  fwrite(peek_info, file.path(QC, "qc_matrix_peek.tsv"), sep = "\t")
}

rpt("\nDone. See ", QC)
