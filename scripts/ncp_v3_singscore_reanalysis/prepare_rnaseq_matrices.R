#!/usr/bin/env Rscript
# Build symbol-indexed TMM/voom matrices from the approved GEO count files.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
ROOT <- dirname(dirname(dirname(script_path)))
DATA <- file.path(ROOT, "data", "ncp_v3_singscore_reanalysis")
META <- file.path(DATA, "metadata")
RAW <- file.path(DATA, "raw")
PROCESSED <- file.path(DATA, "processed")
RESULTS <- file.path(ROOT, "results", "ncp_v3_singscore_reanalysis")

build_one <- function(accession, count_filename) {
  meta_path <- file.path(META, paste0(accession, "_sample_metadata.tsv"))
  meta <- fread(meta_path)
  meta <- meta[analysis_inclusion == TRUE]
  stopifnot(nrow(meta) == if (accession == "GSE112339") 16L else 20L)
  stopifnot(!anyDuplicated(meta$sample_id))
  stopifnot(!anyDuplicated(meta$source_expression_column))

  counts_dt <- fread(file.path(RAW, count_filename), check.names = FALSE)
  raw_gene_id <- as.character(counts_dt[[1]])
  count_columns <- names(counts_dt)[-1]
  if (!all(meta$source_expression_column %in% count_columns)) {
    missing <- setdiff(meta$source_expression_column, count_columns)
    stop(accession, ": count columns missing: ", paste(missing, collapse = ", "))
  }
  source_columns <- meta$source_expression_column
  counts <- as.matrix(counts_dt[, ..source_columns])
  storage.mode(counts) <- "integer"
  rownames(counts) <- make.unique(raw_gene_id, sep = "_dup")
  colnames(counts) <- meta$sample_id

  group <- if (accession == "GSE112339") {
    factor(meta$chamber)
  } else {
    factor(paste(meta$chamber, meta$rhythm, sep = "_"))
  }
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, group = group)
  if (!any(keep)) stop(accession, ": filterByExpr retained no genes")
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  design <- model.matrix(~0 + group)
  voom_obj <- voom(dge, design, plot = FALSE)
  expr <- voom_obj$E

  symbols <- raw_gene_id[keep]
  if (anyDuplicated(symbols)) {
    row_variance <- apply(expr, 1, var, na.rm = TRUE)
    retained <- unlist(
      tapply(seq_along(symbols), symbols, function(idx) idx[which.max(row_variance[idx])]),
      use.names = FALSE
    )
    expr <- expr[retained, , drop = FALSE]
    symbols <- symbols[retained]
  }
  rownames(expr) <- symbols
  expr <- expr[order(rownames(expr)), , drop = FALSE]
  stopifnot(!anyDuplicated(rownames(expr)), all(is.finite(expr)))
  stopifnot(identical(colnames(expr), meta$sample_id))

  out_dir <- file.path(PROCESSED, accession)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(
    data.table(gene_symbol = rownames(expr), as.data.table(expr)),
    file.path(out_dir, "expr_symbol.tsv"),
    sep = "\t"
  )
  data.table(
    dataset = accession,
    n_raw_genes = nrow(counts),
    n_filterByExpr_genes = sum(keep),
    n_symbol_genes = nrow(expr),
    n_samples = ncol(expr),
    normalization = "edgeR TMM + limma voom",
    duplicate_rule = "post-voom maximum variance row per exact symbol"
  )
}

summary <- rbindlist(list(
  build_one("GSE112339", "GSE112339_HTSeq_hg19_Ensembl.txt.gz"),
  build_one("GSE128188", "GSE128188_CountsV2_EdgeR.txt.gz")
))
dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
fwrite(summary, file.path(RESULTS, "rnaseq_preprocessing_summary.tsv"), sep = "\t")

versions <- data.table(
  component = c("R", "edgeR", "limma"),
  version = c(as.character(getRversion()), as.character(packageVersion("edgeR")),
              as.character(packageVersion("limma")))
)
fwrite(versions, file.path(RESULTS, "rnaseq_software_versions.tsv"), sep = "\t")
cat("Wrote RNA-seq symbol matrices for GSE112339 and GSE128188\n")
