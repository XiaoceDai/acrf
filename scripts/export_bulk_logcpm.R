#!/usr/bin/env Rscript
# Export voom-normalized logCPM matrix from GSE112339 bulk data.
#
# Rebuilds the same LA/RA curated 16-sample DGE, filters, TMM-normalizes,
# and runs voom to get fit_full$v$E. Writes that matrix as TSV.
# Also writes the curated sample metadata for downstream plotting.
#
# Inputs (defaults relative to script's working directory):
#   Project_P01/GSE112339/GSE112339_HTSeq_hg19_Ensembl.txt.gz
#   Project_P01/GSE112339/metadata.xlsx
#
# Outputs:
#   Project_P02/bulk/logcpm_voom.tsv
#   Project_P02/bulk/sample_meta.tsv

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(readxl)
})

here <- function(...) file.path("/home/bio/Project_P02", ...)

count_path    <- here("Project_P01/GSE112339/GSE112339_HTSeq_hg19_Ensembl.txt.gz")
metadata_path <- here("Project_P01/GSE112339/metadata.xlsx")
out_dir       <- here("bulk")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

normalize_chamber <- function(x) {
  x <- toupper(trimws(as.character(x)))
  ans <- rep(NA_character_, length(x))
  ans[x %chin% c("LA", "LEFT ATRIUM", "LEFT_ATRIUM", "LEFT-ATRIUM")]  <- "LA"
  ans[x %chin% c("RA", "RIGHT ATRIUM", "RIGHT_ATRIUM", "RIGHT-ATRIUM")] <- "RA"
  ans
}

message("[1/4] Reading count matrix")
counts_dt <- fread(count_path)
gene_id_col <- names(counts_dt)[1]
sample_ids  <- names(counts_dt)[-1]

message("[2/4] Building metadata and enforcing paired LA/RA donors")
meta <- data.table(sample_id = sample_ids)
meta[, chamber := sub("_.*$", "", sample_id)]
meta[, donor   := sub("^[^_]+_", "", sample_id)]
meta <- meta[chamber %chin% c("LA", "RA")]

# Keep donors with both LA and RA
pair_tab <- dcast(meta, donor ~ chamber, value.var = "sample_id",
                  fun.aggregate = length, fill = 0)
valid_donor <- pair_tab[LA == 1 & RA == 1, donor]
meta <- meta[donor %chin% valid_donor]
setorder(meta, donor, chamber)

stopifnot(length(unique(meta$donor)) == 8, nrow(meta) == 16)

meta[, donor   := factor(donor)]
meta[, chamber := factor(chamber, levels = c("RA", "LA"))]

sample_use <- meta$sample_id

raw_gene_id <- as.character(counts_dt[[gene_id_col]])
gene_key <- make.unique(raw_gene_id, sep = "_dup")
gene_map <- data.table(
  gene_key   = gene_key,
  gene_id_raw = raw_gene_id,
  ensembl_id  = ifelse(grepl("^ENSG", raw_gene_id, ignore.case = TRUE),
                       sub("\\..*$", "", raw_gene_id), raw_gene_id),
  gene_symbol = ifelse(grepl("^ENSG", raw_gene_id, ignore.case = TRUE),
                       NA_character_, raw_gene_id)
)

counts_mat <- as.matrix(counts_dt[, ..sample_use])
rownames(counts_mat) <- gene_key
storage.mode(counts_mat) <- "integer"

dge <- DGEList(counts = counts_mat)

message("[3/4] Filtering lowly expressed genes + TMM + voom")
keep <- rowSums(cpm(dge) > 1) >= 4
dge  <- dge[keep, , keep.lib.sizes = FALSE]
dge  <- calcNormFactors(dge, method = "TMM")

design   <- model.matrix(~ donor + chamber, data = meta)
v        <- voom(dge, design, plot = FALSE)
expr_mat <- v$E   # rows = gene_key, cols = sample_id

message("[4/4] Writing outputs")
gene_map <- gene_map[gene_key %chin% rownames(expr_mat)]
expr_df  <- data.table(
  gene_key    = rownames(expr_mat),
  gene_symbol = gene_map[match(rownames(expr_mat), gene_key), gene_symbol]
)
expr_df[, gene_symbol := fifelse(is.na(gene_symbol) | !nzchar(gene_symbol),
                                 gene_key, gene_symbol)]
expr_df <- cbind(expr_df, as.data.table(expr_mat))

fwrite(expr_df, file.path(out_dir, "logcpm_voom.tsv"),
       sep = "\t", na = "NA")
fwrite(meta,    file.path(out_dir, "sample_meta.tsv"),
       sep = "\t", na = "NA")

message("Wrote:")
message(" - ", file.path(out_dir, "logcpm_voom.tsv"),
        "  (", nrow(expr_df), " genes x ", nrow(meta), " samples)")
message(" - ", file.path(out_dir, "sample_meta.tsv"))
