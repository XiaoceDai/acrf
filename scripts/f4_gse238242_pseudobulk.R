#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# F4 - Donor x compartment pseudobulk for GSE238242 (LAA snRNA-seq).
#
# Frozen rules (Figure 5 spec, locked by user):
#   unit         = donor x compartment
#   aggregation  = SUM of RAW counts across cells in (donor, compartment)
#   normalization / log / singscore -> done in F5, NOT here
#
# Cell-type -> compartment mapping (locked):
#   MAIN
#     CM           -> atrialCM
#     FB           -> fibroblast
#     EC           -> endothelial_or_vascular
#     Endocardial  -> endothelial_or_vascular
#     Myeloid      -> immune
#     Lymphoid     -> immune
#     Mast         -> immune
#   SUPP (mural, separate output)
#     PC           -> mural
#     SMC          -> mural
#   EXCLUDED (not emitted)
#     Mesothelial / Neuron / Adipo
#
# Outputs (Project_P5/GSE238242/pseudobulk/):
#   gse238242_pseudobulk_counts_main.rds
#   gse238242_pseudobulk_counts_main.tsv.gz
#   gse238242_pseudobulk_meta_main.tsv
#   gse238242_pseudobulk_summary_main.tsv
#   gse238242_pseudobulk_counts_supp_mural.rds
#   gse238242_pseudobulk_counts_supp_mural.tsv.gz
#   gse238242_pseudobulk_meta_supp_mural.tsv
#   gse238242_pseudobulk_summary_supp_mural.tsv
#   qc_pseudobulk.tsv      (gate failures, MIN_CELLS = 30)
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
GDIR    <- file.path(ROOT, "Project_P5", "GSE238242")
EXT     <- file.path(GDIR, "extracted")
OUT_DIR <- file.path(GDIR, "pseudobulk")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MIN_CELLS <- 30L

CT_MAIN <- list(
  CM           = "atrialCM",
  FB           = "fibroblast",
  EC           = "endothelial_or_vascular",
  Endocardial  = "endothelial_or_vascular",
  Myeloid      = "immune",
  Lymphoid     = "immune",
  Mast         = "immune"
)
CT_SUPP <- list(PC = "mural", SMC = "mural")
CT_EXCLUDED <- c("Mesothelial", "Neuron", "Adipo")

cat(sprintf("ROOT = %s\n", ROOT))

# --------------------------------------------------------------------
# 1. Metadata
# --------------------------------------------------------------------
meta <- fread(
  file.path(GDIR, "GSE238242_snAF.metadata.tsv.gz"),
  showProgress = FALSE, header = FALSE, skip = 1L
)
setnames(meta, c("barcode", "cell_type", "sample", "sex", "Rhythm"))
setnames(meta, "Rhythm", "rhythm")
cat(sprintf("metadata: %d nuclei, %d donors\n",
            nrow(meta), uniqueN(meta$sample)))

meta[, compartment_main := unlist(CT_MAIN)[match(cell_type, names(CT_MAIN))]]
meta[, compartment_supp := unlist(CT_SUPP)[match(cell_type, names(CT_SUPP))]]
meta[, in_main := !is.na(compartment_main)]
meta[, in_supp := !is.na(compartment_supp)]
meta[, in_excl := cell_type %in% CT_EXCLUDED]

cat("\n--- cell_type breakdown ---\n")
print(meta[, .(.N,
               main = sum(in_main),
               supp = sum(in_supp),
               excl = sum(in_excl)), by = cell_type][order(-N)])

# --------------------------------------------------------------------
# 2. Per-sample RNA matrices (genes x cells, raw counts)
# --------------------------------------------------------------------
rna_files <- list.files(EXT, pattern = "RNA\\.counts\\.tsv\\.gz$", full.names = TRUE)
donor_of <- function(p) regmatches(basename(p),
                                   regexpr("CF[0-9]+", basename(p)))
rna_donor <- vapply(rna_files, donor_of, character(1))
stopifnot(uniqueN(rna_donor) == length(rna_donor))
cat(sprintf("\nRNA files (donors): %s\n", paste(rna_donor, collapse = ", ")))

# Container: collect column-aggregates as we read each donor file.
# Output two matrices later: MAIN (gene x pb_main) and SUPP (gene x pb_supp).
all_genes <- NULL
acc_main_cols <- list()  # list of named numeric vectors (length = n_genes)
acc_supp_cols <- list()
meta_main_rows <- list()
meta_supp_rows <- list()

for (k in seq_along(rna_files)) {
  fp <- rna_files[k]
  donor <- rna_donor[k]
  cat(sprintf("\n[%d/%d] %s (donor=%s) ... ",
              k, length(rna_files), basename(fp), donor))
  t0 <- Sys.time()

  # File format: header line = cell barcodes (no entry for gene column);
  # data lines = gene_symbol \t counts ...  -> manual header read.
  hdr <- scan(fp, what = character(), sep = "\t", nlines = 1L,
              quiet = TRUE, quote = "\"")
  m <- fread(fp, sep = "\t", header = FALSE, skip = 1L,
             showProgress = FALSE, col.names = c("gene", hdr),
             check.names = FALSE)
  cells <- hdr

  if (is.null(all_genes)) {
    all_genes <- m$gene
  } else {
    if (!identical(m$gene, all_genes)) {
      # Re-order to match
      stopifnot(setequal(m$gene, all_genes))
      m <- m[match(all_genes, gene)]
    }
  }

  # Align cells with metadata for this donor
  md_d <- meta[sample == donor]
  hit <- md_d$barcode %in% cells
  md_d <- md_d[hit]
  cell_in_file <- intersect(cells, md_d$barcode)
  cat(sprintf("cells: file=%d  metadata=%d  intersect=%d  ",
              length(cells), nrow(md_d), length(cell_in_file)))
  if (length(cell_in_file) == 0L) {
    cat("[skip — no overlap]\n"); next
  }

  # Subset and convert to sparse
  M <- as.matrix(m[, ..cell_in_file])
  storage.mode(M) <- "double"
  rownames(M) <- m$gene
  colnames(M) <- cell_in_file

  # Make a barcode -> compartment lookup ordered like cells in M
  md_lookup <- md_d[match(cell_in_file, barcode)]
  rhythm_v <- md_lookup$rhythm[1]
  sex_v <- md_lookup$sex[1]
  stopifnot(all(md_lookup$rhythm == rhythm_v),
            all(md_lookup$sex == sex_v))

  # MAIN compartments
  for (cmp in unique(na.omit(md_lookup$compartment_main))) {
    cidx <- which(md_lookup$compartment_main == cmp)
    n_c <- length(cidx)
    pb <- rowSums(M[, cidx, drop = FALSE])
    pb_id <- paste(donor, rhythm_v, cmp, sep = "__")
    acc_main_cols[[pb_id]] <- pb
    meta_main_rows[[pb_id]] <- data.table(
      pb_id = pb_id, donor = donor, rhythm = rhythm_v,
      compartment = cmp, n_cells = n_c, sex = sex_v,
      source_sample = donor, source = "GSE238242"
    )
  }

  # SUPP compartments (mural)
  for (cmp in unique(na.omit(md_lookup$compartment_supp))) {
    cidx <- which(md_lookup$compartment_supp == cmp)
    n_c <- length(cidx)
    pb <- rowSums(M[, cidx, drop = FALSE])
    pb_id <- paste(donor, rhythm_v, paste0("supp_", cmp), sep = "__")
    acc_supp_cols[[pb_id]] <- pb
    meta_supp_rows[[pb_id]] <- data.table(
      pb_id = pb_id, donor = donor, rhythm = rhythm_v,
      compartment = paste0("supp_", cmp), n_cells = n_c, sex = sex_v,
      source_sample = donor, source = "GSE238242"
    )
  }
  rm(M); gc(verbose = FALSE)
  cat(sprintf("done (%.1fs)\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# --------------------------------------------------------------------
# 3. Assemble matrices and apply MIN_CELLS gate
# --------------------------------------------------------------------
build_pb <- function(cols_list, meta_list) {
  if (length(cols_list) == 0)
    return(list(M = NULL, meta = data.table(), gate = data.table()))
  meta_dt <- rbindlist(meta_list, use.names = TRUE, fill = TRUE)
  setorder(meta_dt, compartment, rhythm, donor)
  M <- do.call(cbind, cols_list[meta_dt$pb_id])
  rownames(M) <- all_genes
  colnames(M) <- meta_dt$pb_id

  # gate
  gate_fail <- meta_dt[n_cells < MIN_CELLS]
  if (nrow(gate_fail)) {
    keep <- !(meta_dt$pb_id %in% gate_fail$pb_id)
    M <- M[, keep, drop = FALSE]
    meta_dt <- meta_dt[keep]
  }
  list(M = M, meta = meta_dt, gate = gate_fail)
}

cat("\nAssembling MAIN ...\n")
main <- build_pb(acc_main_cols, meta_main_rows)
cat("Assembling SUPP (mural) ...\n")
supp <- build_pb(acc_supp_cols, meta_supp_rows)

write_one <- function(tag, lst) {
  if (is.null(lst$M)) {
    cat(sprintf("  %s: empty, skipping\n", tag)); return()
  }
  M <- lst$M
  meta_out <- lst$meta

  # Sum-of-counts pseudobulks; collapse duplicate gene symbols by SUM
  if (anyDuplicated(rownames(M))) {
    cat(sprintf("  %s: collapsing %d duplicate symbols by SUM\n",
                tag, sum(duplicated(rownames(M)))))
    sym <- rownames(M)
    rs <- rowsum(M, group = sym, reorder = FALSE)
    M <- as.matrix(rs)
  }

  rds <- file.path(OUT_DIR,
                   sprintf("gse238242_pseudobulk_counts_%s.rds", tag))
  tsv <- file.path(OUT_DIR,
                   sprintf("gse238242_pseudobulk_counts_%s.tsv.gz", tag))
  meta_p <- file.path(OUT_DIR,
                      sprintf("gse238242_pseudobulk_meta_%s.tsv", tag))
  sum_p  <- file.path(OUT_DIR,
                      sprintf("gse238242_pseudobulk_summary_%s.tsv", tag))

  M_sp <- as(as(M, "dgeMatrix"), "CsparseMatrix")
  saveRDS(M_sp, rds)

  # tsv.gz
  out <- data.table(gene = rownames(M))
  out <- cbind(out, as.data.table(round(M, 6)))
  fwrite(out, tsv, sep = "\t", compress = "gzip")

  fwrite(meta_out, meta_p, sep = "\t")

  summ <- meta_out[, .(
    n_pb       = .N,
    n_donors   = uniqueN(donor),
    n_AF       = uniqueN(donor[rhythm == "AF"]),
    n_SR       = uniqueN(donor[rhythm == "SR"]),
    median_cells = as.integer(median(n_cells)),
    min_cells  = min(n_cells),
    max_cells  = max(n_cells),
    sum_cells  = sum(n_cells)
  ), by = compartment][order(compartment)]
  fwrite(summ, sum_p, sep = "\t")

  cat(sprintf("\n  wrote %s  (%d pb x %d genes)\n",
              sub(paste0("^", ROOT, "/"), "", tsv),
              ncol(M), nrow(M)))
  cat("  summary:\n")
  print(as.data.frame(summ), row.names = FALSE)
}

write_one("main",       main)
write_one("supp_mural", supp)

# QC: gate failures
gate_all <- rbindlist(list(main$gate, supp$gate), use.names = TRUE, fill = TRUE)
fwrite(gate_all, file.path(OUT_DIR, "qc_pseudobulk.tsv"), sep = "\t")
cat(sprintf("\nQC: %d (donor x compartment) below MIN_CELLS=%d (see qc_pseudobulk.tsv)\n",
            nrow(gate_all), MIN_CELLS))
cat("\nF4 DONE.\n")
