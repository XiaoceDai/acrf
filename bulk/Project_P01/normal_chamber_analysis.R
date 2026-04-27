#!/usr/bin/env Rscript
#
# Example:
# Rscript normal_chamber_analysis.R \
#   --count=GSE112339/GSE112339_HTSeq_hg19_Ensembl.txt.gz \
#   --metadata=GSE112339/metadata.xlsx \
#   --hallmark_gmt=path/to/hallmark.gmt \
#   --reactome_gmt=path/to/reactome.gmt \
#   --out_dir=.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(fgsea)
})

source("normal_chamber_qc.R", local = TRUE)
source("normal_chamber_export.R", local = TRUE)

parse_args <- function(args) {
  out <- list()
  for (a in args) {
    if (!startsWith(a, "--")) next
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[1]
    val <- if (length(kv) >= 2) paste(kv[-1], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

get_arg <- function(args, key, default = NULL) {
  if (!is.null(args[[key]])) args[[key]] else default
}

require_existing_arg <- function(path, flag, label) {
  if (is.na(path) || !nzchar(path)) {
    stop(
      label, " is required. Please pass ", flag, "=<path>.",
      call. = FALSE
    )
  }
  if (!file.exists(path)) {
    stop(
      label, " not found: ", path,
      "\nPlease verify ", flag, "=<path>.",
      call. = FALSE
    )
  }
}

normalize_chamber <- function(x) {
  x <- toupper(trimws(as.character(x)))
  ans <- rep(NA_character_, length(x))
  ans[x %chin% c("LA", "LEFT ATRIUM", "LEFT_ATRIUM", "LEFT-ATRIUM")] <- "LA"
  ans[x %chin% c("RA", "RIGHT ATRIUM", "RIGHT_ATRIUM", "RIGHT-ATRIUM")] <- "RA"
  ans
}

parse_gse112339_count_sample_id <- function(title) {
  title <- trimws(as.character(title))
  title <- sub("\\s*\\[.*$", "", title)
  parts <- strsplit(title, "_", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(NA_character_)

  donor_part <- gsub("-", "_", parts[1], fixed = TRUE)
  chamber_part <- toupper(parts[2])
  if (chamber_part == "SEPTUM") chamber_part <- "IVS"
  paste(chamber_part, donor_part, sep = "_")
}

extract_geo_heart_id <- function(x) {
  x <- as.character(x)
  out <- sub(".*heart id:\\s*", "", x, ignore.case = TRUE)
  out[!grepl("heart id\\s*:", x, ignore.case = TRUE)] <- NA_character_
  out
}

load_gse112339_geo_map <- function() {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop(
      "GEOquery is required to map GSM accessions to GSE112339 count-matrix sample IDs.",
      call. = FALSE
    )
  }

  message("Mapping metadata GSM accessions through GEO GSE112339 annotation")
  gse <- GEOquery::getGEO("GSE112339", GSEMatrix = TRUE, getGPL = FALSE)
  if (length(gse) < 1) stop("Failed to load GEO annotation for GSE112339.", call. = FALSE)

  pd <- as.data.table(Biobase::pData(gse[[1]]))
  if (!all(c("geo_accession", "title", "characteristics_ch1.2") %in% names(pd))) {
    stop("GEO annotation is missing expected columns for GSE112339.", call. = FALSE)
  }

  geo_map <- pd[, .(
    sample_id = as.character(geo_accession),
    count_sample_id = vapply(title, parse_gse112339_count_sample_id, character(1)),
    donor = extract_geo_heart_id(`characteristics_ch1.2`),
    chamber = normalize_chamber(sub(".*_", "", sub("\\s*\\[.*$", "", title)))
  )]
  geo_map <- geo_map[
    !is.na(sample_id) & nzchar(sample_id) &
      !is.na(count_sample_id) & nzchar(count_sample_id)
  ]
  unique(geo_map)
}

read_external_metadata <- function(metadata_path) {
  ext <- tolower(tools::file_ext(metadata_path))
  ext_meta <- switch(
    ext,
    "tsv" = fread(metadata_path),
    "txt" = fread(metadata_path),
    "csv" = fread(metadata_path),
    "xlsx" = {
      if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("metadata is .xlsx but readxl is unavailable.", call. = FALSE)
      }
      as.data.table(readxl::read_xlsx(metadata_path))
    },
    stop("Unsupported metadata extension: ", ext)
  )

  clean_names <- tolower(gsub("[^a-zA-Z0-9]+", "_", names(ext_meta)))
  setnames(ext_meta, names(ext_meta), clean_names)

  sample_col <- if ("sample_id" %in% names(ext_meta)) "sample_id" else if ("geo_accession" %in% names(ext_meta)) "geo_accession" else NA_character_
  donor_col <- if ("donor_id" %in% names(ext_meta)) "donor_id" else if ("donor" %in% names(ext_meta)) "donor" else NA_character_
  chamber_col <- if ("chamber" %in% names(ext_meta)) "chamber" else NA_character_

  if (is.na(sample_col) || is.na(donor_col) || is.na(chamber_col)) {
    stop("Metadata must contain sample_id (or geo_accession), donor_id (or donor), and chamber.", call. = FALSE)
  }

  ext_meta[, .(
    sample_id = as.character(get(sample_col)),
    donor = as.character(get(donor_col)),
    chamber = normalize_chamber(get(chamber_col))
  )]
}

harmonize_metadata_sample_ids <- function(ext_keep, sample_ids) {
  overlap_n <- sum(ext_keep$sample_id %chin% sample_ids)
  has_gsm_ids <- any(grepl("^GSM\\d+$", ext_keep$sample_id, ignore.case = TRUE))
  if (!has_gsm_ids || overlap_n == nrow(ext_keep)) {
    return(ext_keep)
  }

  geo_map <- load_gse112339_geo_map()
  ext_keep <- merge(ext_keep, geo_map, by = "sample_id", all.x = TRUE, sort = FALSE, suffixes = c("_meta", "_geo"))
  ext_keep[, sample_id := fifelse(!is.na(count_sample_id) & count_sample_id %chin% sample_ids, count_sample_id, sample_id)]
  ext_keep[is.na(donor_meta) | !nzchar(donor_meta), donor_meta := donor_geo]
  ext_keep[is.na(chamber_meta) | !nzchar(chamber_meta), chamber_meta := chamber_geo]
  unique(ext_keep[, .(sample_id, donor = donor_meta, chamber = chamber_meta)])
}

read_metadata <- function(metadata_path, sample_ids) {
  meta <- data.table(sample_id = sample_ids)
  meta[, chamber := sub("_.*$", "", sample_id)]
  meta[, donor := sub("^[^_]+_", "", sample_id)]

  if (is.na(metadata_path) || !nzchar(metadata_path) || !file.exists(metadata_path)) {
    return(meta)
  }

  ext_keep <- read_external_metadata(metadata_path)
  ext_keep <- ext_keep[!is.na(sample_id) & nzchar(sample_id)]
  ext_keep <- ext_keep[!duplicated(sample_id)]
  ext_keep <- harmonize_metadata_sample_ids(ext_keep, sample_ids)
  ext_keep <- ext_keep[!duplicated(sample_id)]

  meta <- merge(meta, ext_keep, by = "sample_id", all.x = TRUE, suffixes = c("_infer", "_meta"))
  meta[, donor := fifelse(!is.na(donor_meta) & nzchar(donor_meta), donor_meta, donor_infer)]
  meta[, chamber := fifelse(!is.na(chamber_meta) & nzchar(chamber_meta), chamber_meta, chamber_infer)]
  meta <- meta[, .(sample_id, donor, chamber)]
  meta
}

fit_chamber <- function(dge, meta) {
  design <- model.matrix(~ donor + chamber, data = meta)
  coef_name <- grep("^chamber", colnames(design), value = TRUE)
  if (length(coef_name) != 1) {
    stop("Expected one chamber coefficient, got: ", paste(coef_name, collapse = ", "))
  }
  v <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  tt <- as.data.table(topTable(fit, coef = coef_name, number = Inf, sort.by = "none"), keep.rownames = "gene_key")
  list(v = v, fit = fit, tt = tt, coef_name = coef_name)
}

make_rank_stats <- function(gene_table) {
  rank_dt <- gene_table[, .(
    rank_id = fifelse(!is.na(gene_symbol) & nzchar(gene_symbol), gene_symbol, ensembl_id),
    stat = moderated_t
  )]
  rank_dt <- rank_dt[!is.na(rank_id) & nzchar(rank_id) & is.finite(stat)]
  rank_dt <- rank_dt[order(-abs(stat))]
  rank_dt <- rank_dt[!duplicated(rank_id)]
  stats <- rank_dt$stat
  names(stats) <- rank_dt$rank_id
  sort(stats, decreasing = TRUE)
}

load_pathways <- function(hallmark_gmt, reactome_gmt) {
  if (!file.exists(hallmark_gmt)) stop("Hallmark GMT not found: ", hallmark_gmt)
  if (!file.exists(reactome_gmt)) stop("Reactome GMT not found: ", reactome_gmt)

  h <- fgsea::gmtPathways(hallmark_gmt)
  r <- fgsea::gmtPathways(reactome_gmt)
  names(h) <- paste0("HALLMARK__", names(h))
  names(r) <- paste0("REACTOME__", names(r))
  c(h, r)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  count_path <- get_arg(args, "count", "GSE112339/GSE112339_HTSeq_hg19_Ensembl.txt.gz")
  metadata_path <- get_arg(args, "metadata", "GSE112339/metadata.xlsx")
  hallmark_gmt <- get_arg(args, "hallmark_gmt", "")
  reactome_gmt <- get_arg(args, "reactome_gmt", "")
  out_dir <- get_arg(args, "out_dir", ".")
  top_programs <- as.integer(get_arg(args, "top_programs", "5"))
  gene_fdr_cutoff <- as.numeric(get_arg(args, "gene_fdr", "0.05"))
  path_fdr_cutoff <- as.numeric(get_arg(args, "pathway_fdr", "0.05"))
  top_heatmap_genes_per_side <- as.integer(get_arg(args, "top_heatmap_genes_per_side", "20"))

  require_existing_arg(count_path, "--count", "Count matrix")
  require_existing_arg(metadata_path, "--metadata", "Metadata file")
  require_existing_arg(hallmark_gmt, "--hallmark_gmt", "Hallmark GMT")
  require_existing_arg(reactome_gmt, "--reactome_gmt", "Reactome GMT")

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  qc_dir <- file.path(out_dir, "qc")
  if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)
  fig_dir <- file.path(out_dir, "figures")
  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

  message("[1/8] Reading count matrix")
  counts_dt <- fread(count_path)
  gene_id_col <- names(counts_dt)[1]
  sample_ids <- names(counts_dt)[-1]

  message("[2/8] Building metadata and enforcing paired LA/RA")
  meta <- read_metadata(metadata_path, sample_ids)
  meta <- meta[chamber %chin% c("LA", "RA")]
  meta <- meta[sample_id %chin% sample_ids]
  setorder(meta, donor, chamber)

  pair_tab <- dcast(meta, donor ~ chamber, value.var = "sample_id", fun.aggregate = length, fill = 0)
  valid_donor <- pair_tab[LA == 1 & RA == 1, donor]
  meta <- meta[donor %chin% valid_donor]

  if (length(unique(meta$donor)) != 8) {
    stop("Metadata curation failed: need exactly 8 paired donors; got ", length(unique(meta$donor)))
  }
  if (nrow(meta) != 16) {
    stop("Metadata curation failed: need exactly 16 LA/RA samples; got ", nrow(meta))
  }

  meta[, donor := factor(donor)]
  meta[, chamber := factor(chamber, levels = c("RA", "LA"))]

  sample_use <- meta$sample_id
  miss <- setdiff(sample_use, names(counts_dt))
  if (length(miss) > 0) stop("Samples missing in count matrix: ", paste(miss, collapse = ", "))

  raw_gene_id <- as.character(counts_dt[[gene_id_col]])
  gene_key <- make.unique(raw_gene_id, sep = "_dup")
  gene_map <- data.table(
    gene_key = gene_key,
    gene_id_raw = raw_gene_id,
    ensembl_id = ifelse(grepl("^ENSG", raw_gene_id, ignore.case = TRUE), sub("\\..*$", "", raw_gene_id), raw_gene_id),
    gene_symbol = ifelse(grepl("^ENSG", raw_gene_id, ignore.case = TRUE), NA_character_, raw_gene_id)
  )

  counts_mat <- as.matrix(counts_dt[, ..sample_use])
  rownames(counts_mat) <- gene_key
  storage.mode(counts_mat) <- "integer"

  dge <- DGEList(counts = counts_mat)

  message("[3/8] Minimal QC")
  qc_res <- run_qc_outputs(dge, meta, qc_dir, fig_dir)
  log_cpm_all <- qc_res$log_cpm_all

  message("[4/8] Low-expression filtering and normalization")
  keep <- rowSums(cpm(dge) > 1) >= 4
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  gene_map <- gene_map[gene_key %chin% rownames(dge)]

  message("[5/8] Paired chamber model and full differential ranking")
  fit_full <- fit_chamber(dge, meta)
  gene_out <- merge(fit_full$tt, gene_map, by = "gene_key", all.x = TRUE, sort = FALSE)
  gene_out[, `:=`(
    logFC = logFC,
    moderated_t = t,
    p_value = P.Value,
    fdr = adj.P.Val,
    average_expression = AveExpr,
    direction = ifelse(logFC >= 0, "LA-high", "RA-high")
  )]
  gene_out[, fgsea_id := fifelse(!is.na(gene_symbol) & nzchar(gene_symbol), gene_symbol, ensembl_id)]

  message("[6/8] Leave-one-donor-out stability (genes + programs)")
  donors <- levels(meta$donor)
  lodo_gene <- data.table(gene_key = gene_out$gene_key)
  lodo_prog_list <- vector("list", length(donors))
  names(lodo_prog_list) <- donors

  pathways <- load_pathways(hallmark_gmt, reactome_gmt)
  rank_full <- make_rank_stats(gene_out)
  fgsea_full <- as.data.table(fgsea(pathways = pathways, stats = rank_full, minSize = 15, maxSize = 500))
  fgsea_full[, direction := ifelse(NES >= 0, "LA-biased", "RA-biased")]
  fgsea_full[, leading_edge_genes := vapply(leadingEdge, function(x) paste(x, collapse = ","), character(1))]
  fgsea_full[, leading_edge_n := lengths(leadingEdge)]

  for (d in donors) {
    idx <- meta$donor != d
    meta_i <- droplevels(meta[idx])
    dge_i <- dge[, idx, keep.lib.sizes = FALSE]
    dge_i <- calcNormFactors(dge_i, method = "TMM")
    fit_i <- fit_chamber(dge_i, meta_i)

    t_i <- fit_i$tt[, .(gene_key, t)]
    setnames(t_i, "t", paste0("t_drop_", d))
    lodo_gene <- merge(lodo_gene, t_i, by = "gene_key", all.x = TRUE, sort = FALSE)

    gene_i <- merge(fit_i$tt, gene_map, by = "gene_key", all.x = TRUE, sort = FALSE)
    gene_i[, moderated_t := t]
    rank_i <- make_rank_stats(gene_i)
    fg_i <- as.data.table(fgsea(pathways = pathways, stats = rank_i, minSize = 15, maxSize = 500))
    lodo_prog_list[[d]] <- fg_i[, .(pathway, NES)]
  }

  drop_gene_cols <- grep("^t_drop_", names(lodo_gene), value = TRUE)
  gene_stab <- merge(
    gene_out[, .(gene_key, moderated_t)],
    lodo_gene,
    by = "gene_key",
    all.x = TRUE,
    sort = FALSE
  )
  gene_drop_mat <- as.matrix(gene_stab[, ..drop_gene_cols])
  full_sign <- sign(gene_stab$moderated_t)
  same_sign <- sweep(sign(gene_drop_mat), 1, full_sign, `==`)
  non_zero <- sign(gene_drop_mat) != 0
  gene_out[, stability_consistent_runs := rowSums(same_sign & non_zero, na.rm = TRUE)]
  gene_out[, stable_gene := stability_consistent_runs >= 7]

  prog_stab <- fgsea_full[, .(pathway, NES)]
  for (d in donors) {
    tmp <- copy(lodo_prog_list[[d]])
    setnames(tmp, "NES", paste0("NES_drop_", d))
    prog_stab <- merge(prog_stab, tmp, by = "pathway", all.x = TRUE, sort = FALSE)
  }

  drop_prog_cols <- grep("^NES_drop_", names(prog_stab), value = TRUE)
  prog_drop_mat <- as.matrix(prog_stab[, ..drop_prog_cols])
  prog_sign <- sign(prog_stab$NES)
  prog_same <- sweep(sign(prog_drop_mat), 1, prog_sign, `==`)
  prog_nonzero <- sign(prog_drop_mat) != 0
  prog_stab[, stability_consistent_runs := rowSums(prog_same & prog_nonzero, na.rm = TRUE)]
  prog_stab[, stable_program := stability_consistent_runs >= 6]

  fgsea_full <- merge(
    fgsea_full,
    prog_stab[, .(pathway, stability_consistent_runs, stable_program)],
    by = "pathway",
    all.x = TRUE,
    sort = FALSE
  )

  message("[7/8] Core programs and core genes")
  sig_genes <- unique(gene_out[fdr < gene_fdr_cutoff, fgsea_id])
  sig_prog <- fgsea_full[padj < path_fdr_cutoff]

  setorder(sig_prog, -NES)
  la_keep <- head(sig_prog[direction == "LA-biased"], top_programs)
  setorder(sig_prog, NES)
  ra_keep <- head(sig_prog[direction == "RA-biased"], top_programs)
  prog_out <- rbindlist(list(la_keep, ra_keep), fill = TRUE)
  if (nrow(prog_out) == 0) {
    prog_out <- sig_prog[0]
  }

  core_list <- lapply(prog_out$pathway, function(pw) {
    le <- fgsea_full[pathway == pw, leadingEdge][[1]]
    core <- intersect(le, sig_genes)
    data.table(
      pathway = pw,
      core_genes = paste(core, collapse = ","),
      core_gene_n = length(core)
    )
  })
  core_dt <- if (length(core_list) > 0) rbindlist(core_list) else data.table(pathway = character(), core_genes = character(), core_gene_n = integer())
  prog_out <- merge(prog_out, core_dt, by = "pathway", all.x = TRUE, sort = FALSE)
  prog_out[is.na(core_genes), core_genes := ""]
  prog_out[is.na(core_gene_n), core_gene_n := 0L]

  core_map <- if (nrow(prog_out) > 0) {
    rbindlist(lapply(seq_len(nrow(prog_out)), function(i) {
      pw <- prog_out$pathway[i]
      cgs <- unlist(strsplit(prog_out$core_genes[i], ",", fixed = TRUE))
      cgs <- cgs[nzchar(cgs)]
      if (length(cgs) == 0) return(data.table(fgsea_id = character(), pathway = character()))
      data.table(fgsea_id = cgs, pathway = pw)
    }), fill = TRUE)
  } else {
    data.table(fgsea_id = character(), pathway = character())
  }

  if (nrow(core_map) > 0) {
    core_gene_tag <- core_map[, .(core_programs = paste(unique(pathway), collapse = ";")), by = fgsea_id]
    gene_out <- merge(gene_out, core_gene_tag, by = "fgsea_id", all.x = TRUE, sort = FALSE)
  } else {
    gene_out[, core_programs := NA_character_]
  }
  gene_out[, is_core_gene := !is.na(core_programs) & nzchar(core_programs)]

  message("[8/8] Writing final outputs")
  write_final_outputs(
    gene_out = gene_out,
    prog_out = prog_out,
    expr_mat = fit_full$v$E,
    meta = meta,
    out_dir = out_dir,
    fig_dir = fig_dir,
    gene_fdr_cutoff = gene_fdr_cutoff,
    top_heatmap_genes_per_side = top_heatmap_genes_per_side
  )
}

main()
