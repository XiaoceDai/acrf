#!/usr/bin/env Rscript
# Authoritative sample-level NCP_v3 leading-edge singscore implementation.

suppressPackageStartupMessages({
  library(data.table)
  library(singscore)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1]))
ROOT <- dirname(dirname(dirname(script_path)))
RESULTS <- file.path(ROOT, "results", "ncp_v3_singscore_reanalysis")
DATA <- file.path(ROOT, "data", "ncp_v3_singscore_reanalysis")
META <- file.path(DATA, "metadata")
PROCESSED <- file.path(DATA, "processed")

programs <- fread(file.path(RESULTS, "frozen_program_metadata.tsv"))
stopifnot(nrow(programs) == 30L, uniqueN(programs$program_id) == 30L)
programs[, le_genes := lapply(
  strsplit(leading_edge_genes, ",", fixed = TRUE),
  function(x) unique(trimws(x[nzchar(trimws(x))]))
)]

score_dataset <- function(dataset) {
  meta <- fread(file.path(META, paste0(dataset, "_sample_metadata.tsv")))
  meta <- meta[analysis_inclusion == TRUE]
  expr_dt <- fread(
    file.path(PROCESSED, dataset, "expr_symbol.tsv"),
    check.names = FALSE
  )
  genes <- as.character(expr_dt[[1]])
  expr <- as.matrix(expr_dt[, -1, with = FALSE])
  storage.mode(expr) <- "double"
  rownames(expr) <- genes
  stopifnot(!anyDuplicated(genes), all(is.finite(expr)))
  meta <- meta[match(colnames(expr), sample_id)]
  stopifnot(identical(meta$sample_id, colnames(expr)))
  ranked <- rankGenes(expr)

  rows <- vector("list", nrow(programs))
  for (i in seq_len(nrow(programs))) {
    frozen <- programs$le_genes[[i]]
    measured <- intersect(frozen, rownames(expr))
    score <- rep(NA_real_, ncol(expr))
    if (length(measured) >= 3L) {
      result <- simpleScore(ranked, upSet = measured, centerScore = TRUE)
      score <- result$TotalScore
      stopifnot(identical(rownames(result), colnames(expr)))
    }
    rows[[i]] <- data.table(
      dataset = dataset,
      sample_id = meta$sample_id,
      donor_id = meta$donor_id,
      chamber = meta$chamber,
      rhythm = meta$rhythm,
      program_id = programs$program_id[[i]],
      family_id = as.integer(programs$family_id[[i]]),
      canonical_direction = programs$canonical_direction[[i]],
      direction_sign = as.integer(programs$direction_sign[[i]]),
      n_frozen_leading_edge_genes = length(frozen),
      n_measured_genes = length(measured),
      coverage_fraction = length(measured) / length(frozen),
      program_score = score,
      score_method = paste0(
        "singscore::rankGenes(tiesMethod=min)+",
        "simpleScore(upSet=measured_frozen_LE,centerScore=TRUE)"
      ),
      gene_set_version = programs$gene_set_version[[i]]
    )
  }
  rbindlist(rows)
}

datasets <- c("GSE112339", "GSE128188", "GSE79768", "GSE115574")
scores <- rbindlist(lapply(datasets, score_dataset))
setorder(scores, dataset, sample_id, family_id, program_id)

expected_samples <- c(
  GSE112339 = 16L,
  GSE128188 = 20L,
  GSE79768 = 26L,
  GSE115574 = 56L
)
for (dataset in datasets) {
  target_dataset <- dataset
  observed <- uniqueN(scores[scores[["dataset"]] == target_dataset]$sample_id)
  if (observed != expected_samples[[dataset]]) {
    stop(dataset, ": expected ", expected_samples[[dataset]],
         " scored samples, observed ", observed)
  }
}
stopifnot(nrow(scores) == sum(expected_samples) * 30L)
stopifnot(!any(scores$n_measured_genes > scores$n_frozen_leading_edge_genes))
stopifnot(all(is.na(scores$program_score) == (scores$n_measured_genes < 3L)))

fwrite(scores, file.path(RESULTS, "program_scores_long.tsv"), sep = "\t")
coverage <- unique(scores[, .(
  dataset, program_id, family_id, canonical_direction,
  n_frozen_leading_edge_genes, n_measured_genes, coverage_fraction,
  score_is_missing = is.na(program_score)
)])
fwrite(coverage, file.path(RESULTS, "program_coverage_summary.tsv"), sep = "\t")
cat(sprintf("Wrote %d sample-program scores; missing=%d\n",
            nrow(scores), sum(is.na(scores$program_score))))
