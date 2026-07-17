#!/usr/bin/env Rscript
# Emit authoritative R singscore values for fixed fixtures and healthy data.

suppressPackageStartupMessages({
  library(data.table)
  library(singscore)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1]))
ROOT <- dirname(dirname(dirname(script_path)))
RESULTS <- file.path(ROOT, "results", "ncp_v3_singscore_reanalysis")
DATA <- file.path(ROOT, "data", "ncp_v3_singscore_reanalysis")

score_matrix <- function(dataset, matrix_path, gene_sets) {
  dt <- fread(matrix_path, check.names = FALSE)
  genes <- as.character(dt[[1]])
  mat <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- genes
  stopifnot(!anyDuplicated(rownames(mat)), all(is.finite(mat)))
  ranked <- rankGenes(mat)
  rows <- list()
  idx <- 0L
  for (i in seq_len(nrow(gene_sets))) {
    have <- intersect(gene_sets$genes[[i]], rownames(mat))
    if (length(have) < 3L) next
    result <- simpleScore(ranked, upSet = have, centerScore = TRUE)
    score <- result$TotalScore
    sample_ids <- rownames(result)
    for (j in seq_along(score)) {
      idx <- idx + 1L
      rows[[idx]] <- data.table(
        dataset = dataset,
        sample_id = sample_ids[[j]],
        program_id = gene_sets$program_id[[i]],
        n_measured_genes = length(have),
        n_ranked_genes = nrow(mat),
        r_singscore = as.numeric(score[[j]])
      )
    }
  }
  rbindlist(rows)
}

fixture_sets <- data.table(
  program_id = "fixture_ties_and_missing",
  genes = list(c("G1", "G2", "G3", "G_missing"))
)
programs <- fread(file.path(RESULTS, "frozen_program_metadata.tsv"))
programs[, genes := lapply(
  strsplit(leading_edge_genes, ",", fixed = TRUE),
  function(x) unique(trimws(x[nzchar(trimws(x))]))
)]

scores <- rbindlist(list(
  score_matrix(
    "fixture",
    file.path(RESULTS, "scorer_validation_fixture.tsv"),
    fixture_sets
  ),
  score_matrix(
    "GSE112339",
    file.path(DATA, "processed", "GSE112339", "expr_symbol.tsv"),
    programs[, .(program_id, genes)]
  )
))
fwrite(scores, file.path(RESULTS, "r_scorer_validation_scores.tsv"), sep = "\t")
fwrite(
  data.table(
    component = c("R", "singscore", "matrixStats", "GSEABase"),
    version = c(
      as.character(getRversion()),
      as.character(packageVersion("singscore")),
      as.character(packageVersion("matrixStats")),
      as.character(packageVersion("GSEABase"))
    )
  ),
  file.path(RESULTS, "singscore_software_versions.tsv"),
  sep = "\t"
)
cat(sprintf("Wrote %d authoritative R scorer validation rows\n", nrow(scores)))
