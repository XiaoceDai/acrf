#!/usr/bin/env Rscript
# Score a single genes×samples matrix with singscore 1.28.1 for Fig 3 equivalence.

suppressPackageStartupMessages({
  library(data.table)
  library(singscore)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: validate_fig3_scorer_equivalence.R <expr.tsv> <programs.tsv> <out.tsv>")
}
expr_path <- args[[1]]
programs_path <- args[[2]]
out_path <- args[[3]]

dt <- fread(expr_path, check.names = FALSE)
genes <- as.character(dt[[1]])
mat <- as.matrix(dt[, -1, with = FALSE])
storage.mode(mat) <- "double"
rownames(mat) <- genes
stopifnot(!anyDuplicated(rownames(mat)))

programs <- fread(programs_path)
programs[, genes := lapply(
  strsplit(leading_edge_genes, ",", fixed = TRUE),
  function(x) unique(trimws(x[nzchar(trimws(x))]))
)]

ranked <- rankGenes(mat)
rows <- list()
idx <- 0L
for (i in seq_len(nrow(programs))) {
  have <- intersect(programs$genes[[i]], rownames(mat))
  if (length(have) < 3L) {
    for (sid in colnames(mat)) {
      idx <- idx + 1L
      rows[[idx]] <- data.table(
        sample_id = sid,
        program_id = programs$program_id[[i]],
        n_measured_genes = length(have),
        n_ranked_genes = nrow(mat),
        r_singscore = NA_real_
      )
    }
    next
  }
  result <- simpleScore(ranked, upSet = have, centerScore = TRUE)
  score <- result$TotalScore
  sample_ids <- rownames(result)
  for (j in seq_along(score)) {
    idx <- idx + 1L
    rows[[idx]] <- data.table(
      sample_id = sample_ids[[j]],
      program_id = programs$program_id[[i]],
      n_measured_genes = length(have),
      n_ranked_genes = nrow(mat),
      r_singscore = as.numeric(score[[j]])
    )
  }
}
fwrite(rbindlist(rows), out_path, sep = "\t")
cat(sprintf(
  "singscore %s wrote %s (%d rows)\n",
  as.character(packageVersion("singscore")),
  out_path,
  idx
))
