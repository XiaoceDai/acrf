#!/usr/bin/env Rscript
# Re-export selection-rank fgsea leading edges for the frozen 30 programs.

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1]))
ROOT <- dirname(dirname(dirname(script_path)))
RESULTS <- file.path(ROOT, "results", "ncp_v3_singscore_reanalysis")

DE_PATH <- file.path(ROOT, "bulk", "Project_P01", "bulk_chamber_DE_full.tsv")
GMT_HALL <- file.path(ROOT, "bulk", "Project_P01", "gene_sets",
                      "h.all.v2026.1.Hs.symbols.gmt")
GMT_REAC <- file.path(ROOT, "bulk", "Project_P01", "gene_sets",
                      "c2.cp.reactome.v2026.1.Hs.symbols.gmt")
FROZEN_PATH <- file.path(ROOT, "results", "fig3_v3_program_le.tsv")

de <- fread(DE_PATH)
de[, gene := gsub('"', '', gene_symbol)]
de <- de[!is.na(t_stat) & !is.na(gene) & nzchar(gene)]
ranked <- setNames(de$t_stat, de$gene)
ranked <- ranked[!duplicated(names(ranked))]
ranked <- sort(ranked, decreasing = TRUE)

all_pathways <- c(
  setNames(fgsea::gmtPathways(GMT_HALL),
           paste0("HALLMARK__", names(fgsea::gmtPathways(GMT_HALL)))),
  setNames(fgsea::gmtPathways(GMT_REAC),
           paste0("REACTOME__", names(fgsea::gmtPathways(GMT_REAC))))
)
frozen <- fread(FROZEN_PATH)
stopifnot(nrow(frozen) == 30L, !anyDuplicated(frozen$program_id))
missing_pathways <- setdiff(frozen$program_id, names(all_pathways))
if (length(missing_pathways)) {
  stop("Frozen programs missing from GMT: ", paste(missing_pathways, collapse = ", "))
}

set.seed(42L)
selection_fgsea <- as.data.table(
  fgsea(
    pathways = all_pathways[frozen$program_id],
    stats = ranked,
    minSize = 15L,
    maxSize = 500L,
    nproc = 1
  )
)
selection_fgsea[, selection_leading_edge_genes :=
                  vapply(leadingEdge, paste, collapse = ",", FUN.VALUE = character(1))]
selection_fgsea[, n_selection_le := lengths(leadingEdge)]
selection_fgsea[, leadingEdge := NULL]

split_genes <- function(x) {
  genes <- trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1]])
  sort(unique(genes[nzchar(genes)]))
}

comparison <- merge(
  frozen[, .(program_id, program_order, direction, n_frozen_le = n_LE_total,
             frozen_leading_edge_genes = leading_edge_genes)],
  selection_fgsea[, .(program_id = pathway, NES, padj, n_selection_le,
                      selection_leading_edge_genes)],
  by = "program_id",
  all.x = TRUE,
  sort = FALSE
)
comparison[, `:=`(
  n_intersection = mapply(
    function(a, b) length(intersect(split_genes(a), split_genes(b))),
    frozen_leading_edge_genes, selection_leading_edge_genes
  ),
  n_union = mapply(
    function(a, b) length(union(split_genes(a), split_genes(b))),
    frozen_leading_edge_genes, selection_leading_edge_genes
  ),
  set_identical = mapply(
    function(a, b) identical(split_genes(a), split_genes(b)),
    frozen_leading_edge_genes, selection_leading_edge_genes
  )
)]
comparison[, jaccard := n_intersection / n_union]
setorder(comparison, program_order)

dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
fwrite(
  selection_fgsea[, .(program_id = pathway, NES, padj, n_selection_le,
                      selection_leading_edge_genes)],
  file.path(RESULTS, "selection_time_leading_edges.tsv"),
  sep = "\t"
)
fwrite(comparison, file.path(RESULTS, "leading_edge_comparison.tsv"), sep = "\t")

summary <- data.table(
  n_programs = nrow(comparison),
  n_set_identical = sum(comparison$set_identical),
  min_jaccard = min(comparison$jaccard),
  median_jaccard = median(comparison$jaccard),
  gate_passed = all(comparison$set_identical)
)
fwrite(summary, file.path(RESULTS, "leading_edge_comparison_summary.tsv"), sep = "\t")
cat(sprintf("Leading-edge set identity: %d/%d; min Jaccard %.6f\n",
            summary$n_set_identical, summary$n_programs, summary$min_jaccard))
if (!summary$gate_passed) {
  quit(status = 42L)
}
