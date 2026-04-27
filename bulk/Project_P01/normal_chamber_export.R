clean_program_name <- function(x) {
  x <- sub("^HALLMARK__", "", x)
  x <- sub("^REACTOME__", "", x)
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^REACTOME_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

classify_program <- function(program_name) {
  x <- tolower(program_name)
  if (grepl("oxidative phosphorylation|respiratory electron transport|aerobic respiration|mitochondrial", x)) return("mitochondrial metabolism")
  if (grepl("translation", x)) return("translation")
  if (grepl("heparan sulfate|gag|glycosaminoglycan|ecm|matrix|collagen", x)) return("ECM/GAG remodeling")
  if (grepl("bmp|runx|development|signaling", x)) return("developmental signaling")
  "other"
}

top_feature_genes <- function(x, n = 8L) {
  genes <- unlist(strsplit(as.character(x), ",", fixed = TRUE))
  genes <- genes[nzchar(genes)]
  paste(head(genes, n), collapse = ",")
}

write_final_outputs <- function(gene_out, prog_out, expr_mat, meta, out_dir, fig_dir, gene_fdr_cutoff, top_heatmap_genes_per_side) {
  setorder(gene_out, -moderated_t)
  gene_final <- gene_out[, .(
    ensembl_id,
    gene_symbol,
    logFC = as.numeric(logFC),
    moderated_t = as.numeric(moderated_t),
    p_value = as.numeric(p_value),
    FDR = as.numeric(fdr),
    average_expression = as.numeric(average_expression),
    direction,
    stability_consistent_runs = as.integer(stability_consistent_runs),
    stable_gene = as.logical(stable_gene),
    is_core_gene = as.logical(is_core_gene),
    core_programs
  )]

  gene_core_final <- gene_final[
    stable_gene == TRUE & FDR < gene_fdr_cutoff & abs(logFC) >= 0.5,
    .(
      gene = fifelse(!is.na(gene_symbol) & nzchar(gene_symbol), gene_symbol, ensembl_id),
      gene_id = ensembl_id,
      logFC,
      FDR,
      direction,
      stability = as.integer(stability_consistent_runs),
      core_programs
    )
  ]
  gene_core_final[, abs_logFC := abs(logFC)]
  setorder(gene_core_final, direction, -abs_logFC, FDR)
  gene_core_final[, abs_logFC := NULL]

  gene_core_strict_final <- gene_final[
    stable_gene == TRUE & FDR < gene_fdr_cutoff & abs(logFC) >= 0.5 & is_core_gene == TRUE,
    .(
      gene = fifelse(!is.na(gene_symbol) & nzchar(gene_symbol), gene_symbol, ensembl_id),
      gene_id = ensembl_id,
      logFC,
      FDR,
      direction,
      stability = as.integer(stability_consistent_runs),
      core_programs
    )
  ]
  gene_core_strict_final[, abs_logFC := abs(logFC)]
  setorder(gene_core_strict_final, direction, -abs_logFC, FDR)
  gene_core_strict_final[, abs_logFC := NULL]

  prog_final <- prog_out[, .(
    program_id = pathway,
    program = vapply(pathway, clean_program_name, character(1)),
    program_class = vapply(vapply(pathway, clean_program_name, character(1)), classify_program, character(1)),
    direction,
    NES = as.numeric(NES),
    FDR = as.numeric(padj),
    top_core_genes = vapply(core_genes, top_feature_genes, character(1)),
    leading_edge_genes,
    core_genes,
    leading_edge_n = as.integer(leading_edge_n),
    core_gene_n = as.integer(core_gene_n),
    stability_consistent_runs = as.integer(stability_consistent_runs),
    stable_program = as.logical(stable_program)
  )]

  bulk_de_full <- gene_final[, .(gene_id = ensembl_id, gene_symbol, logFC, FDR, t_stat = moderated_t)]
  bulk_stability <- gene_final[, .(gene_id = ensembl_id, gene_symbol, direction, stability_score = as.integer(stability_consistent_runs))]
  bulk_core_genes <- gene_final[
    stable_gene == TRUE & FDR < gene_fdr_cutoff & abs(logFC) > 0.5 & stability_consistent_runs >= 7,
    .(gene_id = ensembl_id, gene_symbol, direction, logFC, FDR, stability_score = as.integer(stability_consistent_runs))
  ]
  bulk_core_genes[, abs_logFC := abs(logFC)]
  setorder(bulk_core_genes, direction, -abs_logFC, FDR)
  bulk_core_genes[, abs_logFC := NULL]

  bulk_fgsea <- prog_final[, .(
    pathway = program,
    direction,
    NES,
    FDR,
    leading_edge_genes,
    top_core_genes,
    program_class,
    stability_score = as.integer(stability_consistent_runs),
    stable_program
  )]

  prior_gene <- gene_core_final[, .(feature = gene, type = "gene", direction, evidence = "stable")]
  prior_pathway <- prog_final[, .(
    feature = program,
    type = "pathway",
    direction = fifelse(direction == "LA-biased", "LA-high", "RA-high"),
    evidence = fifelse(FDR < 0.001 & stable_program == TRUE, "strong", fifelse(stable_program == TRUE, "moderate", "weak"))
  )]
  normal_chamber_prior <- unique(rbindlist(list(prior_gene, prior_pathway), fill = TRUE))

  plot_figure1b(gene_out, expr_mat, meta, file.path(fig_dir, "figure1B_gse112339_volcano_heatmap.pdf"), gene_fdr_cutoff, top_heatmap_genes_per_side, 0.5)
  plot_figure1b_compact(gene_out, expr_mat, meta, file.path(fig_dir, "figure1B_gse112339_volcano_heatmap_compact.pdf"), gene_fdr_cutoff, 6, 0.5)

  fwrite(gene_final, file.path(out_dir, "normal_chamber_genes.tsv"), sep = "\t", na = "NA")
  fwrite(gene_final, file.path(out_dir, "normal_chamber_genes_full.tsv"), sep = "\t", na = "NA")
  fwrite(gene_core_final, file.path(out_dir, "normal_chamber_genes_core.tsv"), sep = "\t", na = "NA")
  fwrite(gene_core_strict_final, file.path(out_dir, "normal_chamber_genes_core_strict.tsv"), sep = "\t", na = "NA")
  fwrite(prog_final, file.path(out_dir, "normal_chamber_programs.tsv"), sep = "\t", na = "NA")
  fwrite(bulk_de_full, file.path(out_dir, "bulk_chamber_DE_full.tsv"), sep = "\t", na = "NA")
  fwrite(bulk_stability, file.path(out_dir, "bulk_chamber_stability.tsv"), sep = "\t", na = "NA")
  fwrite(bulk_core_genes, file.path(out_dir, "bulk_chamber_core_genes.tsv"), sep = "\t", na = "NA")
  fwrite(bulk_fgsea, file.path(out_dir, "bulk_chamber_fgsea.tsv"), sep = "\t", na = "NA")
  fwrite(normal_chamber_prior, file.path(out_dir, "NormalChamberPrior.tsv"), sep = "\t", na = "NA")

  message("Done. Outputs:")
  for (f in c(
    "normal_chamber_genes.tsv",
    "normal_chamber_genes_full.tsv",
    "normal_chamber_genes_core.tsv",
    "normal_chamber_genes_core_strict.tsv",
    "normal_chamber_programs.tsv",
    "bulk_chamber_DE_full.tsv",
    "bulk_chamber_stability.tsv",
    "bulk_chamber_core_genes.tsv",
    "bulk_chamber_fgsea.tsv",
    "NormalChamberPrior.tsv"
  )) {
    message(" - ", file.path(out_dir, f))
  }
  message(" - ", file.path(fig_dir, "figure1A_gse112339_pca_la_vs_ra.pdf"))
  message(" - ", file.path(fig_dir, "figure1B_gse112339_volcano_heatmap.pdf"))
  message(" - ", file.path(fig_dir, "figure1B_gse112339_volcano_heatmap_compact.pdf"))
}
