#!/usr/bin/env Rscript
# S3 + S4 — pseudobulk chamber (LA vs RA) DE and pre-ranked FGSEA.
#
# For each top cell type:
#   1) read pseudobulk/{ct}/counts.tsv + samples.tsv
#   2) edgeR filterByExpr + TMM + voom + limma  with design ~ donor + chamber
#   3) write {short}_chamber_DE.tsv
#   4) pre-ranked FGSEA on moderated_t against Hallmark + Reactome
#   5) write {short}_chamber_fgsea.tsv

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(fgsea)
  library(data.table)
})

set.seed(42)

PB_ROOT  <- "/home/bio/Project_P02/pseudobulk"
OUT_ROOT <- "/home/bio/Project_P02/results"
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

GMT <- list(
  Hallmark = "/home/bio/Project_P02/h.all.v2026.1.Hs.symbols.gmt",
  Reactome = "/home/bio/Project_P02/c2.cp.reactome.v2026.1.Hs.symbols.gmt"
)

CT_MAP <- list(
  atrial_cardiomyocyte    = "atrialCM",
  fibroblast              = "fibro",
  endothelial_or_vascular = "endo",
  immune                  = "immune"
)

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
              paste0(..., collapse = " ")))
}

# ---- load gene sets once
log_msg("loading gene sets...")
gs_hallmark <- gmtPathways(GMT$Hallmark)
gs_reactome <- gmtPathways(GMT$Reactome)
gene_sets <- c(gs_hallmark, gs_reactome)
source_tag <- c(rep("Hallmark", length(gs_hallmark)),
                rep("Reactome", length(gs_reactome)))
names(source_tag) <- names(gene_sets)
log_msg("  Hallmark:", length(gs_hallmark),
        " Reactome:", length(gs_reactome),
        " total:", length(gene_sets))

run_ct <- function(ct_long, ct_short) {
  log_msg("========", ct_long, "->", ct_short, "========")
  ct_dir <- file.path(PB_ROOT, ct_long)
  counts <- as.matrix(
    read.table(file.path(ct_dir, "counts.tsv"),
               sep = "\t", header = TRUE, row.names = 1,
               check.names = FALSE)
  )
  samples <- read.table(file.path(ct_dir, "samples.tsv"),
                        sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE)

  # align columns
  counts <- counts[, samples$sample_id]
  stopifnot(all(colnames(counts) == samples$sample_id))
  log_msg("  counts:", nrow(counts), "genes x", ncol(counts), "samples")

  samples$donor   <- factor(samples$donor)
  samples$chamber <- factor(samples$chamber, levels = c("LA", "RA"))  # RA coef = RA vs LA

  design <- model.matrix(~ donor + chamber, data = samples)
  cham_col <- grep("chamber", colnames(design), value = TRUE)
  stopifnot(length(cham_col) == 1L)
  log_msg("  design rank =", qr(design)$rank, "/", ncol(design),
          " contrast coef =", cham_col)

  # ---- edgeR + voom + limma
  y <- DGEList(counts = counts)
  keep <- filterByExpr(y, design = design)
  log_msg("  filterByExpr kept", sum(keep), "/", length(keep), "genes")
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  v <- voom(y, design = design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  tt <- topTable(fit, coef = cham_col, number = Inf, sort.by = "none",
                 confint = FALSE)
  tt$gene <- rownames(tt)
  tt$direction <- ifelse(tt$logFC > 0, "RA_up", "LA_up")
  de_out <- data.frame(
    gene         = tt$gene,
    logFC        = tt$logFC,
    AveExpr      = tt$AveExpr,
    moderated_t  = tt$t,
    P.Value      = tt$P.Value,
    FDR          = tt$adj.P.Val,
    direction    = tt$direction,
    B            = tt$B,
    stringsAsFactors = FALSE
  )
  de_out <- de_out[order(de_out$FDR, de_out$P.Value), ]
  de_path <- file.path(OUT_ROOT, paste0(ct_short, "_chamber_DE.tsv"))
  write.table(de_out, de_path, sep = "\t", row.names = FALSE, quote = FALSE)
  log_msg("  wrote", de_path,
          "  [n_sig_FDR05 =", sum(de_out$FDR < 0.05, na.rm = TRUE),
          ", n_sig_FDR10 =", sum(de_out$FDR < 0.10, na.rm = TRUE), "]")

  # ---- FGSEA pre-ranked on moderated_t
  ranks <- de_out$moderated_t
  names(ranks) <- de_out$gene
  ranks <- ranks[is.finite(ranks)]
  # deduplicate just in case (gene symbols should be unique already)
  if (any(duplicated(names(ranks)))) {
    ranks <- tapply(ranks, names(ranks), max)  # keep largest-magnitude t
  }
  ranks <- sort(ranks, decreasing = TRUE)
  log_msg("  ranks: n =", length(ranks),
          " range = [", signif(min(ranks), 3), ",", signif(max(ranks), 3), "]")

  fres <- fgsea(pathways = gene_sets,
                stats    = ranks,
                minSize  = 10,
                maxSize  = 500,
                eps      = 0)
  fres$direction <- ifelse(fres$NES > 0, "RA_up", "LA_up")
  fres$source    <- source_tag[fres$pathway]
  fres$leading_edge_genes <- vapply(
    fres$leadingEdge, function(v) paste(v, collapse = ","), character(1)
  )

  fres_out <- as.data.frame(fres[, c("pathway", "source", "size",
                                     "NES", "ES", "pval", "padj",
                                     "direction", "leading_edge_genes")])
  fres_out <- fres_out[order(fres_out$padj, -abs(fres_out$NES)), ]
  fg_path <- file.path(OUT_ROOT, paste0(ct_short, "_chamber_fgsea.tsv"))
  write.table(fres_out, fg_path, sep = "\t", row.names = FALSE, quote = FALSE)
  log_msg("  wrote", fg_path,
          "  [padj<0.05:", sum(fres_out$padj < 0.05, na.rm = TRUE),
          ", padj<0.25:", sum(fres_out$padj < 0.25, na.rm = TRUE), "]")

  invisible(list(de = de_out, fgsea = fres_out))
}

results <- list()
for (ct_long in names(CT_MAP)) {
  results[[ct_long]] <- run_ct(ct_long, CT_MAP[[ct_long]])
}

log_msg("all done.")
