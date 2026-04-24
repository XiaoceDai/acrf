#!/usr/bin/env Rscript
# S7 - Gene-level attribution of bulk LA/RA core genes to cell types.
#
# Approach:
#   * Split bulk core genes into two sets: bulk_LA_core and bulk_RA_core.
#   * For each of the 4 cell types, rank sc DE genes by moderated_t (descending).
#   * Run fgsea with the two bulk sets against each cell-type rank.
#   * NOTE on sign: our sc DE moderated_t follows a RA-vs-LA contrast,
#     so t-descending rank puts RA_up genes on top and LA_up at bottom.
#     Hence we expect bulk_LA_core to enrich at the bottom (NES < 0) and
#     bulk_RA_core at the top (NES > 0) if the sc signal agrees with bulk.
#   * A bulk LA-high gene is "carried by" cell type C if:
#       - fgsea(bulk_LA_core, C).NES < 0 AND padj < 0.25
#       - AND the gene is in fgsea leadingEdge for (bulk_LA_core, C).
#     Mirror rule for RA-high genes (NES > 0, padj < 0.25, in leadingEdge).
#   * Labels:
#       0 carriers  -> Unresolved
#       1 carrier   -> that cell type display name
#       >=2 carriers-> Mixed
#
# Output: results/normal_chamber_gene_attribution.tsv
# Columns: gene, bulk_direction, dominant_cell_type, attribution_score, note

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

set.seed(42)

PADJ_THRESH <- 0.25
CTS         <- c("atrialCM", "fibro", "endo", "immune")
CT_DISPLAY  <- c(atrialCM = "atrialCM",
                 fibro    = "fibroblast",
                 endo     = "endothelial_or_vascular",
                 immune   = "immune")

cat("[S7] Reading bulk core genes ...\n")
bulk <- fread("normal_chamber_genes_core.tsv")
bulk[, gene := as.character(gene)]
bulk[, direction := as.character(direction)]

bulk_LA <- unique(bulk[direction == "LA-high", gene])
bulk_RA <- unique(bulk[direction == "RA-high", gene])
cat(sprintf("  bulk_LA_core: %d genes\n  bulk_RA_core: %d genes\n",
            length(bulk_LA), length(bulk_RA)))

gene_sets <- list(bulk_LA_core = bulk_LA, bulk_RA_core = bulk_RA)

run_fgsea_for_ct <- function(ct) {
  de <- fread(sprintf("results/%s_chamber_DE.tsv", ct))
  de <- de[!is.na(moderated_t)]
  de <- de[!duplicated(gene)]
  ranks <- setNames(de$moderated_t, de$gene)
  ranks <- sort(ranks, decreasing = TRUE)

  fg <- fgsea(pathways   = gene_sets,
              stats      = ranks,
              minSize    = 1,
              maxSize    = length(ranks),
              nPermSimple = 10000)
  list(fg = fg, ranks = ranks, de = de)
}

cat("[S7] Running fgsea per cell type ...\n")
ct_res <- lapply(CTS, run_fgsea_for_ct)
names(ct_res) <- CTS

for (ct in CTS) {
  fg <- ct_res[[ct]]$fg
  cat(sprintf("  %-8s | ", ct))
  for (i in seq_len(nrow(fg))) {
    cat(sprintf("%s NES=%+0.2f padj=%.2e  ",
                fg$pathway[i], fg$NES[i], fg$padj[i]))
  }
  cat("\n")
}

# Per-cell-type "is this gene carried by CT?"
# We compute a wide indicator for each bulk core gene.
cat("[S7] Determining per-gene carriers ...\n")
attr <- data.table(gene = bulk$gene, bulk_direction = bulk$direction)

for (ct in CTS) {
  fg <- ct_res[[ct]]$fg
  de <- ct_res[[ct]]$de

  row_LA <- fg[pathway == "bulk_LA_core"]
  row_RA <- fg[pathway == "bulk_RA_core"]

  le_LA <- if (nrow(row_LA) > 0) unlist(row_LA$leadingEdge) else character(0)
  le_RA <- if (nrow(row_RA) > 0) unlist(row_RA$leadingEdge) else character(0)

  # sc DE uses RA-vs-LA contrast: positive t = RA_up, negative t = LA_up.
  # With t-descending rank, bulk_LA_core should enrich at the bottom (NES<0)
  # and bulk_RA_core at the top (NES>0).
  LA_valid <- nrow(row_LA) > 0 && !is.na(row_LA$NES) &&
              row_LA$NES < 0 && !is.na(row_LA$padj) && row_LA$padj < PADJ_THRESH
  RA_valid <- nrow(row_RA) > 0 && !is.na(row_RA$NES) &&
              row_RA$NES > 0 && !is.na(row_RA$padj) && row_RA$padj < PADJ_THRESH

  # carrier indicator per gene
  attr[, (paste0("in_LE_", ct)) := ifelse(
    bulk_direction == "LA-high",
    gene %in% le_LA & LA_valid,
    gene %in% le_RA & RA_valid
  )]

  # attach moderated_t from the cell type (for attribution_score / note)
  tmap <- setNames(de$moderated_t, de$gene)
  attr[, (paste0("t_", ct)) := tmap[gene]]
}

# Count carriers
carrier_cols <- paste0("in_LE_", CTS)
attr[, n_carriers := rowSums(as.matrix(.SD), na.rm = TRUE), .SDcols = carrier_cols]

# Build dominant_cell_type, attribution_score, note
assign_label <- function(row_idx) {
  carriers <- CTS[sapply(carrier_cols, function(c) isTRUE(attr[[c]][row_idx]))]
  ts_vec <- sapply(CTS, function(ct) attr[[paste0("t_", ct)]][row_idx])

  if (length(carriers) == 0) {
    dominant <- "Unresolved"
    score <- 0
    # Note: list the strongest concordant (direction-matching) t across CTs, if any
    dir <- attr$bulk_direction[row_idx]
    # sc t uses RA-vs-LA contrast: LA_up gene => t<0, RA_up => t>0
    sign_req <- if (dir == "LA-high") -1 else 1
    concordant <- ts_vec[!is.na(ts_vec) & sign(ts_vec) == sign_req]
    if (length(concordant) > 0) {
      best <- names(concordant)[which.max(abs(concordant))]
      note <- sprintf("best_concordant=%s(t=%+0.2f)",
                      CT_DISPLAY[best], concordant[best])
    } else {
      note <- "no_concordant_signal"
    }
  } else {
    cts_t <- ts_vec[carriers]
    # strongest |t| among carriers
    score <- max(abs(cts_t), na.rm = TRUE)
    if (length(carriers) == 1) {
      dominant <- CT_DISPLAY[carriers]
    } else {
      dominant <- "Mixed"
    }
    parts <- sapply(carriers, function(ct) {
      sprintf("%s(t=%+0.2f)", CT_DISPLAY[ct], cts_t[ct])
    })
    note <- paste(parts, collapse = "; ")
  }
  list(dominant = dominant, score = score, note = note)
}

labels <- lapply(seq_len(nrow(attr)), assign_label)
attr[, dominant_cell_type := vapply(labels, `[[`, character(1), "dominant")]
attr[, attribution_score  := vapply(labels, `[[`, numeric(1),   "score")]
attr[, note                := vapply(labels, `[[`, character(1), "note")]

out <- attr[, .(gene, bulk_direction, dominant_cell_type,
                attribution_score, note)]
outpath <- "results/normal_chamber_gene_attribution.tsv"
fwrite(out, outpath, sep = "\t", quote = FALSE)
cat(sprintf("[S7] Wrote %s (%d rows)\n", outpath, nrow(out)))

cat("\n[S7] Label counts:\n")
print(out[, .N, by = .(bulk_direction, dominant_cell_type)][order(bulk_direction, -N)])
