#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Figure 1c — replot only.
# Reuses cached selection tables from fig1c_dual_method_selection.R:
#   results/fig1c_pathway_selection_long.tsv
#   results/fig1c_ncp_panel_v2.tsv
# Re-renders figures/fig1c_dual_method_selection.pdf with:
#   - bigger fonts / clearer layout
#   - no numeric labels on NCP markers
#   - RA -> blue (#4393C3), LA -> red (#D6604D)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

find_root <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("--file=", "", grep("--file=", ca, value = TRUE)[1])
  if (length(f) && !is.na(f) && nzchar(f)) {
    d <- normalizePath(dirname(f), mustWork = FALSE)
    if (basename(d) == "scripts") return(dirname(d))
    return(d)
  }
  cwd <- getwd(); if (basename(cwd) == "scripts") return(dirname(cwd)); cwd
}
ROOT <- find_root()

LONG_TSV  <- file.path(ROOT, "results", "fig1c_pathway_selection_long.tsv")
PANEL_TSV <- file.path(ROOT, "results", "fig1c_ncp_panel_v2.tsv")
FIG_OUT   <- file.path(ROOT, "figures", "fig1c_dual_method_selection.pdf")

stopifnot(file.exists(LONG_TSV), file.exists(PANEL_TSV))

out  <- fread(LONG_TSV)
ncp  <- fread(PANEL_TSV)

# Tag panel membership in `out` if missing
if (!"ncp_panel" %in% names(out)) out[, ncp_panel := pathway %in% ncp$pathway]
if (!"program_order" %in% names(out)) {
  out <- merge(out,
               ncp[, .(pathway, program_order)],
               by = "pathway", all.x = TRUE)
}

panel_pts <- out[ncp_panel == TRUE]
setorder(panel_pts, program_order)

# Locked display convention (updated):
#   RA -> blue, LA -> red
COL_RA <- "#4393C3"
COL_LA <- "#D6604D"

PATHWAY_SHORT <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"          = "OXPHOS",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT"     = "ETC",
  "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT" = "Aerobic resp / ETC",
  "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE"        = "TCA cycle",
  "REACTOME_COMPLEX_I_BIOGENESIS"               = "Complex I biogenesis",
  "REACTOME_COMPLEX_III_ASSEMBLY"               = "Complex III assembly",
  "REACTOME_COMPLEX_IV_ASSEMBLY"                = "Complex IV assembly",
  "REACTOME_MITOCHONDRIAL_TRANSLATION"          = "Mito translation",
  "REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT"       = "Mito protein import",
  "REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION"  = "Mito protein degradation",
  "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION" = "Mito FA \u03b2-oxidation",
  "REACTOME_MITOCHONDRIAL_TRNA_AMINOACYLATION"  = "Mito tRNA aminoacyl.",
  "REACTOME_MITOCHONDRIAL_CALCIUM_ION_TRANSPORT" = "Mito Ca\u00b2\u207a transport",
  "REACTOME_TRANSLATION"                        = "Translation (cyto)",
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION"  = "Euk translation elong.",
  "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION"  = "Euk translation init.",
  "REACTOME_RIBOSOME_ASSOCIATED_QUALITY_CONTROL" = "Ribosome QC (RQC)",
  "REACTOME_TRNA_AMINOACYLATION"                = "tRNA aminoacylation",
  "REACTOME_PROTEIN_LOCALIZATION"               = "Protein localization",
  "REACTOME_HOST_INTERACTIONS_OF_HIV_FACTORS"   = "HIV host interactions",
  "REACTOME_SARS_COV_2_MODULATES_HOST_TRANSLATION_MACHINERY" = "SARS-CoV-2 host translation",
  "REACTOME_TRIGLYCERIDE_METABOLISM"            = "TG metabolism",
  "HALLMARK_FATTY_ACID_METABOLISM"              = "FA metabolism",
  "REACTOME_BRANCHED_CHAIN_AMINO_ACID_CATABOLISM" = "BCAA catabolism",
  "HALLMARK_ADIPOGENESIS"                       = "Adipogenesis",
  "HALLMARK_PEROXISOME"                         = "Peroxisome",
  "REACTOME_PEROXISOMAL_PROTEIN_IMPORT"         = "Peroxisomal protein import",
  "REACTOME_HEPARAN_SULFATE_HEPARIN_HS_GAG_METABOLISM" = "HS-GAG metabolism",
  "REACTOME_HS_GAG_BIOSYNTHESIS"                = "HS-GAG biosynthesis",
  "REACTOME_DISEASES_ASSOCIATED_WITH_GLYCOSAMINOGLYCAN_METABOLISM" = "GAG metabolism diseases",
  "REACTOME_RUNX2_REGULATES_BONE_DEVELOPMENT"   = "RUNX2 / bone dev.",
  "REACTOME_SIGNALING_BY_BMP"                   = "BMP signalling",
  "HALLMARK_TGF_BETA_SIGNALING"                 = "TGF-\u03b2 signalling",
  "REACTOME_CLASS_B_2_SECRETIN_FAMILY_RECEPTORS" = "Class B/2 Secretin receptors"
)

short_label <- function(s) {
  s2 <- sub("^HALLMARK__HALLMARK_", "HALLMARK_", s)
  s2 <- sub("^REACTOME__REACTOME_", "REACTOME_", s2)
  if (s2 %in% names(PATHWAY_SHORT)) return(PATHWAY_SHORT[[s2]])
  s2 <- sub("^HALLMARK_", "", s2)
  s2 <- sub("^REACTOME_", "", s2)
  s2 <- gsub("_", " ", s2)
  if (nchar(s2) > 32) s2 <- paste0(substr(s2, 1, 30), "\u2026")
  s2
}

cat("re-plotting Fig 1c (cached) ...\n")
COL_SINGLE <- 89 / 25.4

cairo_pdf(FIG_OUT, width = COL_SINGLE * 2.25, height = COL_SINGLE * 1.22,
          pointsize = 8)
layout(matrix(c(1, 2), nrow = 1), widths = c(1.12, 1.08))

# ---- Left panel : scatter ----
par(mar = c(4.3, 4.4, 2.8, 0.8), mgp = c(2.6, 0.65, 0), tck = -0.02,
    cex.axis = 1.05, cex.lab = 1.15)

stab_colors <- colorRampPalette(c("#CCCCCC", "#999999",
                                   "#666666", "#333333"))(9)
xrng <- range(out$NES, na.rm = TRUE)
yrng <- range(out$gsva_logFC, na.rm = TRUE)

plot(out$NES, out$gsva_logFC,
     pch = 21, cex = 0.72,
     bg = stab_colors[pmin(out$stability_consistent_runs + 1, 9)],
     col = "white", lwd = 0.3,
     xlim = xrng, ylim = yrng,
     xlab = "signed GSEA NES (LA-biased > 0)",
     ylab = "GSVA paired logFC (LA - RA, log2)",
     main = "")
abline(h = 0, col = "#888888", lwd = 0.4)
abline(v = 0, col = "#888888", lwd = 0.4)

title(main = sprintf("c  Dual-method selection (K = %d)", nrow(panel_pts)),
      cex.main = 1.02, font.main = 1, line = 0.5)

# Highlight the K NCP panel members (no numeric labels)
points(panel_pts$NES, panel_pts$gsva_logFC, pch = 21,
       cex = 1.55, bg = ifelse(panel_pts$NES > 0, COL_LA, COL_RA),
       col = "black", lwd = 0.7)

legend("topleft", inset = c(0.0, -0.02),
       legend = c("0/8", "4/8", "8/8",
                  "NCP (LA)", "NCP (RA)"),
       pch = c(21, 21, 21, 21, 21),
       pt.bg = c(stab_colors[1], stab_colors[5], stab_colors[9],
                 COL_LA, COL_RA),
       pt.cex = c(0.95, 0.95, 0.95, 1.35, 1.35),
       col = c("white", "white", "white", "black", "black"),
       pt.lwd = 0.3, bty = "n", cex = 0.86, y.intersp = 1.05)

# ---- Right panel : pathway label column (no numeric markers) ----
par(mar = c(4.3, 0.1, 2.8, 0.5))
plot.new()
K <- nrow(panel_pts)
plot.window(xlim = c(0, 1), ylim = c(0, K + 1))
title(main = sprintf("NCP pathways (K = %d)", K),
      cex.main = 1.02, font.main = 1, line = 0.5)

if (K > 18) {
  half <- ceiling(K / 2)
  for (i in seq_len(K)) {
    in_left <- i <= half
    col_x   <- if (in_left) 0.0 else 0.50
    row_i   <- if (in_left) i else (i - half)
    y <- (half + 0.5) - (row_i - 1) * (half / max(K - half, 1))
    is_la <- panel_pts$NES[i] > 0
    points(col_x + 0.04, y, pch = 21, cex = 1.40,
           bg = ifelse(is_la, COL_LA, COL_RA),
           col = "black", lwd = 0.6)
    text(col_x + 0.09, y,
         labels = short_label(panel_pts$pathway[i]),
         pos = 4, cex = 0.82, col = "black", font = 1, offset = 0)
  }
} else {
  for (i in seq_len(K)) {
    y <- (K + 0.5) - (i - 1) * 0.95
    is_la <- panel_pts$NES[i] > 0
    points(0.04, y, pch = 21, cex = 1.55,
           bg = ifelse(is_la, COL_LA, COL_RA),
           col = "black", lwd = 0.7)
    text(0.10, y, labels = short_label(panel_pts$pathway[i]),
         pos = 4, cex = 0.94, col = "black", font = 1, offset = 0)
  }
}

dev.off()
cat(sprintf("wrote %s\n",
            sub(paste0("^", ROOT, "/"), "", FIG_OUT)))
