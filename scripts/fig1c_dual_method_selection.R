#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Figure 1c — Curated discovery primitive: dual-method (GSEA + GSVA)
# pathway selection with cross-resampling stability gate.
#
# Forward analytical operation (no downstream results referenced):
#   Input  : healthy GSE112339 paired LA / RA log-CPM (8 donors x 2 chambers)
#            + paired-limma DE summary (logFC, t-stat, FDR per gene)
#            + Hallmark + Reactome MSigDB GMT files
#   Step 1 : signed GSEA on the t-statistic ranking      -> NES per pathway
#   Step 2 : GSVA per-sample scoring + paired limma      -> logFC per pathway
#   Step 3 : 8 leave-one-donor-out fgsea resamplings      -> stability count
#   Output : combined per-pathway table with NES, GSVA logFC, stability;
#            scatter plot (NES x GSVA logFC; colour = stability;
#            highlight 10 retained programs that pass the stability gate
#            stability_consistent_runs >= 7).
#
# Files written:
#   results/fig1c_pathway_selection_long.tsv      (all candidate pathways)
#   results/fig1c_pathway_selection_retained.tsv  (10 retained programs)
#   figures/fig1c_dual_method_selection.pdf       (Fig 1c PDF panel)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(GSVA)
  library(limma)
  library(edgeR)
})

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
ROOT  <- find_root()

GMT_HALL  <- file.path(ROOT, "bulk", "Project_P01", "gene_sets",
                       "h.all.v2026.1.Hs.symbols.gmt")
GMT_REAC  <- file.path(ROOT, "bulk", "Project_P01", "gene_sets",
                       "c2.cp.reactome.v2026.1.Hs.symbols.gmt")
DE_PATH   <- file.path(ROOT, "bulk", "Project_P01", "bulk_chamber_DE_full.tsv")
COUNTS_PATH <- file.path(ROOT, "bulk", "Project_P01", "GSE112339",
                         "GSE112339_HTSeq_hg19_Ensembl.txt.gz")

OUT_LONG  <- file.path(ROOT, "results", "fig1c_pathway_selection_long.tsv")
OUT_PANEL <- file.path(ROOT, "results", "fig1c_ncp_panel_v2.tsv")
FIG_OUT   <- file.path(ROOT, "figures", "fig1c_dual_method_selection.pdf")

STABILITY_GATE <- 7L     # retained programs require >= 7 of 8 LOO resamplings
MIN_PATHWAY    <- 15L
MAX_PATHWAY    <- 500L
SET_SEED       <- 42L
TOP_N_TO_SHOW  <- 80L    # top candidates to retain in the panel

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
cat("[1/5] loading inputs ...\n")
de <- fread(DE_PATH)
de[, gene := gsub('"', '', gene_symbol)]
de <- de[!is.na(t_stat) & !is.na(gene) & nzchar(gene)]
ranked <- setNames(de$t_stat, de$gene)
ranked <- ranked[!duplicated(names(ranked))]
ranked <- sort(ranked, decreasing = TRUE)

raw <- fread(COUNTS_PATH)
gene_col <- names(raw)[1]
keep_cols <- c(gene_col, grep("^(LA|RA)_", names(raw), value = TRUE))
raw <- raw[, ..keep_cols]
mat <- as.matrix(raw[, -1, with = FALSE])
rownames(mat) <- raw[[gene_col]]
mat <- mat[!duplicated(rownames(mat)), ]
storage.mode(mat) <- "integer"
cat(sprintf("  raw counts: %d genes x %d samples\n", nrow(mat), ncol(mat)))

dge <- DGEList(counts = mat)
keep <- filterByExpr(dge, min.count = 5)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
dge  <- calcNormFactors(dge, method = "TMM")
log_cpm <- cpm(dge, log = TRUE, prior.count = 1)
cat(sprintf("  log-CPM: %d genes x %d samples\n",
            nrow(log_cpm), ncol(log_cpm)))

# Build sample table
sample_id <- colnames(log_cpm)
chamber <- substr(sample_id, 1, 2)
donor   <- sub("^(LA|RA)_", "", sample_id)
sample_meta <- data.table(sample_id = sample_id,
                           chamber = chamber, donor = donor)
cat(sprintf("  paired donors: %d (LA n=%d, RA n=%d)\n",
            length(unique(donor)),
            sum(chamber == "LA"), sum(chamber == "RA")))

pathways <- c(
  setNames(fgsea::gmtPathways(GMT_HALL),
           paste0("HALLMARK__", names(fgsea::gmtPathways(GMT_HALL)))),
  setNames(fgsea::gmtPathways(GMT_REAC),
           paste0("REACTOME__", names(fgsea::gmtPathways(GMT_REAC))))
)
cat(sprintf("  pathways loaded: %d (Hallmark + Reactome)\n", length(pathways)))

# ---------------------------------------------------------------------------
# 2. Step 1 - full fgsea
# ---------------------------------------------------------------------------
cat("[2/5] step 1 - signed GSEA on full ranking ...\n")
set.seed(SET_SEED)
fg_full <- as.data.table(fgsea(pathways = pathways, stats = ranked,
                                minSize = MIN_PATHWAY, maxSize = MAX_PATHWAY,
                                nproc = 4))
fg_full[, direction := ifelse(NES >= 0, "LA-biased", "RA-biased")]
cat(sprintf("  fgsea full: %d pathways tested, %d with padj < 0.10\n",
            nrow(fg_full), sum(fg_full$padj < 0.10, na.rm = TRUE)))

# Top by |NES|, padj < 0.25 to focus the panel
candidates <- fg_full[padj < 0.25][order(-abs(NES))][seq_len(min(TOP_N_TO_SHOW, .N))]
cat(sprintf("  retained for panel: %d candidates\n", nrow(candidates)))

# ---------------------------------------------------------------------------
# 3. Step 2 - GSVA per-sample scoring + paired limma
# ---------------------------------------------------------------------------
cat("[3/5] step 2 - GSVA per-sample scoring (this may take a minute) ...\n")
# Restrict GSVA to candidate pathways for speed
candidate_paths <- pathways[candidates$pathway]
param <- gsvaParam(exprData = log_cpm, geneSets = candidate_paths,
                   minSize = MIN_PATHWAY, maxSize = MAX_PATHWAY)
gsva_scores <- gsva(param, verbose = FALSE)
cat(sprintf("  GSVA: %d pathways x %d samples\n",
            nrow(gsva_scores), ncol(gsva_scores)))

# paired limma: chamber effect with donor as block
design <- model.matrix(~ chamber + donor, data = sample_meta)
fit <- lmFit(gsva_scores, design)
fit <- eBayes(fit)
gsva_top <- topTable(fit, coef = "chamberRA", number = Inf, sort.by = "none")
gsva_dt  <- data.table(pathway = rownames(gsva_top),
                        # negate because chamberRA coef is RA - LA;
                        # we want LA - RA to align with NES sign convention
                        gsva_logFC = -gsva_top$logFC,
                        gsva_p     = gsva_top$P.Value,
                        gsva_padj  = gsva_top$adj.P.Val)
cat(sprintf("  GSVA paired limma: %d pathways with padj < 0.10\n",
            sum(gsva_dt$gsva_padj < 0.10, na.rm = TRUE)))

# ---------------------------------------------------------------------------
# 4. Step 3 - leave-one-donor-out fgsea x 8 -> stability count
# ---------------------------------------------------------------------------
cat("[4/5] step 3 - leave-one-donor-out fgsea x 8 ...\n")
unique_donors <- sort(unique(donor))
stability <- integer(length(candidates$pathway))
names(stability) <- candidates$pathway

# Recompute paired DE per LOO subset to get rank vector
make_rank_loo <- function(drop_donor) {
  keep_idx <- which(sample_meta$donor != drop_donor)
  m <- log_cpm[, keep_idx, drop = FALSE]
  meta <- sample_meta[keep_idx]
  dsg <- model.matrix(~ chamber + donor, data = meta)
  fit_loo <- lmFit(m, dsg)
  fit_loo <- eBayes(fit_loo)
  tt <- topTable(fit_loo, coef = "chamberRA", number = Inf, sort.by = "none")
  # negate so positive = LA-biased (matches NES convention)
  rk <- -tt$t
  names(rk) <- rownames(tt)
  rk <- rk[!is.na(rk) & !duplicated(names(rk))]
  sort(rk, decreasing = TRUE)
}

set.seed(SET_SEED)
for (k in seq_along(unique_donors)) {
  d_drop <- unique_donors[k]
  rk_loo <- make_rank_loo(d_drop)
  fg_loo <- as.data.table(fgsea(pathways = pathways[candidates$pathway],
                                 stats = rk_loo,
                                 minSize = MIN_PATHWAY, maxSize = MAX_PATHWAY,
                                 nproc = 4))
  # Stability requires GSEA significance at padj < 0.05 AND same sign as the
  # full-sample NES; this matches the project's original stability convention.
  hits <- fg_loo[padj < 0.05 & sign(NES) == sign(candidates$NES[match(pathway, candidates$pathway)]),
                 pathway]
  stability[hits] <- stability[hits] + 1L
  cat(sprintf("  LOO drop %s : %d / %d candidates reproduced\n",
              d_drop, length(hits), length(candidates$pathway)))
}
candidates[, stability_consistent_runs := stability[pathway]]
candidates[, retained := stability_consistent_runs >= STABILITY_GATE]
cat(sprintf("  retained at stability >= %d: %d pathways\n",
            STABILITY_GATE, sum(candidates$retained)))

# ---------------------------------------------------------------------------
# 5. Merge GSVA + Stability into one table; write outputs
# ---------------------------------------------------------------------------
cat("[5/5] merging and writing outputs ...\n")
out <- merge(candidates[, .(pathway, NES, padj, direction,
                             stability_consistent_runs, retained)],
             gsva_dt, by = "pathway", all.x = TRUE)
out[, abs_NES := abs(NES)]
setorder(out, -abs_NES)
out[, abs_NES := NULL]
fwrite(out, OUT_LONG, sep = "\t")
cat(sprintf("  wrote %s (%d candidate rows)\n",
            sub(paste0("^", ROOT, "/"), "", OUT_LONG), nrow(out)))

# ---------------------------------------------------------------------------
# 5b. Build the NEW NCP panel from this fresh analysis.  No inheritance of
#     old project locks — the panel size K is data-driven.
#
#     Selection criteria:
#       (i)   GSEA padj < 0.01
#       (ii)  GSVA padj < 0.01
#       (iii) sign(NES) == sign(GSVA logFC)            (same direction)
#       (iv)  stability_consistent_runs == 8           (8 / 8 LOO stable)
#       (v)   Jaccard redundancy clustering at distance >= 0.5 on gene-set
#             leading-edge / member overlap; LA and RA branches are
#             clustered separately so the chamber direction is preserved.
#       (vi)  per cluster, the representative is the pathway with the
#             largest |NES|.
#
#     The resulting K is the data-derived NCP panel size; all downstream
#     analyses are intended to be re-run on this panel.
# ---------------------------------------------------------------------------
NEW_PANEL_GATE <- list(gsea_padj = 0.01, gsva_padj = 0.01, stab = 8L,
                        jaccard_height = 0.5)

retained <- out[padj < NEW_PANEL_GATE$gsea_padj &
                gsva_padj < NEW_PANEL_GATE$gsva_padj &
                stability_consistent_runs >= NEW_PANEL_GATE$stab &
                sign(NES) == sign(gsva_logFC)]
cat(sprintf("  dual-method + stability gate: %d / %d candidates pass\n",
            nrow(retained), nrow(out)))

# Jaccard distance matrix per direction
jaccard_dist <- function(gs_list) {
  n <- length(gs_list)
  J <- matrix(0, n, n, dimnames = list(names(gs_list), names(gs_list)))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    a <- gs_list[[i]]; b <- gs_list[[j]]
    J[i, j] <- 1 - length(intersect(a, b)) / length(union(a, b))
  }
  J
}

cluster_per_direction <- function(d_sub) {
  if (nrow(d_sub) < 2) {
    d_sub[, cluster_id := 1L]
    return(d_sub)
  }
  gs <- pathways[d_sub$pathway]
  J <- jaccard_dist(gs)
  hc <- hclust(as.dist(J), method = "average")
  d_sub[, cluster_id := cutree(hc, h = NEW_PANEL_GATE$jaccard_height)]
  d_sub
}

retained <- rbindlist(list(
  cluster_per_direction(retained[direction == "LA-biased"])[
                             , cluster_id := paste0("LA_", cluster_id)],
  cluster_per_direction(retained[direction == "RA-biased"])[
                             , cluster_id := paste0("RA_", cluster_id)]
))

retained[, abs_NES := abs(NES)]
ncp_panel <- retained[, .SD[which.max(abs_NES)], by = cluster_id]
setorder(ncp_panel, -abs_NES)
ncp_panel[, program_order := seq_len(.N)]

cat(sprintf("  Jaccard clustering at h=%.2f, LA + RA separately\n",
            NEW_PANEL_GATE$jaccard_height))
cat(sprintf("  -> %d clusters = NEW NCP panel (K = %d)\n",
            nrow(ncp_panel), nrow(ncp_panel)))

# Write the NEW NCP panel (data-driven; the panel used by all downstream
# analyses going forward).
fwrite(ncp_panel[, .(program_order, cluster_id, pathway, direction,
                     NES, padj, gsva_logFC, gsva_padj,
                     stability_consistent_runs)],
       OUT_PANEL, sep = "\t")
cat(sprintf("  wrote %s (%d NCP_v2 panel members)\n",
            sub(paste0("^", ROOT, "/"), "", OUT_PANEL), nrow(ncp_panel)))

# Tag the retained panel back into the candidate table
out[, ncp_panel := pathway %in% ncp_panel$pathway]
out <- merge(out, ncp_panel[, .(pathway, program_order, cluster_id)],
             by = "pathway", all.x = TRUE)

# ---------------------------------------------------------------------------
# 6. Plot
# ---------------------------------------------------------------------------
cat("plotting Fig 1c ...\n")

# Locked display convention (updated):
#   RA -> blue, LA -> red
COL_RA <- "#4393C3"
COL_LA <- "#D6604D"
COL_SINGLE <- 89 / 25.4

# Wider / taller layout for readability in downstream montage assembly.
pdf(FIG_OUT, width = COL_SINGLE * 2.25, height = COL_SINGLE * 1.22,
    useDingbats = FALSE, pointsize = 8)
layout(matrix(c(1, 2), nrow = 1), widths = c(1.12, 1.08))

# ---- Left panel : scatter ----
par(mar = c(4.3, 4.4, 2.8, 0.8), mgp = c(2.6, 0.65, 0), tck = -0.02,
    cex.axis = 0.95, cex.lab = 1.05)

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
      cex.main = 0.92, font.main = 1, line = 0.5)

# Highlight the K NCP panel members (data-driven), without numeric IDs.
panel_pts <- out[ncp_panel == TRUE]
setorder(panel_pts, program_order)
points(panel_pts$NES, panel_pts$gsva_logFC, pch = 21,
       cex = 1.55, bg = ifelse(panel_pts$NES > 0, COL_LA, COL_RA),
       col = "black", lwd = 0.7)

# Stability + NCP panel legend in left panel
legend("topleft", inset = c(0.0, -0.02),
       legend = c("0/8", "4/8", "8/8",
                  "NCP (LA)", "NCP (RA)"),
       pch = c(21, 21, 21, 21, 21),
       pt.bg = c(stab_colors[1], stab_colors[5], stab_colors[9],
                 COL_LA, COL_RA),
       pt.cex = c(0.95, 0.95, 0.95, 1.35, 1.35),
       col = c("white", "white", "white", "black", "black"),
       pt.lwd = 0.3, bty = "n", cex = 0.76, y.intersp = 1.05)

# ---- Right panel : pathway label column (no numeric markers) -----------
par(mar = c(4.3, 0.1, 2.8, 0.5))
plot.new()
K <- nrow(panel_pts)
plot.window(xlim = c(0, 1), ylim = c(0, K + 1))
title(main = sprintf("NCP pathways (K = %d)", K),
      cex.main = 0.92, font.main = 1, line = 0.5)

# Pathway abbreviation table (parallel to scripts/_pathway_labels.py)
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

# Two-column layout when K > 18
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
         pos = 4, cex = 0.74, col = "black", font = 1, offset = 0)
  }
} else {
  for (i in seq_len(K)) {
    y <- (K + 0.5) - (i - 1) * 0.95
    is_la <- panel_pts$NES[i] > 0
    points(0.04, y, pch = 21, cex = 1.55,
           bg = ifelse(is_la, COL_LA, COL_RA),
           col = "black", lwd = 0.7)
    text(0.10, y, labels = short_label(panel_pts$pathway[i]),
         pos = 4, cex = 0.86, col = "black", font = 1, offset = 0)
  }
}

dev.off()
cat(sprintf("wrote %s\n",
            sub(paste0("^", ROOT, "/"), "", FIG_OUT)))
cat("\nDONE.\n")
