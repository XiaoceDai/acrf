#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# G3b — Extended Fig 5D ROI singscores against:
#         (a) multiple NCP_v2 programs covering all 4 evidence-class families
#             from Fig 5C (CONVERGENT, BULK-LED, DURATION-LED) + an LA OXPHOS
#             control;
#         (b) multiple canonical cell-type marker panels (fibroblast,
#             lymphatic endothelial, atrial cardiomyocyte, macrophage,
#             vascular endothelial), using small explicit hand-curated panels
#             — NOT data-derived signatures.
#
# Inputs:
#   Project_P5/GSE261363/GSE261363_All_Data_WTA_1_7_22.xlsx
#   results/fig3_ncpv2_program_le.tsv     (NCP_v2 LE genes)
#   Project_P5/GSE261363/roi_metadata.tsv (from g3, has all metadata cols)
#
# Outputs:
#   Project_P5/GSE261363/roi_scores_v2.tsv     long format (gsm × signature)
#   Project_P5/GSE261363/roi_scores_v2_wide.tsv wide (gsm × signature)
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl); library(data.table); library(singscore)
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
ROOT <- find_root()
G    <- file.path(ROOT, "Project_P5", "GSE261363")

# ---- 1. Load Q3 expression matrix ----
xls  <- file.path(G, "GSE261363_All_Data_WTA_1_7_22.xlsx")
expr <- as.data.table(read_excel(xls, sheet = "Q3_5percfilter"))
expr_cols <- setdiff(names(expr), "TargetName")
expr_mat  <- as.matrix(expr[, ..expr_cols])
rownames(expr_mat) <- expr$TargetName
cat(sprintf("WTA matrix:  %d genes  ×  %d ROIs\n",
            nrow(expr_mat), ncol(expr_mat)))

# ---- 2. Load ROI metadata (keeping series_matrix → expr column order) ----
md <- fread(file.path(G, "roi_metadata.tsv"))
md_in <- md[match(expr_cols, md$title)]
stopifnot(all(md_in$title == expr_cols))

# ---- 3. NCP_v2 program panel + LE genes (R3b output) ----
le_tbl <- fread(file.path(ROOT, "results", "fig3_ncpv2_program_le.tsv"))

# Programs to score (cover the four classes in Fig 5C + 1 LA-control)
PROGRAM_TARGETS <- list(
  list(id = "REACTOME__REACTOME_DISEASES_ASSOCIATED_WITH_GLYCOSAMINOGLYCAN_METABOLISM",
       short = "F10  GAG metabolism",          family = "F10", role = "CONVERGENT (RA)"),
  list(id = "REACTOME__REACTOME_SIGNALING_BY_BMP",
       short = "F11  BMP signalling",          family = "F11", role = "CONVERGENT (RA)"),
  list(id = "REACTOME__REACTOME_HEPARAN_SULFATE_HEPARIN_HS_GAG_METABOLISM",
       short = "F9   HS-GAG metabolism",       family = "F9",  role = "BULK-LED (RA)"),
  list(id = "REACTOME__REACTOME_MITOCHONDRIAL_CALCIUM_ION_TRANSPORT",
       short = "F2   Mito Ca²⁺ transport",     family = "F2",  role = "DURATION-LED (LA)"),
  list(id = "HALLMARK__HALLMARK_OXIDATIVE_PHOSPHORYLATION",
       short = "F7   OXPHOS (LA control)",     family = "F7",  role = "trend (LA)")
)

# Resolve LE genes per target program
le_lookup <- setNames(strsplit(le_tbl$leading_edge_genes, ",", fixed = TRUE),
                       le_tbl$program_id)

# ---- 4. Canonical cell-type marker panels (small, hand-curated, declared) ----
CELL_PANELS <- list(
  fibroblast            = c("PDGFRA", "COL1A1", "COL1A2", "COL3A1", "COL5A1",
                              "DCN", "LUM", "POSTN", "FAP", "FBN1", "MMP2", "BGN"),
  lymphatic_endothelial = c("PROX1", "PDPN", "FLT4", "LYVE1", "CCL21",
                              "CLDN5", "TFF3"),
  atrialCM              = c("MYH6", "MYL7", "NPPA", "NPPB", "TNNT2",
                              "ACTC1", "RYR2", "PLN", "ATP2A2", "GJA5",
                              "KCNJ3", "KCNA5"),
  macrophage            = c("CD68", "CD163", "MARCO", "C1QA", "C1QB",
                              "MRC1", "MS4A6A", "TYROBP", "LYZ"),
  vascular_endothelial  = c("PECAM1", "VWF", "CDH5", "KDR", "FLT1",
                              "TIE1", "CLDN5", "NOTCH4")
)
cat("\nCanonical cell-type panels (n_genes × n_in_WTA):\n")
for (n in names(CELL_PANELS)) {
  in_wta <- sum(CELL_PANELS[[n]] %in% rownames(expr_mat))
  cat(sprintf("  %-22s   total=%2d   in_WTA=%2d\n",
              n, length(CELL_PANELS[[n]]), in_wta))
}

# ---- 5. Score everything ----
rk <- rankGenes(expr_mat)
score_set <- function(rk, gs) {
  gs_in <- intersect(gs, rownames(rk))
  if (length(gs_in) < 3L) {
    return(list(score = rep(NA_real_, ncol(rk)),
                n_in = length(gs_in)))
  }
  s <- simpleScore(rk, upSet = gs_in,
                   knownDirection = TRUE, centerScore = TRUE)
  list(score = s$TotalScore, n_in = length(gs_in))
}

rows <- list()
# programs
for (p in PROGRAM_TARGETS) {
  gs <- if (!is.null(le_lookup[[p$id]]))
          unique(trimws(le_lookup[[p$id]])) else character(0)
  res <- score_set(rk, gs)
  rows[[length(rows) + 1]] <- data.table(
    signature_type = "program",
    signature_id   = p$id,
    signature_lab  = p$short,
    family         = p$family,
    evidence_role  = p$role,
    n_genes_total  = length(gs),
    n_genes_in_WTA = res$n_in,
    gsm            = md_in$gsm,
    title          = md_in$title,
    score          = res$score)
}
# cell-type panels
for (n in names(CELL_PANELS)) {
  gs <- CELL_PANELS[[n]]
  res <- score_set(rk, gs)
  rows[[length(rows) + 1]] <- data.table(
    signature_type = "cell_type",
    signature_id   = n,
    signature_lab  = n,
    family         = NA_character_,
    evidence_role  = NA_character_,
    n_genes_total  = length(gs),
    n_genes_in_WTA = res$n_in,
    gsm            = md_in$gsm,
    title          = md_in$title,
    score          = res$score)
}
long <- rbindlist(rows)
long <- merge(long,
              md_in[, .(gsm, patient, rhythm, tissue, area, marker)],
              by = "gsm", sort = FALSE)
long[, marker_norm := fifelse(marker %in% c("DNA", "DNA+"), "cardiac",
                       fifelse(marker == "FABP4+", "adipose", marker))]
long[, eligible_for_5D := (tissue == "Atrial appendage") &
                          (marker_norm == "cardiac") &
                          (!is.na(rhythm)) &
                          (area %in% c("distal", "border"))]

fwrite(long, file.path(G, "roi_scores_v2.tsv"), sep = "\t")
cat(sprintf("\nwrote roi_scores_v2.tsv  (%d rows)\n", nrow(long)))

# wide format for convenience
wide <- dcast(long[, .(gsm, signature_lab, score)],
              gsm ~ signature_lab, value.var = "score")
fwrite(wide, file.path(G, "roi_scores_v2_wide.tsv"), sep = "\t")
cat(sprintf("wrote roi_scores_v2_wide.tsv  (%d rows × %d cols)\n",
            nrow(wide), ncol(wide)))

# Brief summary on eligible AA cardiac ROIs
cat("\n=== AA cardiac ROI (n = 24) AF − SR Δ summary per signature ===\n")
elig <- long[eligible_for_5D == TRUE]
af_sr <- elig[, .(
    score_af  = mean(score[rhythm == "AF"], na.rm = TRUE),
    score_sr  = mean(score[rhythm == "SR"], na.rm = TRUE),
    delta_af_minus_sr = mean(score[rhythm == "AF"], na.rm = TRUE)
                        - mean(score[rhythm == "SR"], na.rm = TRUE),
    n_af = sum(rhythm == "AF" & !is.na(score)),
    n_sr = sum(rhythm == "SR" & !is.na(score))
  ), by = .(signature_type, signature_lab)]
print(af_sr[order(-abs(delta_af_minus_sr))])

cat("\nG3b DONE.\n")
