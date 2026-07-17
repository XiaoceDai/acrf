#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# G3b v3 — Extended Fig 5d v3 ROI singscores against:
#         (a) 6 NCP axis families covering the v3 Atlas top-priority + a
#             metabolic anchor (F2 Mito Ca²⁺) + an adipocyte-carrier metabolic
#             family (F4 BCAA).  Each family signature = UNION of LE genes
#             across its NCP pathway programs (rank-singscore semantics);
#         (b) the same canonical hand-curated cell-type marker panels as v2
#             (fibroblast, lymphatic_endothelial, atrialCM, macrophage,
#             vascular_endothelial) — explicitly NOT data-derived.
#
# Inputs:
#   Project_P5/GSE261363/GSE261363_All_Data_WTA_1_7_22.xlsx  (Q3-normalised)
#   results/fig3_v3_program_le.tsv                            (NCP LE)
#   Project_P5/GSE261363/roi_metadata.tsv                     (from g3, all metadata)
#
# Outputs:
#   Project_P5/GSE261363/roi_scores_v3.tsv       long format (gsm × signature)
#   Project_P5/GSE261363/roi_scores_v3_wide.tsv  wide  (gsm × signature)
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
cat(sprintf("WTA matrix:  %d genes  x  %d ROIs\n",
            nrow(expr_mat), ncol(expr_mat)))

# ---- 2. Load ROI metadata ----
md <- fread(file.path(G, "roi_metadata.tsv"))
md_in <- md[match(expr_cols, md$title)]
stopifnot(all(md_in$title == expr_cols))

# ---- 3. NCP axis family panel (union LE per family) ----
le_tbl <- fread(file.path(ROOT, "results", "fig3_v3_program_le.tsv"))
# le_tbl columns: program_id, program_order, direction, family_id,
#                 n_full_pathway, n_LE_total, leading_edge_genes

# Family targets: 6 chosen from atlas Top + biology anchors
FAMILY_TARGETS <- list(
  list(family_id = 15, short = "F15  NTRK2 / TrkB",
       role = "Triple-positive (***)", chamber_driver = "LA-driven"),
  list(family_id = 11, short = "F11  BMP signalling",
       role = "Single-axis L1 (***)", chamber_driver = "LA-driven"),
  list(family_id = 9,  short = "F9   HS-GAG metabolism",
       role = "Double-positive L1+L3 (***)", chamber_driver = "LA-driven"),
  list(family_id = 10, short = "F10  TGF-beta signalling",
       role = "Double-positive L1+L3 (***)", chamber_driver = "LA-driven"),
  list(family_id = 2,  short = "F2   Mito Ca2+ transport",
       role = "Single-axis L3 (Strong-trend)", chamber_driver = "Bilateral"),
  list(family_id = 4,  short = "F4   BCAA catabolism",
       role = "Silent (atlas Top, Adipo carrier)", chamber_driver = "Bilateral")
)

# Build family-level LE union
family_le <- function(fid) {
  rows <- le_tbl[family_id == fid]
  unique(unlist(strsplit(rows$leading_edge_genes, ",", fixed = TRUE)))
}

# ---- 4. Canonical cell-type marker panels (unchanged from v2) ----
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
                              "TIE1", "CLDN5", "NOTCH4"),
  adipocyte             = c("ADIPOQ", "PLIN1", "PLIN4", "FABP4", "LEP",
                              "PPARG", "CEBPA", "GPD1", "LIPE", "CIDEA")
)

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

cat("\n--- Family signatures (NCP LE union) ---\n")
rows <- list()
for (p in FAMILY_TARGETS) {
  gs <- family_le(p$family_id)
  res <- score_set(rk, gs)
  cat(sprintf("  family_id=%2d  %-32s  n_LE=%4d  in_WTA=%4d\n",
              p$family_id, p$short, length(gs), res$n_in))
  rows[[length(rows) + 1]] <- data.table(
    signature_type = "family",
    signature_id   = sprintf("F%02d", p$family_id),
    signature_lab  = p$short,
    family         = sprintf("F%d", p$family_id),
    evidence_role  = p$role,
    chamber_driver = p$chamber_driver,
    n_genes_total  = length(gs),
    n_genes_in_WTA = res$n_in,
    gsm            = md_in$gsm,
    title          = md_in$title,
    score          = res$score)
}

cat("\n--- Cell-type panels ---\n")
for (n in names(CELL_PANELS)) {
  gs <- CELL_PANELS[[n]]
  res <- score_set(rk, gs)
  cat(sprintf("  %-22s  n_total=%2d  in_WTA=%2d\n",
              n, length(gs), res$n_in))
  rows[[length(rows) + 1]] <- data.table(
    signature_type = "cell_type",
    signature_id   = n,
    signature_lab  = n,
    family         = NA_character_,
    evidence_role  = NA_character_,
    chamber_driver = NA_character_,
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
long[, eligible_for_5d := (tissue == "Atrial appendage") &
                          (marker_norm == "cardiac") &
                          (!is.na(rhythm)) &
                          (area %in% c("distal", "border"))]

fwrite(long, file.path(G, "roi_scores_v3.tsv"), sep = "\t")
cat(sprintf("\nwrote roi_scores_v3.tsv  (%d rows)\n", nrow(long)))

wide <- dcast(long[, .(gsm, signature_lab, score)],
              gsm ~ signature_lab, value.var = "score")
fwrite(wide, file.path(G, "roi_scores_v3_wide.tsv"), sep = "\t")
cat(sprintf("wrote roi_scores_v3_wide.tsv  (%d rows x %d cols)\n",
            nrow(wide), ncol(wide)))

cat("\n=== AA cardiac ROI (n = 24) AF - SR Delta summary per family ===\n")
elig <- long[eligible_for_5d == TRUE & signature_type == "family"]
af_sr <- elig[, .(
    score_af  = mean(score[rhythm == "AF"], na.rm = TRUE),
    score_sr  = mean(score[rhythm == "SR"], na.rm = TRUE),
    delta_af_minus_sr = mean(score[rhythm == "AF"], na.rm = TRUE)
                        - mean(score[rhythm == "SR"], na.rm = TRUE),
    n_af = sum(rhythm == "AF" & !is.na(score)),
    n_sr = sum(rhythm == "SR" & !is.na(score))
  ), by = .(signature_lab)]
print(af_sr[order(-abs(delta_af_minus_sr))])

cat("\nG3b v3 DONE.\n")
