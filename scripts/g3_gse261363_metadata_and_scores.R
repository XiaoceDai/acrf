#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# G3 — Minimal Figure 5D inset:
#       1) ROI metadata (atrial appendage vs EAT, distal vs border, AF / SR)
#       2) Two ROI-level singscores (frozen):
#            - ECM proteoglycans     (LE genes from normal_chamber_programs.tsv)
#            - canonical fibroblast  (small explicit panel, NOT a new signature)
#       3) Minimal QC: distributions of the two scores.
#
# Subordinate-only (FIGURE5_WORKFLOW.md §5D):
#   - no causal claim
#   - no new cell-type inference
#   - no new program panel
#   - no AF vs SR primary statistic
#
# Inputs:
#   Project_P5/GSE261363/GSE261363_All_Data_WTA_1_7_22.xlsx
#   Project_P5/GSE261363/GSE261363_series_matrix.txt.gz
#   normal_chamber_programs.tsv
#
# Outputs:
#   Project_P5/GSE261363/roi_metadata.tsv
#   Project_P5/GSE261363/roi_scores.tsv
#   Project_P5/GSE261363/qc/g3_score_distributions.pdf
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
QC   <- file.path(G, "qc")
dir.create(QC, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Series-matrix metadata ----
parse_series <- function(path) {
  L <- readLines(gzfile(path), warn = FALSE)
  pick <- function(prefix) {
    rows <- grep(prefix, L, value = TRUE, fixed = TRUE)
    rows
  }
  parse_row <- function(line) {
    toks <- strsplit(line, "\t", fixed = TRUE)[[1]]
    gsub('^"|"$', '', toks[-1])
  }
  out <- list(
    title  = parse_row(pick("!Sample_title")[1]),
    gsm    = parse_row(pick("!Sample_geo_accession")[1]),
    source = parse_row(pick("!Sample_source_name_ch1")[1])
  )
  ch <- pick("!Sample_characteristics_ch1")
  for (k in ch) {
    vals <- parse_row(k)
    tag  <- sub("^\\s*([^:]+):.*$", "\\1", vals[1])
    val  <- sub("^[^:]+:\\s*", "", vals)
    out[[tag]] <- val
  }
  as.data.table(out)
}
md <- parse_series(file.path(G, "GSE261363_series_matrix.txt.gz"))
setnames(md, c("title", "gsm", "source", "tissue", "slide_id",
               "roi_number", "area"))
# control well = first sample (slide_id == "control")
md[, is_control := slide_id == "control"]
md[, rhythm := fifelse(grepl("AF$", slide_id), "AF",
                fifelse(grepl("SR$", slide_id), "SR", NA_character_))]
md[, patient := sub("-(AF|SR)$", "", slide_id)]
# segment marker (DNA+ / FABP4+) is encoded in title after the second '|'
md[, marker := vapply(strsplit(title, "|", fixed = TRUE),
                      function(x) trimws(x[length(x)]), character(1))]

cat("ROI metadata (incl. control):\n")
print(md[, .N, by = .(rhythm, tissue, area, marker)][order(rhythm, tissue, area, marker)])

# ---- 2. Q3-normalized expression matrix from xlsx ----
xls <- file.path(G, "GSE261363_All_Data_WTA_1_7_22.xlsx")
sp  <- as.data.table(read_excel(xls, sheet = "SegmentProperties"))
expr <- as.data.table(read_excel(xls, sheet = "Q3_5percfilter"))

# match column names (e.g. "P161_P167|001|DNA+") to GSM via SegmentProperties → series_matrix title
# SegmentProperties has SegmentDisplayName matching expr column names (without leading "X").
# series_matrix has Sample_title that should match SegmentDisplayName.
title_map <- md[, .(title, gsm)]
expr_cols <- setdiff(names(expr), "TargetName")
hit <- match(expr_cols, title_map$title)
cat(sprintf("\nWTA matrix: %d genes x %d ROIs (control excluded if present)\n",
            nrow(expr), length(expr_cols)))
cat(sprintf("Matched %d / %d expr columns to series_matrix titles\n",
            sum(!is.na(hit)), length(expr_cols)))

expr_mat <- as.matrix(expr[, ..expr_cols])
rownames(expr_mat) <- expr$TargetName

# Restrict metadata to ROIs in matrix
md_in <- md[match(expr_cols, md$title)]
stopifnot(all(md_in$title == expr_cols))

# ---- 3. Frozen gene sets ----
prog <- read.delim(file.path(ROOT, "normal_chamber_programs.tsv"),
                   sep = "\t", stringsAsFactors = FALSE)
for (col in c("program_id", "leading_edge_genes"))
  if (col %in% names(prog)) prog[[col]] <- gsub('^"|"$', '', prog[[col]])
ecm_id <- "REACTOME__REACTOME_DISEASES_ASSOCIATED_WITH_GLYCOSAMINOGLYCAN_METABOLISM"
ecm_le <- unique(trimws(strsplit(
  prog$leading_edge_genes[prog$program_id == ecm_id], ",", fixed = TRUE)[[1]]))

# Canonical fibroblast marker panel — small, declared, NOT data-derived
FIBRO_MARKERS <- c("PDGFRA", "COL1A1", "COL1A2", "COL3A1", "COL5A1",
                   "DCN", "LUM", "POSTN", "FAP", "FBN1", "MMP2", "BGN")

cat(sprintf("\nECM proteoglycans LE genes: %d (matched in WTA: %d)\n",
            length(ecm_le), sum(ecm_le %in% rownames(expr_mat))))
cat(sprintf("Canonical fibroblast markers: %d (matched in WTA: %d)\n",
            length(FIBRO_MARKERS),
            sum(FIBRO_MARKERS %in% rownames(expr_mat))))

# ---- 4. ROI-level singscores ----
# Q3 values are already on log-ish positive scale; rank-based scoring is
# scale-invariant within column.
rk <- rankGenes(expr_mat)
score_set <- function(rk, gs, label) {
  gs_in <- intersect(gs, rownames(rk))
  if (length(gs_in) < 3L) {
    return(data.table(score = rep(NA_real_, ncol(rk)),
                      n_in = length(gs_in)))
  }
  s <- simpleScore(rk, upSet = gs_in, knownDirection = TRUE,
                   centerScore = TRUE)
  data.table(score = s$TotalScore, n_in = length(gs_in))
}
ecm_sc   <- score_set(rk, ecm_le, "ECM_proteoglycans")
fibro_sc <- score_set(rk, FIBRO_MARKERS, "fibroblast_canonical")

scores <- data.table(
  gsm = md_in$gsm,
  title = md_in$title,
  patient = md_in$patient,
  rhythm = md_in$rhythm,
  tissue = md_in$tissue,
  area = md_in$area,
  marker = md_in$marker,
  ecm_score = ecm_sc$score,
  fibro_score = fibro_sc$score,
  n_ecm_genes_in = ecm_sc$n_in,
  n_fibro_genes_in = fibro_sc$n_in
)

# ---- 5. Mark the 5D candidate ROI set:
#   atrial appendage tissue (AA), DNA+ marker, AF or SR slide,
#   area in {distal, border}.  These are the ROIs we will use for the
#   ECM × fibroblast sanity scatter.
# SR slides label the cardiac segment as "DNA" (no "+"); AF slides as "DNA+".
# Treat both as the cardiac/atrial segment for the 5D inset eligibility.
scores[, marker_norm := fifelse(marker %in% c("DNA", "DNA+"), "cardiac",
                         fifelse(marker == "FABP4+", "adipose",
                                 marker))]
scores[, eligible_for_5D := (tissue == "Atrial appendage") &
                            (marker_norm == "cardiac") &
                            (!is.na(rhythm)) &
                            (area %in% c("distal", "border"))]

fwrite(scores,  file.path(G, "roi_scores.tsv"),  sep = "\t")
fwrite(md,      file.path(G, "roi_metadata.tsv"), sep = "\t")
cat("\nwrote roi_metadata.tsv  (", nrow(md), " rows; ",
    sum(!md$is_control), " non-control ROIs)\n", sep = "")
cat("wrote roi_scores.tsv   (", nrow(scores), " rows; ",
    sum(scores$eligible_for_5D), " eligible_for_5D)\n", sep = "")

cat("\n--- ROI eligibility breakdown ---\n")
print(scores[, .N, by = .(tissue, area, marker, rhythm)][order(tissue, area, marker, rhythm)])

cat("\n--- Score summary (non-NA) ---\n")
print(scores[!is.na(ecm_score),
             .(n = .N,
               ecm_min = round(min(ecm_score), 3),
               ecm_max = round(max(ecm_score), 3),
               fibro_min = round(min(fibro_score), 3),
               fibro_max = round(max(fibro_score), 3))])

cat("\n--- Score distribution by group (eligible only) ---\n")
elig <- scores[eligible_for_5D == TRUE]
print(elig[, .(n = .N,
               ecm_med = round(median(ecm_score), 3),
               fibro_med = round(median(fibro_score), 3)),
           by = .(rhythm, area)][order(rhythm, area)])

# ---- 6. QC distribution PDF ----
pdf(file.path(QC, "g3_score_distributions.pdf"), width = 6.5, height = 5.0)
op <- par(mfrow = c(2, 2), mar = c(3.5, 3.5, 2.4, 0.5),
          mgp = c(2, 0.6, 0), cex = 0.8, las = 1)

# (1) histogram by tissue
breaks <- pretty(scores$ecm_score, n = 18)
hist(scores$ecm_score[scores$tissue == "Atrial appendage"],
     breaks = breaks, col = "#7570B355", border = "white",
     main = "ECM proteoglycans score", xlab = "singscore", ylab = "ROIs")
hist(scores$ecm_score[scores$tissue == "Epicardial AT"],
     breaks = breaks, col = "#E6AB0255", border = "white", add = TRUE)
legend("topleft", c("Atrial appendage", "Epicardial AT"),
       fill = c("#7570B388", "#E6AB0288"), border = NA, bty = "n", cex = 0.8)

breaks <- pretty(scores$fibro_score, n = 18)
hist(scores$fibro_score[scores$tissue == "Atrial appendage"],
     breaks = breaks, col = "#7570B355", border = "white",
     main = "Fibroblast canonical score", xlab = "singscore", ylab = "ROIs")
hist(scores$fibro_score[scores$tissue == "Epicardial AT"],
     breaks = breaks, col = "#E6AB0255", border = "white", add = TRUE)

# (3) ECM vs fibroblast scatter — eligible ROIs only
cols <- ifelse(elig$rhythm == "AF", "#B35806", "#542788")
pchs <- ifelse(elig$area  == "border", 17, 16)
plot(elig$fibro_score, elig$ecm_score,
     col = cols, pch = pchs, cex = 1.4,
     xlab = "fibroblast (canonical) score",
     ylab = "ECM proteoglycans score",
     main = "AA ROIs (DNA+) — fibroblast vs ECM")
abline(h = 0, v = 0, col = "#bbbbbb", lty = 3)
legend("topleft", c("AF", "SR"), col = c("#B35806", "#542788"),
       pch = 16, bty = "n", cex = 0.8)
legend("bottomright", c("border", "distal"), pch = c(17, 16),
       col = "#444444", bty = "n", cex = 0.8)
rho <- suppressWarnings(cor(elig$fibro_score, elig$ecm_score,
                            method = "spearman", use = "complete.obs"))
mtext(sprintf("Spearman ρ = %.2f  (n=%d, descriptive)", rho, nrow(elig)),
      side = 3, line = -1.4, cex = 0.7, adj = 0.02, col = "#444444")

# (4) all ROIs, marker-coloured to confirm cardiac vs adipose separation
cols2 <- ifelse(scores$marker_norm == "cardiac", "#1B7837",
         ifelse(scores$marker_norm == "adipose", "#762A83", "#aaaaaa"))
plot(scores$fibro_score, scores$ecm_score,
     col = cols2, pch = 16, cex = 1.1,
     xlab = "fibroblast (canonical) score",
     ylab = "ECM proteoglycans score",
     main = "All ROIs - separation by GeoMx segment marker")
legend("topleft", c("cardiac (DNA / DNA+)", "adipose (FABP4+)"),
       col = c("#1B7837", "#762A83"), pch = 16, bty = "n", cex = 0.8)
abline(h = 0, v = 0, col = "#bbbbbb", lty = 3)

par(op); dev.off()
cat("wrote ", file.path(QC, "g3_score_distributions.pdf"), "\n", sep = "")

cat("\nG3 DONE.\n")
