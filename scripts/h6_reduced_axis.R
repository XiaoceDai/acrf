#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# H6 — Reduced-axis readout (Figure S9).
#
# Goal: show that the dual-axis substrate framework can be approximately
#       reproduced by a small (5-gene) reduced panel per axis.  Not a
#       classifier, no ROC; framework-feasibility check only.
#
# Locked panels (drawn from already-used component genes in Figs 3E / 6C
# / Fig 1 OXPHOS LE):
#
#   matrix_reduced    = PRELP, OMD, FMOD, CHST3, B4GALT1
#   metabolic_reduced = IDH2, NDUFA8, NDUFS2, COX5A, UQCRC1
#
# Outputs (under results/):
#   figure_s9_reduced_panels.tsv          locked gene panels
#   figure_s9_reduced_scores_long.tsv     donor x cohort x axis (full / reduced)
#   figure_s9_reduced_correlation.tsv     per-cohort correlation full vs reduced
#   figure_s9_reduced_centroids.tsv       cohort-level SR / AF centroids on reduced
#   figure_s9_reduced_duration.tsv        LAA duration ρ for full vs reduced metab
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(singscore)
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
P1 <- file.path(ROOT, "Project_P1")
RESULTS <- file.path(ROOT, "results"); dir.create(RESULTS, showWarnings = FALSE)
COHORTS <- c("GSE128188", "GSE79768", "GSE115574")

# ---------------------------------------------------------------------------
# Locked reduced panels
# ---------------------------------------------------------------------------
PANELS <- list(
  matrix    = c("PRELP", "OMD", "FMOD", "CHST3", "B4GALT1"),
  metabolic = c("IDH2", "NDUFA8", "NDUFS2", "COX5A", "UQCRC1")
)
PANEL_DIR <- list(matrix = -1L, metabolic = +1L)  # baseline_dir per axis

# Persist panel definitions
fwrite(rbindlist(lapply(names(PANELS), function(ax) {
  data.table(axis = ax, gene = PANELS[[ax]], baseline_dir = PANEL_DIR[[ax]])
})), file.path(RESULTS, "figure_s9_reduced_panels.tsv"), sep = "\t")

# ---------------------------------------------------------------------------
# Load expression matrices + healthy bulk
# ---------------------------------------------------------------------------
load_mat <- function(expr_path) {
  expr <- fread(expr_path, sep = "\t", header = TRUE)
  sym <- expr[[1]]; M <- as.matrix(expr[, -1, with = FALSE])
  rownames(M) <- sym
  if (anyDuplicated(sym)) {
    v <- apply(M, 1, var, na.rm = TRUE)
    keep <- tapply(seq_along(sym), sym, function(idx) idx[which.max(v[idx])])
    M <- M[unlist(keep, use.names = FALSE), , drop = FALSE]
  }
  M
}

cat("loading matrices ...\n")
mats <- list()
for (cohort in COHORTS) {
  Mp <- load_mat(file.path(P1, cohort, "expr_symbol.tsv"))
  meta_p <- fread(file.path(P1, cohort, "sample_metadata.tsv"))
  mats[[cohort]] <- list(rk = rankGenes(Mp), meta = meta_p)
  cat(sprintf("  %s: %d genes x %d samples\n", cohort,
              nrow(Mp), ncol(Mp)))
}
# healthy bulk for baseline magnitude
bulk_expr <- read.delim(file.path(ROOT, "bulk", "logcpm_voom.tsv"),
                        check.names = FALSE)
bulk_meta <- read.delim(file.path(ROOT, "bulk", "sample_meta.tsv"))
for (col in c("sample_id", "chamber", "donor"))
  if (col %in% colnames(bulk_meta))
    bulk_meta[[col]] <- gsub('^"|"$', '', bulk_meta[[col]])
sym <- bulk_expr$gene_symbol
samp_cols <- setdiff(colnames(bulk_expr), c("gene_key", "gene_symbol"))
Mb <- as.matrix(bulk_expr[, samp_cols, drop = FALSE])
rownames(Mb) <- sym
if (anyDuplicated(sym)) {
  v <- apply(Mb, 1, var, na.rm = TRUE)
  keep <- tapply(seq_along(sym), sym, function(idx) idx[which.max(v[idx])])
  Mb <- Mb[unlist(keep, use.names = FALSE), , drop = FALSE]
}
mats[["GSE112339"]] <- list(rk = rankGenes(Mb), meta = bulk_meta)

# ---------------------------------------------------------------------------
# Score one panel via singscore + per-donor aCAI (paired cohort)
# ---------------------------------------------------------------------------
score_panel <- function(panel_genes, ranks_obj) {
  hit <- intersect(panel_genes, rownames(ranks_obj$rk))
  if (length(hit) < 3L) return(rep(NA_real_, ncol(ranks_obj$rk)))
  s <- simpleScore(ranks_obj$rk, upSet = hit,
                   knownDirection = TRUE, centerScore = TRUE)
  setNames(s$TotalScore, rownames(s))
}

acai_per_donor <- function(scores, meta, baseline_dir) {
  meta <- as.data.table(meta)
  meta[, score := scores[match(meta$gsm, names(scores))]]
  by_donor <- meta[!is.na(score), .(
    LA = mean(score[chamber %in% c("LA", "left atrium", "left_atrium")],
              na.rm = TRUE),
    RA = mean(score[chamber %in% c("RA", "right atrium", "right_atrium")],
              na.rm = TRUE),
    rhythm = unique(rhythm[1])
  ), by = donor_id]
  by_donor <- by_donor[!is.na(LA) & !is.na(RA)]
  by_donor[, CAI := LA - RA]
  by_donor[, aCAI := baseline_dir * CAI]
  by_donor
}

bulk_mag <- function(scores_bulk, meta_bulk) {
  meta_bulk <- as.data.table(meta_bulk)
  meta_bulk[, score := scores_bulk[match(sample_id, names(scores_bulk))]]
  by_d <- meta_bulk[!is.na(score), .(
    LA = mean(score[chamber == "LA"], na.rm = TRUE),
    RA = mean(score[chamber == "RA"], na.rm = TRUE)
  ), by = donor]
  abs(mean(by_d$LA - by_d$RA, na.rm = TRUE))
}

# ---------------------------------------------------------------------------
# Compute reduced-axis aCAI_rel per donor per cohort
# ---------------------------------------------------------------------------
all_scores <- list()
for (ax in names(PANELS)) {
  bd <- PANEL_DIR[[ax]]
  s_bulk <- score_panel(PANELS[[ax]], mats[["GSE112339"]])
  bm <- bulk_mag(s_bulk, mats[["GSE112339"]]$meta)
  for (cohort in COHORTS) {
    sc <- score_panel(PANELS[[ax]], mats[[cohort]])
    bd_per <- acai_per_donor(sc, mats[[cohort]]$meta, bd)
    bd_per[, axis := ax]
    bd_per[, cohort := cohort]
    bd_per[, aCAI_rel_reduced := aCAI / bm]
    all_scores[[paste(ax, cohort, sep = "_")]] <- bd_per
  }
}
red <- rbindlist(all_scores, use.names = TRUE)
setnames(red, "donor_id", "donor")

# ---------------------------------------------------------------------------
# Pull full-axis composite from Fig 6 results
# ---------------------------------------------------------------------------
full <- fread(file.path(RESULTS, "figure6_axis_composite_long.tsv"))
full <- full[!is.na(composite),
             .(donor = donor_id, cohort, rhythm,
               axis = ifelse(axis == "STRUCTURAL", "matrix",
                      ifelse(axis == "METABOLIC", "metabolic", axis)),
               aCAI_rel_full = composite)]
# Some donor IDs in fig6 axis composite may not match (e.g. 'donor_id' types).
# Use string matching directly.
mrg <- merge(red[, .(donor, cohort, axis, rhythm, aCAI_rel_reduced)],
             full, by = c("donor", "cohort", "axis", "rhythm"),
             all = FALSE)
mrg <- mrg[!is.na(aCAI_rel_full) & !is.na(aCAI_rel_reduced)]
fwrite(mrg, file.path(RESULTS, "figure_s9_reduced_scores_long.tsv"), sep = "\t")
cat(sprintf("\nMerged reduced+full rows: %d\n", nrow(mrg)))

# ---------------------------------------------------------------------------
# Correlation: reduced vs full (per axis × cohort)
# ---------------------------------------------------------------------------
cor_tbl <- mrg[, {
  if (.N < 4) return(list(n = .N, spearman_rho = NA_real_,
                          pearson_r = NA_real_))
  r_s <- suppressWarnings(cor(aCAI_rel_reduced, aCAI_rel_full, method = "spearman"))
  r_p <- suppressWarnings(cor(aCAI_rel_reduced, aCAI_rel_full, method = "pearson"))
  list(n = .N, spearman_rho = r_s, pearson_r = r_p)
}, by = .(axis, cohort)]
fwrite(cor_tbl, file.path(RESULTS, "figure_s9_reduced_correlation.tsv"),
       sep = "\t")
cat("\n--- correlation reduced vs full ---\n")
print(cor_tbl)

# ---------------------------------------------------------------------------
# Centroid shift on reduced axis (per cohort, SR vs AF)
# ---------------------------------------------------------------------------
cent_red <- red[!is.na(aCAI_rel_reduced),
                .(mean = mean(aCAI_rel_reduced)),
                by = .(axis, cohort, rhythm)]
cent_red <- dcast(cent_red, axis + cohort ~ rhythm, value.var = "mean")
cent_red[, delta_AF_minus_SR := AF - SR]
fwrite(cent_red, file.path(RESULTS, "figure_s9_reduced_centroids.tsv"),
       sep = "\t")
cat("\n--- reduced centroids ---\n"); print(cent_red)

# ---------------------------------------------------------------------------
# LAA duration ρ on reduced metabolic panel (vs full Fig 4 OXPHOS / METAB)
# ---------------------------------------------------------------------------
laa_path  <- file.path(ROOT, "Project_P4", "GSE41177", "expr_symbol_LAA.tsv")
meta_path <- file.path(ROOT, "Project_P4", "GSE41177", "sample_metadata.tsv")
if (file.exists(laa_path) && file.exists(meta_path)) {
  M_laa <- load_mat(laa_path)
  rk <- rankGenes(M_laa)
  meta_laa <- fread(meta_path)
  rho_rows <- list()
  for (ax in names(PANELS)) {
    sc <- score_panel(PANELS[[ax]], list(rk = rk, meta = meta_laa))
    md <- meta_laa[tissue == "LAA" & rhythm == "AF" &
                   !is.na(duration_months) & keep_for_primary == TRUE]
    md[, score := sc[match(gsm, names(sc))]]
    md <- md[!is.na(score)]
    if (nrow(md) >= 6) {
      r <- suppressWarnings(cor(md$score, md$duration_months,
                                 method = "spearman"))
      rho_rows[[ax]] <- data.table(axis = ax, n_AF = nrow(md),
                                   reduced_rho_LAA = r)
    }
  }
  rho_red <- rbindlist(rho_rows, use.names = TRUE)
  # Reference: real OXPHOS full ρ from Fig 4 primary table
  full_dur <- fread(file.path(ROOT, "Project_P4", "GSE41177", "results",
                              "duration_trend_laa_primary.tsv"))
  full_dur <- full_dur[, .(program_label, full_rho_LAA = statistic)]
  rho_full_ox <- full_dur[program_label == "OXPHOS", full_rho_LAA]
  rho_red[axis == "metabolic", full_rho_LAA_reference :=
            rho_full_ox]
  fwrite(rho_red, file.path(RESULTS, "figure_s9_reduced_duration.tsv"),
         sep = "\t")
  cat("\n--- reduced LAA duration ρ ---\n"); print(rho_red)
}

cat("\nH6 DONE.\n")
