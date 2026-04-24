#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# F5 — TMM + logCPM (per compartment) + singscore for the 6 frozen support
#       programs from Figure 5 (3 STRUCTURAL + 3 METABOLIC).
#
# Frozen rules (FIGURE5_WORKFLOW.md):
#   - main compartments only: atrialCM, fibroblast, endothelial_or_vascular, immune
#   - donor-level outputs (no cell-level)
#   - this script: program score only, NO regulator / ATAC layer
#
# Inputs:
#   Project_P5/GSE238242/pseudobulk/gse238242_pseudobulk_counts_main.rds
#   Project_P5/GSE238242/pseudobulk/gse238242_pseudobulk_meta_main.tsv
#   results/figure5_program_sets.tsv          (frozen 6 program IDs)
#   normal_chamber_programs.tsv               (Fig 1 prior; LE genes)
#
# Outputs (Project_P5/GSE238242/results/):
#   figure5_program_scores_pseudobulk.tsv     long: donor x compartment x program
#   figure5_program_delta_main.tsv            mean_AF / mean_SR / delta_score
#   figure5_program_support_check.tsv         pass/fail vs Fig 2-4 expectation
#   figure5_program_coverage.tsv              per-(program, compartment) LE coverage
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(edgeR); library(singscore)
  library(GSEABase)
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
PB   <- file.path(ROOT, "Project_P5", "GSE238242", "pseudobulk")
OUT  <- file.path(ROOT, "Project_P5", "GSE238242", "results")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MAIN_CTS <- c("atrialCM", "fibroblast", "endothelial_or_vascular", "immune")

# --------------------------------------------------------------------
# 1. Frozen program panel (6 programs) + LE gene sets from Fig 1 prior
# --------------------------------------------------------------------
panel_in <- fread(file.path(ROOT, "results", "figure5_program_sets.tsv"))
panel <- panel_in[set_name %in% c("STRUCTURAL", "METABOLIC")]
stopifnot(nrow(panel) == 6L)

prog_in <- read.delim(file.path(ROOT, "normal_chamber_programs.tsv"),
                      sep = "\t", stringsAsFactors = FALSE)
for (col in c("program_id", "program", "direction", "leading_edge_genes")) {
  if (col %in% colnames(prog_in)) prog_in[[col]] <- gsub('^"|"$', '', prog_in[[col]])
}

idx <- match(panel$program_id, prog_in$program_id)
stopifnot(!anyNA(idx))
panel[, baseline_direction := ifelse(prog_in$direction[idx] == "LA-biased", 1L, -1L)]
panel[, le_str := prog_in$leading_edge_genes[idx]]
panel[, le_genes := lapply(le_str,
                           function(s) unique(trimws(strsplit(s, ",", fixed = TRUE)[[1]])))]
panel[, n_LE_total := vapply(le_genes, length, integer(1))]

cat("--- Figure 5 frozen panel ---\n")
print(panel[, .(set_name, program_id, program_label,
                baseline_direction, n_LE_total)])

# --------------------------------------------------------------------
# 2. Pseudobulks + meta
# --------------------------------------------------------------------
M    <- as.matrix(readRDS(file.path(PB, "gse238242_pseudobulk_counts_main.rds")))
meta <- fread(file.path(PB, "gse238242_pseudobulk_meta_main.tsv"))
stopifnot(all(meta$pb_id %in% colnames(M)))
M <- M[, meta$pb_id]
stopifnot(all(meta$compartment %in% MAIN_CTS))
cat(sprintf("\nloaded %d genes x %d pb (%d donors)\n",
            nrow(M), ncol(M), uniqueN(meta$donor)))

# --------------------------------------------------------------------
# 3. Per-compartment TMM + logCPM, then singscore each program
# --------------------------------------------------------------------
score_one <- function(M_lcpm, le, prog_id, prog_label, compartment) {
  rk    <- rankGenes(M_lcpm)                       # genes x samples ranks
  up_in <- intersect(le, rownames(M_lcpm))
  if (length(up_in) < 3L) {
    return(data.table(
      compartment = compartment,
      program_id  = prog_id,
      program_label = prog_label,
      n_LE_matched = length(up_in),
      gsm_or_donor = colnames(M_lcpm),
      score = NA_real_
    ))
  }
  s <- simpleScore(rk, upSet = up_in, knownDirection = TRUE,
                   centerScore = TRUE)
  data.table(
    compartment = compartment,
    program_id  = prog_id,
    program_label = prog_label,
    n_LE_matched = length(up_in),
    gsm_or_donor = rownames(s),
    score = s$TotalScore
  )
}

scores_long <- list()
coverage    <- list()
for (cmp in MAIN_CTS) {
  m_meta <- meta[compartment == cmp]
  m_pb   <- M[, m_meta$pb_id, drop = FALSE]
  cat(sprintf("\n[%s] %d pb ... ", cmp, ncol(m_pb)))

  y <- DGEList(counts = m_pb)
  keep <- filterByExpr(y, group = m_meta$rhythm, min.count = 5)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  lcpm <- cpm(y, log = TRUE, prior.count = 1)
  cat(sprintf("genes after filter = %d\n", nrow(lcpm)))

  for (i in seq_len(nrow(panel))) {
    res <- score_one(lcpm, panel$le_genes[[i]],
                     panel$program_id[i], panel$program_label[i], cmp)
    # join donor / rhythm using pb_id == gsm_or_donor here
    res2 <- merge(res, m_meta[, .(pb_id, donor, rhythm, n_cells)],
                  by.x = "gsm_or_donor", by.y = "pb_id", sort = FALSE)
    setnames(res2, "gsm_or_donor", "pb_id")
    scores_long[[paste(cmp, panel$program_id[i], sep = "::")]] <- res2

    coverage[[paste(cmp, panel$program_id[i], sep = "::")]] <- data.table(
      compartment = cmp, program_id = panel$program_id[i],
      program_label = panel$program_label[i],
      n_LE_total = panel$n_LE_total[i],
      n_LE_matched = unique(res$n_LE_matched)
    )
  }
}
sc <- rbindlist(scores_long)
cv <- rbindlist(coverage)
sc <- sc[, .(donor, rhythm, compartment, program_id, program_label,
             score, n_cells, n_LE_matched)]

fwrite(sc, file.path(OUT, "figure5_program_scores_pseudobulk.tsv"), sep = "\t")
fwrite(cv, file.path(OUT, "figure5_program_coverage.tsv"), sep = "\t")
cat(sprintf("\nwrote figure5_program_scores_pseudobulk.tsv  (%d rows)\n", nrow(sc)))

# --------------------------------------------------------------------
# 4. Δscore (AF − SR) per (program, compartment)
# --------------------------------------------------------------------
delta <- sc[, .(
  n_AF      = uniqueN(donor[rhythm == "AF"]),
  n_SR      = uniqueN(donor[rhythm == "SR"]),
  mean_AF   = mean(score[rhythm == "AF"], na.rm = TRUE),
  mean_SR   = mean(score[rhythm == "SR"], na.rm = TRUE),
  sd_AF     = sd(score[rhythm == "AF"], na.rm = TRUE),
  sd_SR     = sd(score[rhythm == "SR"], na.rm = TRUE)
), by = .(program_id, program_label, compartment)]
delta[, delta_score := mean_AF - mean_SR]

# Attach axis (STRUCTURAL/METABOLIC) and locked support compartment
panel_axis <- panel[, .(program_id, axis = set_name)]
delta <- merge(delta, panel_axis, by = "program_id", sort = FALSE)
setorder(delta, axis, program_label, compartment)

fwrite(delta, file.path(OUT, "figure5_program_delta_main.tsv"), sep = "\t")
cat(sprintf("wrote figure5_program_delta_main.tsv  (%d rows)\n", nrow(delta)))

# --------------------------------------------------------------------
# 5. Support check vs Fig 2-4 expectation
#    (locked: STRUCTURAL programs must lift in FIBROBLAST;
#             METABOLIC programs must drop in ATRIALCM)
# --------------------------------------------------------------------
EXPECT <- list(STRUCTURAL = list(compartment = "fibroblast", sign = +1L),
               METABOLIC  = list(compartment = "atrialCM",   sign = -1L))

# Thresholds (singscore is centred; ranges typically |Δ| < 0.10).
# strong   |Δ| ≥ 0.030 in expected direction
# moderate |Δ| ∈ [0.010, 0.030) in expected direction
# weak     |Δ| < 0.010
# absent   |Δ| ≥ 0.010 in opposite direction
classify <- function(d, expect_sign) {
  if (is.na(d)) return("NA")
  ad <- abs(d); s <- sign(d)
  if (ad < 0.010) return("weak")
  if (s == expect_sign) {
    if (ad >= 0.030) "strong" else "moderate"
  } else {
    "absent"
  }
}

support_rows <- list()
for (ax_name in names(EXPECT)) {
  cmp_exp <- EXPECT[[ax_name]]$compartment
  sgn_exp <- EXPECT[[ax_name]]$sign
  axis_target <- ax_name
  sub <- delta[axis == axis_target & compartment == cmp_exp]
  for (i in seq_len(nrow(sub))) {
    d   <- sub$delta_score[i]
    lvl <- classify(d, sgn_exp)
    support_rows[[paste(ax_name, sub$program_id[i], sep = "::")]] <- data.table(
      axis = ax_name,
      program_id = sub$program_id[i],
      program_label = sub$program_label[i],
      compartment_expected = cmp_exp,
      expected_sign = sgn_exp,
      delta_observed = d,
      sign_matches_expected = (sign(d) == sgn_exp) & abs(d) >= 0.010,
      support_level = lvl
    )
  }
}
support <- rbindlist(support_rows)
fwrite(support, file.path(OUT, "figure5_program_support_check.tsv"), sep = "\t")
cat(sprintf("wrote figure5_program_support_check.tsv  (%d rows)\n", nrow(support)))

cat("\n--- support check ---\n")
print(support)

cat("\nF5 DONE.\n")
