#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# H10 — Sex-stratified sensitivity for the dual-axis chamber-rewiring meta.
#
# For each cohort that has BOTH male and female donors with BOTH SR and AF,
# we re-run Fig 6A axis composite + Hedges' g + REML / Hartung–Knapp pool
# on the male-only and female-only subsets, then compare to the
# all-donor pooled g.
#
# Cohorts with reliable per-donor sex:
#   GSE128188 — all male; no sex stratification
#   GSE79768  — M=5 / F=8; can stratify
#   GSE115574 — sex NA in GEO; cannot stratify
#   GSE41177  — M=8 / F=8 (LAA only)
#
# Output:
#   results/sex_sensitivity_axis_long.tsv
#   results/sex_sensitivity_axis_pooled.tsv
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(metafor)
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
RES  <- file.path(ROOT, "results")

hedges_g <- function(x_AF, x_SR) {
  nA <- length(x_AF); nS <- length(x_SR)
  if (nA < 2 || nS < 2) return(list(g = NA_real_, var = NA_real_,
                                    n_AF = nA, n_SR = nS))
  m_A <- mean(x_AF); s_A <- sd(x_AF)
  m_S <- mean(x_SR); s_S <- sd(x_SR)
  df <- nA + nS - 2
  s_p <- sqrt(((nA - 1) * s_A^2 + (nS - 1) * s_S^2) / df)
  if (!is.finite(s_p) || s_p == 0)
    return(list(g = NA_real_, var = NA_real_, n_AF = nA, n_SR = nS))
  d <- (m_A - m_S) / s_p
  J <- 1 - 3 / (4 * df - 1)
  g <- J * d
  v <- (nA + nS) / (nA * nS) + g^2 / (2 * (nA + nS))
  list(g = g, var = v, n_AF = nA, n_SR = nS)
}

reml_pool <- function(g_vec, v_vec) {
  ok <- !is.na(g_vec) & !is.na(v_vec) & v_vec > 0
  if (sum(ok) < 2) return(list(g = NA_real_, lo = NA_real_, hi = NA_real_,
                               I2 = NA_real_, p = NA_real_))
  fit <- tryCatch(rma(yi = g_vec[ok], vi = v_vec[ok],
                      method = "REML", test = "knha"),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(g = NA_real_, lo = NA_real_, hi = NA_real_,
                                I2 = NA_real_, p = NA_real_))
  list(g = as.numeric(fit$beta), lo = as.numeric(fit$ci.lb),
       hi = as.numeric(fit$ci.ub),
       I2 = as.numeric(fit$I2), p = fit$pval)
}

# ---------------------------------------------------------------------------
# Load axis composite (already donor-level) + per-donor sex from cohorts
# ---------------------------------------------------------------------------
comp <- fread(file.path(RES, "figure6_axis_composite_long.tsv"))

# Per-cohort donor → sex map
read_sex <- function(path) {
  if (!file.exists(path)) return(data.table())
  md <- fread(path)
  if (!"donor_id" %in% names(md)) return(data.table())
  md <- md[, .(donor_id = as.character(donor_id), sex = sex)]
  md <- unique(md, by = "donor_id")
  md
}

sex_map <- list(
  GSE128188 = read_sex(file.path(ROOT, "Project_P1", "GSE128188", "sample_metadata.tsv")),
  GSE79768  = read_sex(file.path(ROOT, "Project_P1", "GSE79768",  "sample_metadata.tsv")),
  GSE115574 = read_sex(file.path(ROOT, "Project_P1", "GSE115574", "sample_metadata.tsv"))
)
for (c in names(sex_map)) {
  if (nrow(sex_map[[c]]) == 0) next
  sex_map[[c]][, cohort := c]
}
sex_long <- rbindlist(sex_map[lengths(sex_map) > 0], fill = TRUE)
sex_long[, sex := tolower(sex)]
sex_long[sex %in% c("male", "m"),    sex_norm := "M"]
sex_long[sex %in% c("female", "f"),  sex_norm := "F"]

cat("Sex coverage by cohort:\n")
print(sex_long[, .N, by = .(cohort, sex_norm)])

# Merge into composite
comp[, donor_id := as.character(donor_id)]
mrg <- merge(comp, sex_long[, .(cohort, donor_id, sex_norm)],
             by = c("cohort", "donor_id"), all.x = TRUE)
mrg <- mrg[!is.na(composite)]

# ---------------------------------------------------------------------------
# Stratified per-cohort g
# ---------------------------------------------------------------------------
strata <- c("ALL", "M", "F")
rec_long <- list()
for (axis_lab in c("STRUCTURAL", "METABOLIC")) {
  for (.coh in unique(mrg$cohort)) {
    for (st in strata) {
      sub <- mrg[axis == axis_lab & cohort == .coh]
      if (st != "ALL") sub <- sub[sex_norm == st]
      g <- hedges_g(sub$composite[sub$rhythm == "AF"],
                    sub$composite[sub$rhythm == "SR"])
      rec_long[[paste(axis_lab, .coh, st, sep = "_")]] <- data.table(
        axis = axis_lab, cohort = .coh, stratum = st,
        n_AF = g$n_AF, n_SR = g$n_SR,
        g = g$g, var_g = g$var,
        se_g = if (is.na(g$var)) NA_real_ else sqrt(g$var)
      )
    }
  }
}
out <- rbindlist(rec_long, use.names = TRUE)
out[, ci_lo := g - 1.96 * se_g]
out[, ci_hi := g + 1.96 * se_g]
fwrite(out, file.path(RES, "sex_sensitivity_axis_long.tsv"), sep = "\t")
cat("\n--- per-(axis, cohort, stratum) Hedges' g ---\n")
print(out[, .(axis, cohort, stratum, n_AF, n_SR,
              g = round(g, 3),
              ci_lo = round(ci_lo, 3),
              ci_hi = round(ci_hi, 3))])

# ---------------------------------------------------------------------------
# Pooled g per (axis, stratum) across cohorts
# ---------------------------------------------------------------------------
pooled_rows <- list()
for (axis_lab in c("STRUCTURAL", "METABOLIC")) {
  for (st in strata) {
    sub <- out[axis == axis_lab & stratum == st]
    pool <- reml_pool(sub$g, sub$var_g)
    pooled_rows[[paste(axis_lab, st, sep = "_")]] <- data.table(
      axis = axis_lab, stratum = st,
      n_studies = sum(!is.na(sub$g) & !is.na(sub$var_g)),
      pooled_g = pool$g, ci_lo = pool$lo, ci_hi = pool$hi,
      I2 = pool$I2, p_HK = pool$p
    )
  }
}
pooled <- rbindlist(pooled_rows, use.names = TRUE)
fwrite(pooled, file.path(RES, "sex_sensitivity_axis_pooled.tsv"), sep = "\t")
cat("\n--- pooled g by (axis, stratum) ---\n")
print(pooled[, .(axis, stratum, n_studies,
                 g = round(pooled_g, 3),
                 ci = sprintf("[%+.2f, %+.2f]", ci_lo, ci_hi),
                 I2 = round(I2, 1),
                 p_HK = ifelse(is.na(p_HK), NA_character_,
                               format.pval(p_HK, 3)))])

cat("\nH10 DONE.\n")
