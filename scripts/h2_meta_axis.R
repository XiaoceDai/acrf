#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# H2 — Per-cohort Hedges' g + REML/Hartung–Knapp random-effects meta on the
#       2 axis composites computed by H1.
#
# Axis-level REML / Hartung–Knapp pooling on cohort g estimates.
#   - g = J(df) * (mean_AF - mean_SR) / s_pooled        (AF - SR convention)
#   - Var(g) = (nAF+nSR)/(nAF*nSR) + g^2/(2*(nAF+nSR))
#   - Pooled: rma(method="REML", test="knha")
#
# Outputs:
#   results/figure6_meta_axis.tsv         per-(axis, cohort) g/SE + pooled
#   results/figure6_meta_axis_loco.tsv    leave-one-cohort-out (supp)
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(metafor)
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
IN <- file.path(ROOT, "results", "figure6_axis_composite_long.tsv")
OUT_MAIN <- file.path(ROOT, "results", "figure6_meta_axis.tsv")
OUT_LOCO <- file.path(ROOT, "results", "figure6_meta_axis_loco.tsv")

cat(sprintf("ROOT = %s\n", ROOT))
comp <- fread(IN, showProgress = FALSE)
comp <- comp[!is.na(composite)]

hedges_g <- function(x_AF, x_SR) {
  nA <- length(x_AF); nS <- length(x_SR)
  if (nA < 2 || nS < 2) return(list(g = NA_real_, var = NA_real_,
                                    n_AF = nA, n_SR = nS,
                                    mean_AF = mean(x_AF, na.rm = TRUE),
                                    sd_AF = sd(x_AF, na.rm = TRUE),
                                    mean_SR = mean(x_SR, na.rm = TRUE),
                                    sd_SR = sd(x_SR, na.rm = TRUE)))
  m_A <- mean(x_AF); s_A <- sd(x_AF)
  m_S <- mean(x_SR); s_S <- sd(x_SR)
  df <- nA + nS - 2
  s_p <- sqrt(((nA - 1) * s_A^2 + (nS - 1) * s_S^2) / df)
  if (!is.finite(s_p) || s_p == 0) return(list(g = NA_real_, var = NA_real_,
                                                n_AF = nA, n_SR = nS,
                                                mean_AF = m_A, sd_AF = s_A,
                                                mean_SR = m_S, sd_SR = s_S))
  d <- (m_A - m_S) / s_p
  J <- 1 - 3 / (4 * df - 1)
  g <- J * d
  v <- (nA + nS) / (nA * nS) + g^2 / (2 * (nA + nS))
  list(g = g, var = v, n_AF = nA, n_SR = nS,
       mean_AF = m_A, sd_AF = s_A, mean_SR = m_S, sd_SR = s_S)
}

axes <- sort(unique(comp$axis))
all_rows <- list()
loco_rows <- list()

for (ax in axes) {
  per_cohort <- comp[axis == ax, .(
    list_AF = list(composite[rhythm == "AF"]),
    list_SR = list(composite[rhythm == "SR"])
  ), by = cohort]

  per_cohort_g <- per_cohort[, {
    res <- hedges_g(unlist(list_AF), unlist(list_SR))
    list(
      n_AF = res$n_AF, n_SR = res$n_SR,
      mean_AF = res$mean_AF, sd_AF = res$sd_AF,
      mean_SR = res$mean_SR, sd_SR = res$sd_SR,
      g = res$g, var_g = res$var, se_g = sqrt(res$var)
    )
  }, by = cohort]
  per_cohort_g[, axis := ax]
  per_cohort_g[, ci_lo := g - 1.96 * se_g]
  per_cohort_g[, ci_hi := g + 1.96 * se_g]
  per_cohort_g[, type := "cohort"]
  per_cohort_g[, n_studies := NA_integer_]
  per_cohort_g[, p_HK := NA_real_]
  per_cohort_g[, tau2 := NA_real_]
  per_cohort_g[, I2   := NA_real_]
  per_cohort_g[, Q := NA_real_]
  per_cohort_g[, df_Q := NA_integer_]

  rma_input <- per_cohort_g[!is.na(g)]
  if (nrow(rma_input) >= 2) {
    fit <- rma(yi = rma_input$g, vi = rma_input$var_g,
               method = "REML", test = "knha")
    pool <- data.table(
      cohort = "RE pooled (HK)",
      n_AF = sum(rma_input$n_AF), n_SR = sum(rma_input$n_SR),
      mean_AF = NA_real_, sd_AF = NA_real_,
      mean_SR = NA_real_, sd_SR = NA_real_,
      g = as.numeric(fit$beta), var_g = as.numeric(fit$se)^2,
      se_g = as.numeric(fit$se),
      axis = ax,
      ci_lo = as.numeric(fit$ci.lb), ci_hi = as.numeric(fit$ci.ub),
      type = "pooled", n_studies = fit$k,
      p_HK = fit$pval, tau2 = fit$tau2, I2 = fit$I2,
      Q = fit$QE, df_Q = fit$k - 1
    )
    setcolorder(pool, names(per_cohort_g))
    all_rows[[ax]] <- rbind(per_cohort_g, pool)

    loco <- leave1out(fit)
    loco_rows[[ax]] <- data.table(
      axis = ax,
      excluded_cohort = rma_input$cohort,
      g = loco$estimate, se = loco$se, ci_lo = loco$ci.lb, ci_hi = loco$ci.ub,
      tau2 = loco$tau2, I2 = loco$I2, p = loco$pval
    )
  } else {
    all_rows[[ax]] <- per_cohort_g
  }
}

out <- rbindlist(all_rows, use.names = TRUE)
fwrite(out, OUT_MAIN, sep = "\t")

if (length(loco_rows)) {
  fwrite(rbindlist(loco_rows, use.names = TRUE), OUT_LOCO, sep = "\t")
} else {
  fwrite(data.table(), OUT_LOCO, sep = "\t")
}

cat("\n--- 6A meta result ---\n")
print(out[, .(axis, cohort, n_AF, n_SR, g = round(g, 3),
              ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3),
              I2 = round(I2, 1), tau2 = round(tau2, 3),
              p_HK = ifelse(is.na(p_HK), NA_character_, format.pval(p_HK, 3)))])

cat("\nH2 DONE.\n")
