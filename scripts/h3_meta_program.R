#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# H3 — Per-program random-effects meta on aCAI_rel (AF − SR), 10 frozen
#       chamber programs.  Inputs and method follow H2 §2/§3:
#         g  = J(df) * (mean_AF - mean_SR) / s_pooled   (Hedges, AF - SR)
#         RE pooled: rma(method="REML", test="knha")
# Spec: FIGURE6_WORKFLOW.md §5 (program-level only).
#
# Outputs:
#   results/figure6_meta_program.tsv
#   results/figure6_meta_program_loco.tsv
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
LONG <- file.path(ROOT, "Project_P1", "program_scores_long.tsv")
PANEL <- file.path(ROOT, "results", "figure3_program_panel.tsv")
SETS <- file.path(ROOT, "results", "figure5_program_sets.tsv")
OUT_MAIN <- file.path(ROOT, "results", "figure6_meta_program.tsv")
OUT_LOCO <- file.path(ROOT, "results", "figure6_meta_program_loco.tsv")

cat(sprintf("ROOT = %s\n", ROOT))

long <- fread(LONG, showProgress = FALSE)
panel <- fread(PANEL, showProgress = FALSE)
sets <- fread(SETS, showProgress = FALSE)
struct_ids <- sets[set_name == "STRUCTURAL", program_id]
metab_ids  <- sets[set_name == "METABOLIC",  program_id]

# Axis membership lookup
axis_of <- function(pid) {
  if (pid %in% struct_ids) return("STRUCTURAL")
  if (pid %in% metab_ids)  return("METABOLIC")
  "panel-only"
}

panel <- panel[order(program_order)]
panel[, axis := vapply(program_id, axis_of, character(1))]
cat("\nPanel × axis assignment:\n")
print(panel[, .(program_order, program_label, axis)])

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

all_rows <- list()
loco_rows <- list()

for (i in seq_len(nrow(panel))) {
  pid <- panel$program_id[i]
  plab <- panel$program_label[i]
  pord <- panel$program_order[i]
  axis_lab <- panel$axis[i]

  d <- long[program_id == pid & !is.na(aCAI_rel)]
  if (nrow(d) == 0) next

  per_cohort <- d[, .(
    list_AF = list(aCAI_rel[rhythm == "AF"]),
    list_SR = list(aCAI_rel[rhythm == "SR"])
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
  per_cohort_g[, program_id := pid]
  per_cohort_g[, program_label := plab]
  per_cohort_g[, program_order := pord]
  per_cohort_g[, axis := axis_lab]
  per_cohort_g[, ci_lo := g - 1.96 * se_g]
  per_cohort_g[, ci_hi := g + 1.96 * se_g]
  per_cohort_g[, type := "cohort"]
  per_cohort_g[, n_studies := NA_integer_]
  per_cohort_g[, p_HK := NA_real_]
  per_cohort_g[, tau2 := NA_real_]
  per_cohort_g[, I2 := NA_real_]
  per_cohort_g[, Q := NA_real_]
  per_cohort_g[, df_Q := NA_integer_]

  rma_input <- per_cohort_g[!is.na(g)]
  if (nrow(rma_input) >= 2) {
    fit <- tryCatch(
      rma(yi = rma_input$g, vi = rma_input$var_g,
          method = "REML", test = "knha"),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      pool <- data.table(
        cohort = "RE pooled (HK)",
        n_AF = sum(rma_input$n_AF), n_SR = sum(rma_input$n_SR),
        mean_AF = NA_real_, sd_AF = NA_real_,
        mean_SR = NA_real_, sd_SR = NA_real_,
        g = as.numeric(fit$beta), var_g = as.numeric(fit$se)^2,
        se_g = as.numeric(fit$se),
        program_id = pid, program_label = plab,
        program_order = pord, axis = axis_lab,
        ci_lo = as.numeric(fit$ci.lb), ci_hi = as.numeric(fit$ci.ub),
        type = "pooled", n_studies = fit$k,
        p_HK = fit$pval, tau2 = fit$tau2, I2 = fit$I2,
        Q = fit$QE, df_Q = fit$k - 1
      )
      setcolorder(pool, names(per_cohort_g))
      all_rows[[pid]] <- rbind(per_cohort_g, pool)

      loco <- leave1out(fit)
      loco_rows[[pid]] <- data.table(
        program_id = pid, program_label = plab, program_order = pord,
        excluded_cohort = rma_input$cohort,
        g = loco$estimate, se = loco$se,
        ci_lo = loco$ci.lb, ci_hi = loco$ci.ub,
        tau2 = loco$tau2, I2 = loco$I2, p = loco$pval
      )
    } else {
      all_rows[[pid]] <- per_cohort_g
    }
  } else {
    all_rows[[pid]] <- per_cohort_g
  }
}

out <- rbindlist(all_rows, use.names = TRUE)
fwrite(out, OUT_MAIN, sep = "\t")

if (length(loco_rows)) {
  fwrite(rbindlist(loco_rows, use.names = TRUE), OUT_LOCO, sep = "\t")
}

cat("\n--- Pooled (HK) summary, ordered by program_order ---\n")
pool <- out[type == "pooled"]
pool <- pool[order(program_order)]
print(pool[, .(program_order, program_label, axis,
               n_AF, n_SR,
               g = round(g, 3),
               ci_lo = round(ci_lo, 3),
               ci_hi = round(ci_hi, 3),
               I2 = round(I2, 1),
               tau2 = round(tau2, 3),
               p_HK = ifelse(is.na(p_HK), NA_character_, format.pval(p_HK, 3)))])

cat("\nH3 DONE.\n")
