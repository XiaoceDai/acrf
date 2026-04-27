#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# H9 — Leave-one-program-out (LOPO) robustness of the dual-axis meta.
#
# For STRUCTURAL and METABOLIC axes (dual-axis support programs):
#   For each axis, in turn remove one of its 3 frozen support programs and
#   recompute donor-level axis composite + cohort-level Hedges' g + REML/
#   Hartung-Knapp pool, exactly as in Fig 6A (h2_meta_axis.R).
#
# Output:
#   results/figure_s8d_lopo_axis.tsv
#       (axis, removed_program, remaining_programs, k_remain,
#        pooled_g, ci_lo, ci_hi, I2, tau2, p_HK,
#        direction_preserved, robustness_label)
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
COHORTS <- c("GSE128188", "GSE79768", "GSE115574")

# ---- shared utilities ------------------------------------------------------
hedges_g <- function(x_AF, x_SR) {
  nA <- length(x_AF); nS <- length(x_SR)
  if (nA < 2 || nS < 2) return(list(g = NA_real_, var = NA_real_))
  m_A <- mean(x_AF); s_A <- sd(x_AF)
  m_S <- mean(x_SR); s_S <- sd(x_SR)
  df <- nA + nS - 2
  s_p <- sqrt(((nA - 1) * s_A^2 + (nS - 1) * s_S^2) / df)
  if (!is.finite(s_p) || s_p == 0) return(list(g = NA_real_, var = NA_real_))
  d <- (m_A - m_S) / s_p
  J <- 1 - 3 / (4 * df - 1)
  g <- J * d
  v <- (nA + nS) / (nA * nS) + g^2 / (2 * (nA + nS))
  list(g = g, var = v)
}

# ---- inputs ----------------------------------------------------------------
long <- fread(file.path(ROOT, "Project_P1", "program_scores_long.tsv"))
sets <- fread(file.path(RES, "figure5_program_sets.tsv"))

axes <- list(
  STRUCTURAL = list(ids = sets[set_name == "STRUCTURAL", program_id],
                    names = sets[set_name == "STRUCTURAL", program_label]),
  METABOLIC  = list(ids = sets[set_name == "METABOLIC",  program_id],
                    names = sets[set_name == "METABOLIC",  program_label])
)

# Map ID → label for pretty output
id_to_label <- setNames(sets$program_label, sets$program_id)

# ---- core: composite + pooled g for an arbitrary set of program IDs --------
axis_pool <- function(ids) {
  sub <- long[program_id %in% ids,
              .(composite = mean(aCAI_rel, na.rm = TRUE)),
              by = .(cohort, donor_id, rhythm)]
  rec <- list()
  for (.coh in COHORTS) {
    ss <- sub[cohort == .coh]
    g <- hedges_g(ss$composite[ss$rhythm == "AF"],
                  ss$composite[ss$rhythm == "SR"])
    rec[[.coh]] <- data.table(cohort = .coh, g = g$g, var = g$var)
  }
  cc <- rbindlist(rec, use.names = TRUE)
  ok <- !is.na(cc$g) & !is.na(cc$var) & cc$var > 0
  if (sum(ok) < 2) return(list(g = NA_real_, ci_lo = NA_real_,
                               ci_hi = NA_real_, I2 = NA_real_,
                               tau2 = NA_real_, p_HK = NA_real_,
                               n_AF = sum(cc$g != 0, na.rm = TRUE)))
  fit <- tryCatch(rma(yi = cc$g[ok], vi = cc$var[ok],
                      method = "REML", test = "knha"),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(g = NA_real_, ci_lo = NA_real_,
                                ci_hi = NA_real_, I2 = NA_real_,
                                tau2 = NA_real_, p_HK = NA_real_,
                                n_AF = NA_integer_))
  list(g = as.numeric(fit$beta), ci_lo = as.numeric(fit$ci.lb),
       ci_hi = as.numeric(fit$ci.ub),
       I2 = as.numeric(fit$I2), tau2 = as.numeric(fit$tau2),
       p_HK = fit$pval, n_AF = NA_integer_)
}

# ---- run LOPO --------------------------------------------------------------
rows <- list()
for (axis_lab in names(axes)) {
  ids_full <- axes[[axis_lab]]$ids
  cat(sprintf("\n=== LOPO axis %s  (n_programs = %d) ===\n",
              axis_lab, length(ids_full)))
  # full axis (reference)
  full <- axis_pool(ids_full)
  rec <- data.table(
    axis = axis_lab, removed_program = "(none)",
    removed_label = "FULL axis",
    remaining_programs = paste(id_to_label[ids_full], collapse = " + "),
    k_remain = length(ids_full),
    pooled_g = full$g, ci_lo = full$ci_lo, ci_hi = full$ci_hi,
    I2 = full$I2, tau2 = full$tau2, p_HK = full$p_HK
  )
  rows[[paste0(axis_lab, "_FULL")]] <- rec
  cat(sprintf("  FULL  g=%+.3f  CI=[%+.2f, %+.2f]  I²=%.0f%%  p_HK=%.3f\n",
              full$g, full$ci_lo, full$ci_hi, full$I2, full$p_HK))

  for (k in seq_along(ids_full)) {
    rm_id <- ids_full[k]
    keep <- ids_full[-k]
    res <- axis_pool(keep)
    rec <- data.table(
      axis = axis_lab, removed_program = rm_id,
      removed_label = id_to_label[rm_id],
      remaining_programs = paste(id_to_label[keep], collapse = " + "),
      k_remain = length(keep),
      pooled_g = res$g, ci_lo = res$ci_lo, ci_hi = res$ci_hi,
      I2 = res$I2, tau2 = res$tau2, p_HK = res$p_HK
    )
    rows[[paste(axis_lab, rm_id, sep = "_")]] <- rec
    cat(sprintf("  -%-25s g=%+.3f  CI=[%+.2f, %+.2f]  I²=%.0f%%  p_HK=%.3f\n",
                substr(id_to_label[rm_id], 1, 24),
                res$g, res$ci_lo, res$ci_hi,
                ifelse(is.na(res$I2), 0, res$I2),
                ifelse(is.na(res$p_HK), NA, res$p_HK)))
  }
}
out <- rbindlist(rows, use.names = TRUE)

# ---- direction-preservation tag --------------------------------------------
full_signs <- out[removed_program == "(none)",
                  .(full_g = pooled_g, full_sign = sign(pooled_g)),
                  by = axis]
out <- merge(out, full_signs, by = "axis", all.x = TRUE)
out[, direction_preserved :=
      ifelse(is.na(pooled_g) | is.na(full_sign), NA,
             sign(pooled_g) == full_sign)]

# robustness label per axis
robustness <- out[removed_program != "(none)",
                  .(n_total = .N,
                    n_pres  = sum(direction_preserved, na.rm = TRUE)),
                  by = axis]
robustness[, robustness_label := ifelse(
  n_pres == n_total, "robust",
  ifelse(n_pres == n_total - 1, "partially-sensitive",
         "single-program-dependent"))]
out <- merge(out, robustness[, .(axis, robustness_label)], by = "axis",
             all.x = TRUE)

fwrite(out, file.path(RES, "figure_s8d_lopo_axis.tsv"), sep = "\t")
cat("\n--- robustness summary ---\n")
print(unique(out[, .(axis, robustness_label)]))
cat("\n--- full table ---\n")
print(out[, .(axis, removed_label, k_remain,
              g = round(pooled_g, 3),
              ci = sprintf("[%+.2f,%+.2f]", ci_lo, ci_hi),
              I2 = round(I2, 0),
              dir_preserved = direction_preserved,
              robust = robustness_label)])

cat("\nH9 DONE.\n")
