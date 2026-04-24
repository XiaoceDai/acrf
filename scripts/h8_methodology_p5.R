#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# H8 — Method-validity P5 supplement.
#
# Three sub-modules (selectable via CLI):
#   sim     — simulation power curves for the axis-meta pipeline
#   boot    — cohort-donor bootstrap CI for STRUCTURAL / METABOLIC pooled g
#   classic — classical per-cohort limma chamber × rhythm comparison vs
#             our aCAI_rel-based pooled g (direction concordance)
#
# Usage:
#   Rscript scripts/h8_methodology_p5.R [sim] [boot] [classic]
#   default = boot + classic (sim is slowest)
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(metafor); library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
do_run <- if (length(args)) tolower(args) else c("boot", "classic")
if ("all" %in% do_run) do_run <- c("sim", "boot", "classic")

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
RES <- file.path(ROOT, "results"); dir.create(RES, showWarnings = FALSE)
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

reml_pool_g <- function(g_vec, v_vec, knha = TRUE) {
  ok <- !is.na(g_vec) & !is.na(v_vec) & v_vec > 0
  if (sum(ok) < 2) return(NA_real_)
  fit <- tryCatch(rma(yi = g_vec[ok], vi = v_vec[ok],
                      method = "REML", test = if (knha) "knha" else "z"),
                  error = function(e) NULL)
  if (is.null(fit)) NA_real_ else as.numeric(fit$beta)
}

# Composite aCAI_rel per donor per axis from the existing long table
load_long <- function() {
  long <- fread(file.path(P1, "program_scores_long.tsv"))
  sets <- fread(file.path(RES, "figure5_program_sets.tsv"))
  list(long = long, sets = sets)
}

axis_composite_long <- function(long, sets) {
  S <- sets[set_name == "STRUCTURAL", program_id]
  M <- sets[set_name == "METABOLIC",  program_id]
  rec <- list()
  for (axis_lab in c("STRUCTURAL", "METABOLIC")) {
    ids <- if (axis_lab == "STRUCTURAL") S else M
    sub <- long[program_id %in% ids,
                .(composite = mean(aCAI_rel, na.rm = TRUE)),
                by = .(cohort, donor_id, rhythm)]
    sub[, axis := axis_lab]
    rec[[axis_lab]] <- sub
  }
  rbindlist(rec, use.names = TRUE)
}

axis_pooled_g <- function(comp) {
  rec <- list()
  for (axis_lab in c("STRUCTURAL", "METABOLIC")) {
    sub <- comp[axis == axis_lab]
    gs <- numeric(length(COHORTS)); vs <- numeric(length(COHORTS))
    for (k in seq_along(COHORTS)) {
      .coh <- COHORTS[k]
      ss <- sub[cohort == .coh]
      g <- hedges_g(ss$composite[ss$rhythm == "AF"],
                    ss$composite[ss$rhythm == "SR"])
      gs[k] <- g$g; vs[k] <- g$var
    }
    rec[[axis_lab]] <- list(axis = axis_lab,
                            g = reml_pool_g(gs, vs))
  }
  rec
}

# ===========================================================================
# Module: BOOTSTRAP — cohort-donor bootstrap for STRUCTURAL / METABOLIC g
# ===========================================================================
mod_bootstrap <- function(B = 2000L, seed = 20260420L) {
  cat("\n=== P5c: cohort-donor bootstrap of axis pooled g ===\n")
  set.seed(seed)
  ll <- load_long()
  comp <- axis_composite_long(ll$long, ll$sets)

  rec <- matrix(NA_real_, nrow = B, ncol = 2,
                dimnames = list(NULL, c("STRUCTURAL", "METABOLIC")))
  donors_per_cohort <- comp[, unique(donor_id), by = cohort]

  t0 <- Sys.time()
  for (b in seq_len(B)) {
    boot_long <- list()
    for (.coh in COHORTS) {
      donors <- donors_per_cohort[cohort == .coh, V1]
      if (length(donors) < 3) next
      bs_donors <- sample(donors, length(donors), replace = TRUE)
      # gather their rows; preserve donor multiplicity by tagging
      sub_cohort <- comp[cohort == .coh & axis %in% c("STRUCTURAL", "METABOLIC")]
      take <- list()
      for (i in seq_along(bs_donors)) {
        rows <- sub_cohort[donor_id == bs_donors[i]]
        if (nrow(rows) == 0) next
        rows <- copy(rows); rows[, donor_id := paste0(bs_donors[i], "_b", i)]
        take[[i]] <- rows
      }
      boot_long[[.coh]] <- rbindlist(take, use.names = TRUE)
    }
    bp <- rbindlist(boot_long, use.names = TRUE)
    res <- axis_pooled_g(bp)
    rec[b, 1] <- res[["STRUCTURAL"]]$g
    rec[b, 2] <- res[["METABOLIC"]]$g
    if (b %% 500 == 0) cat(sprintf("  b=%d (%.1fs)\n", b,
                                   as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }

  # real values
  real <- axis_pooled_g(comp)
  out <- data.table(
    axis = c("STRUCTURAL", "METABOLIC"),
    real_g = c(real[["STRUCTURAL"]]$g, real[["METABOLIC"]]$g),
    boot_mean = colMeans(rec, na.rm = TRUE),
    boot_sd = apply(rec, 2, sd, na.rm = TRUE),
    boot_q025 = apply(rec, 2, quantile, 0.025, na.rm = TRUE),
    boot_q975 = apply(rec, 2, quantile, 0.975, na.rm = TRUE),
    n_boot = B
  )
  fwrite(out, file.path(RES, "method_p5_bootstrap_axis.tsv"), sep = "\t")
  fwrite(as.data.table(rec)[, b := .I],
         file.path(RES, "method_p5_bootstrap_axis_distribution.tsv.gz"),
         sep = "\t", compress = "gzip")
  cat("\n--- bootstrap summary ---\n"); print(out)
}

# ===========================================================================
# Module: CLASSIC — per-cohort limma per-chamber AF vs SR + concordance
#                   with our framework's per-program direction.
# ===========================================================================
mod_classic <- function() {
  cat("\n=== P5b: classical limma comparison ===\n")
  ll <- load_long()
  long <- ll$long
  sets <- ll$sets

  # Reference: per-program per-cohort donor-level Welch t on aCAI_rel
  # (matches our framework, but using the simpler t-statistic instead
  # of Hedges' g — useful as a sanity replicate).
  ref <- long[, {
    t_test <- tryCatch(t.test(aCAI_rel[rhythm == "AF"],
                              aCAI_rel[rhythm == "SR"], var.equal = FALSE),
                       error = function(e) NULL)
    if (is.null(t_test)) list(t_stat = NA_real_, t_p = NA_real_,
                              dir_framework = NA_integer_) else
      list(t_stat = unname(t_test$statistic),
           t_p   = t_test$p.value,
           dir_framework = sign(unname(t_test$statistic)))
  }, by = .(cohort, program_id, program_label)]
  cat("framework Welch direction (sign of mean_AF − mean_SR aCAI_rel):\n")
  print(ref[, .(cohort, program_label,
                t_stat = round(t_stat, 2),
                t_p = signif(t_p, 2),
                dir = dir_framework)])

  # Classical: per cohort, per chamber (LA + RA), classical limma AF vs SR
  # at the gene level, then per-program Welch on logFC means.  This is the
  # "naive" approach; we test concordance of program-level direction with
  # the framework's t direction.
  classical_rows <- list()
  for (cohort in COHORTS) {
    expr <- fread(file.path(P1, cohort, "expr_symbol.tsv"),
                  sep = "\t", header = TRUE)
    sym <- expr[[1]]; M <- as.matrix(expr[, -1, with = FALSE])
    rownames(M) <- sym
    if (anyDuplicated(sym)) {
      v <- apply(M, 1, var, na.rm = TRUE)
      keep <- tapply(seq_along(sym), sym, function(idx) idx[which.max(v[idx])])
      M <- M[unlist(keep, use.names = FALSE), , drop = FALSE]
    }
    meta <- fread(file.path(P1, cohort, "sample_metadata.tsv"))
    meta <- meta[match(colnames(M), gsm)]

    # per-chamber AF vs SR limma
    for (chamber_lab in c("LA", "RA")) {
      sel <- meta$chamber == chamber_lab & !is.na(meta$rhythm)
      if (sum(sel) < 4) next
      Ms <- M[, sel]; ms <- meta[sel]
      design <- model.matrix(~ 0 + ms$rhythm)
      colnames(design) <- gsub("ms\\$rhythm", "", colnames(design))
      if (!all(c("AF", "SR") %in% colnames(design))) next
      fit <- lmFit(Ms, design)
      cont <- makeContrasts(AFvsSR = AF - SR, levels = design)
      fit2 <- eBayes(contrasts.fit(fit, cont))
      tt <- topTable(fit2, number = Inf, sort.by = "none")
      tt$gene <- rownames(tt)
      classical_rows[[paste(cohort, chamber_lab, sep = "_")]] <-
        data.table(cohort = cohort, chamber = chamber_lab,
                   gene = tt$gene, logFC = tt$logFC,
                   t_stat = tt$t, P = tt$P.Value)
    }
  }
  classical <- rbindlist(classical_rows, use.names = TRUE)
  fwrite(classical, file.path(RES, "method_p5_classical_de.tsv.gz"),
         sep = "\t", compress = "gzip")

  # Aggregate per-program: for each (cohort, program), program-level
  # logFC = mean of LE-gene logFC in the LA-arm minus mean in the RA-arm.
  prog <- read.delim(file.path(ROOT, "normal_chamber_programs.tsv"),
                     sep = "\t", stringsAsFactors = FALSE)
  for (col in c("program_id", "leading_edge_genes"))
    if (col %in% colnames(prog)) prog[[col]] <- gsub('^"|"$', '', prog[[col]])
  panel <- fread(file.path(RES, "figure3_program_panel.tsv"))
  prog <- merge(prog[, c("program_id", "direction", "leading_edge_genes")],
                panel[, .(program_id, program_label, baseline_direction)],
                by = "program_id", sort = FALSE)
  prog$le_genes <- lapply(prog$leading_edge_genes,
                          function(s) unique(trimws(strsplit(s, ",", fixed = TRUE)[[1]])))

  prog_rows <- list()
  for (i in seq_len(nrow(prog))) {
    le <- prog$le_genes[[i]]; bdir <- prog$baseline_direction[i]
    pid <- prog$program_id[i]; plab <- prog$program_label[i]
    for (.coh in COHORTS) {
      la_lfc <- classical[cohort == .coh & chamber == "LA" & gene %in% le, mean(logFC, na.rm = TRUE)]
      ra_lfc <- classical[cohort == .coh & chamber == "RA" & gene %in% le, mean(logFC, na.rm = TRUE)]
      if (is.nan(la_lfc) || is.nan(ra_lfc)) next
      classical_dir <- bdir * (la_lfc - ra_lfc)  # aligned to NCP direction
      prog_rows[[paste(.coh, pid, sep = "_")]] <- data.table(
        cohort = .coh, program_id = pid, program_label = plab,
        baseline_dir = bdir,
        mean_logFC_LA = la_lfc, mean_logFC_RA = ra_lfc,
        classical_aligned = classical_dir,
        sign_classical = sign(classical_dir)
      )
    }
  }
  classical_prog <- rbindlist(prog_rows, use.names = TRUE)
  classical_prog <- merge(classical_prog,
                          ref[, .(cohort, program_id,
                                  framework_aligned = t_stat,
                                  sign_framework = dir_framework)],
                          by = c("cohort", "program_id"), all.x = TRUE)
  classical_prog[, concordant := sign_classical == sign_framework]
  fwrite(classical_prog,
         file.path(RES, "method_p5_classical_vs_framework.tsv"), sep = "\t")

  # summary concordance per cohort
  conc <- classical_prog[, .(n_concord = sum(concordant, na.rm = TRUE),
                             n_total   = sum(!is.na(concordant)),
                             pct       = round(100 * mean(concordant, na.rm = TRUE), 1)),
                         by = cohort]
  cat("\n--- classical-vs-framework direction concordance ---\n")
  print(conc)
  cat("\n--- per-program concordance summary ---\n")
  print(classical_prog[, .(cohort, program_label,
                           classical = round(classical_aligned, 3),
                           framework = round(framework_aligned, 2),
                           concord = concordant)])
}

# ===========================================================================
# Module: SIMULATION — power curve for axis-pooled g detection
# ===========================================================================
mod_simulation <- function(reps = 200L, n_genes = 200L,
                            n_LE = 30L, sigma_gene = 1.0,
                            seed = 20260421L) {
  cat("\n=== P5a: simulation power curves ===\n")
  set.seed(seed)
  effects <- c(0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.70, 1.00)

  # Cohort donor structure (mirrors paired-chamber cohorts)
  cohort_n <- list(GSE128188 = c(SR = 5, AF = 5),
                   GSE79768  = c(SR = 6, AF = 7),
                   GSE115574 = c(SR = 15, AF = 13))

  sim_one <- function(eff) {
    # For each cohort, simulate donors with paired LA/RA samples.
    # Healthy program direction: LE genes are LA-up by sigma_gene * 0.5.
    # AF rewiring: aCAI of the program shifts by eff (negative = attenuation).
    rec <- list()
    for (cohort in names(cohort_n)) {
      ns <- cohort_n[[cohort]]
      donors <- c(rep("SR", ns["SR"]), rep("AF", ns["AF"]))
      # Simulate per-donor LA - RA difference of mean LE gene expression
      cai <- sapply(donors, function(rh) {
        base <- 0.5  # healthy LA-RA baseline difference per donor
        af_shift <- if (rh == "AF") -eff else 0  # attenuation when AF
        rnorm(1, mean = base + af_shift, sd = 0.40)
      })
      bulk_mag <- 0.5
      acai_rel <- cai / bulk_mag
      g <- hedges_g(acai_rel[donors == "AF"], acai_rel[donors == "SR"])
      rec[[cohort]] <- data.table(cohort = cohort, g = g$g, var = g$var)
    }
    cc <- rbindlist(rec, use.names = TRUE)
    reml_pool_g(cc$g, cc$var)
  }

  rows <- list()
  for (eff in effects) {
    pooled <- replicate(reps, sim_one(eff))
    detect_rate <- mean(abs(pooled) > 0.30, na.rm = TRUE)
    rows[[as.character(eff)]] <- data.table(
      effect_size = eff,
      mean_pooled_g = mean(pooled, na.rm = TRUE),
      sd_pooled_g = sd(pooled, na.rm = TRUE),
      detect_rate = detect_rate,
      reps = reps
    )
    cat(sprintf("  eff=%.2f  mean g = %+.3f  detect rate (|g|>0.3) = %.2f\n",
                eff, mean(pooled, na.rm = TRUE), detect_rate))
  }
  out <- rbindlist(rows, use.names = TRUE)
  fwrite(out, file.path(RES, "method_p5_simulation.tsv"), sep = "\t")
}

# ----- run requested modules ------------------------------------------------
if ("boot"    %in% do_run) mod_bootstrap()
if ("classic" %in% do_run) mod_classic()
if ("sim"     %in% do_run) mod_simulation()
cat("\nH8 DONE.\n")
