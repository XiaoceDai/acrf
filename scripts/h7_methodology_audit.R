#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# H7 — Methodology audit for the chamber-aware rewiring framework.
#
# Four internal-consistency tests (all on already-loaded data):
#   A1  NCP self-recovery:   on GSE112339 healthy bulk,
#                            aCAI_rel(program) should be ≈ +1.
#   A2  Direction stability: fraction of healthy donors with aCAI > 0
#                            (signed in NCP direction).
#   A3  Cross-method:        Spearman ρ of singscore vs mean LE-expression
#                            score on the healthy bulk (per program).
#   A4  Cross-method robustness: replace singscore with mean-z scoring in
#                                Fig 2 + Fig 6A pipeline; compare pooled g
#                                per axis, both methods.
#
# Output:
#   results/methodology_audit_program.tsv   (A1+A2+A3, per program)
#   results/methodology_audit_axis_g.tsv    (A4, per axis × method)
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(singscore); library(metafor)
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
RES <- file.path(ROOT, "results"); dir.create(RES, showWarnings = FALSE)
COHORTS <- c("GSE128188", "GSE79768", "GSE115574")

# ---- programs + LE genes ---------------------------------------------------
prog <- read.delim(file.path(ROOT, "normal_chamber_programs.tsv"),
                   sep = "\t", stringsAsFactors = FALSE)
for (col in c("program_id", "leading_edge_genes"))
  if (col %in% colnames(prog)) prog[[col]] <- gsub('^"|"$', '', prog[[col]])
panel <- fread(file.path(RES, "figure3_program_panel.tsv"))
prog <- merge(prog[, c("program_id", "direction", "leading_edge_genes")],
              panel[, .(program_id, program_label, program_order,
                        baseline_direction)],
              by = "program_id", sort = FALSE)
setorder(prog, program_order)
prog$le_genes <- lapply(prog$leading_edge_genes,
                        function(s) unique(trimws(strsplit(s, ",", fixed = TRUE)[[1]])))
prog$n_LE <- vapply(prog$le_genes, length, integer(1))

# ---- load matrices ---------------------------------------------------------
load_mat <- function(path) {
  expr <- fread(path, sep = "\t", header = TRUE)
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
  Mc <- load_mat(file.path(P1, cohort, "expr_symbol.tsv"))
  mc <- fread(file.path(P1, cohort, "sample_metadata.tsv"))
  mats[[cohort]] <- list(M = Mc, rk = rankGenes(Mc), meta = mc)
  cat(sprintf("  %s: %d genes x %d samples\n", cohort, nrow(Mc), ncol(Mc)))
}
# healthy bulk
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
mats[["GSE112339"]] <- list(M = Mb, rk = rankGenes(Mb), meta = bulk_meta)

# ---- scoring functions -----------------------------------------------------
singscore_one <- function(rk, gene_set) {
  hit <- intersect(gene_set, rownames(rk))
  if (length(hit) < 3) return(rep(NA_real_, ncol(rk)))
  s <- simpleScore(rk, upSet = hit, knownDirection = TRUE,
                   centerScore = TRUE)
  setNames(s$TotalScore, rownames(s))
}

mean_z_one <- function(M, gene_set) {
  hit <- intersect(gene_set, rownames(M))
  if (length(hit) < 3) return(rep(NA_real_, ncol(M)))
  Z <- t(scale(t(M[hit, , drop = FALSE]), center = TRUE, scale = TRUE))
  Z[!is.finite(Z)] <- NA_real_
  colMeans(Z, na.rm = TRUE)
}

# ---- per-donor aCAI helper -------------------------------------------------
acai_donors <- function(scores, meta, baseline_dir,
                        sample_col = "gsm",
                        donor_col = "donor_id",
                        chamber_col = "chamber",
                        rhythm_col = NA) {
  meta <- as.data.table(meta)
  meta[, score := scores[match(meta[[sample_col]], names(scores))]]
  agg_cols <- c(donor_col)
  if (!is.na(rhythm_col) && rhythm_col %in% names(meta))
    agg_cols <- c(agg_cols, rhythm_col)
  by_donor <- meta[!is.na(score),
                   .(LA = mean(score[get(chamber_col) %in% c("LA", "left atrium", "left_atrium")],
                               na.rm = TRUE),
                     RA = mean(score[get(chamber_col) %in% c("RA", "right atrium", "right_atrium")],
                               na.rm = TRUE),
                     rhythm = if (!is.na(rhythm_col)) get(rhythm_col)[1] else NA_character_),
                   by = agg_cols]
  setnames(by_donor, donor_col, "donor_id")
  by_donor <- by_donor[!is.na(LA) & !is.na(RA)]
  by_donor[, CAI := LA - RA]
  by_donor[, aCAI := baseline_dir * CAI]
  by_donor
}

# ---- A1-A3 per program (healthy bulk) --------------------------------------
audit_program_rows <- list()
cat("\n=== A1 + A2 + A3 (per program, healthy bulk) ===\n")
for (i in seq_len(nrow(prog))) {
  pid <- prog$program_id[i]
  bd  <- prog$baseline_direction[i]
  le  <- prog$le_genes[[i]]
  s_ss   <- singscore_one(mats[["GSE112339"]]$rk, le)
  s_mz   <- mean_z_one(mats[["GSE112339"]]$M, le)
  donors_ss <- acai_donors(s_ss, mats[["GSE112339"]]$meta, bd,
                           sample_col = "sample_id",
                           donor_col = "donor",
                           chamber_col = "chamber")
  donors_mz <- acai_donors(s_mz, mats[["GSE112339"]]$meta, bd,
                           sample_col = "sample_id",
                           donor_col = "donor",
                           chamber_col = "chamber")
  bulk_mag_ss <- abs(mean(donors_ss$aCAI, na.rm = TRUE))
  bulk_mag_mz <- abs(mean(donors_mz$aCAI, na.rm = TRUE))
  # NCP self-recovery: aCAI_rel on healthy donors should be ≈ +1
  acai_rel_ss <- donors_ss$aCAI / bulk_mag_ss
  acai_rel_mz <- donors_mz$aCAI / bulk_mag_mz
  # Direction stability: fraction of donors with aCAI > 0
  frac_pos_ss <- mean(donors_ss$aCAI > 0, na.rm = TRUE)
  # Cross-method correlation on healthy donors (aCAI numbers, signed)
  ok <- !is.na(donors_ss$aCAI) & !is.na(donors_mz$aCAI)
  rho_methods <- if (sum(ok) >= 4)
    suppressWarnings(cor(donors_ss$aCAI[ok], donors_mz$aCAI[ok],
                         method = "spearman")) else NA_real_

  audit_program_rows[[pid]] <- data.table(
    program_id = pid, program_label = prog$program_label[i],
    program_order = prog$program_order[i],
    baseline_dir = bd, n_LE = prog$n_LE[i],
    A1_mean_aCAIrel_healthy_ss = mean(acai_rel_ss, na.rm = TRUE),
    A1_sd_aCAIrel_healthy_ss   = sd(acai_rel_ss, na.rm = TRUE),
    A2_frac_donors_aCAIpos_ss  = frac_pos_ss,
    A3_singscore_vs_meanZ_rho  = rho_methods,
    bulk_mag_singscore = bulk_mag_ss,
    bulk_mag_meanZ = bulk_mag_mz
  )
}
prog_audit <- rbindlist(audit_program_rows, use.names = TRUE)
fwrite(prog_audit, file.path(RES, "methodology_audit_program.tsv"),
       sep = "\t")
cat("\n--- per-program audit ---\n")
print(prog_audit[, .(program_label,
                     A1 = round(A1_mean_aCAIrel_healthy_ss, 2),
                     A1_sd = round(A1_sd_aCAIrel_healthy_ss, 2),
                     A2_dir = round(A2_frac_donors_aCAIpos_ss, 2),
                     A3_methods = round(A3_singscore_vs_meanZ_rho, 2))])

# ---- A4: cross-method robustness on Fig 2 / Fig 6A axis pooled g -----------
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

reml_pool <- function(g_vec, v_vec) {
  ok <- !is.na(g_vec) & !is.na(v_vec) & v_vec > 0
  if (sum(ok) < 2) return(list(g = NA_real_, I2 = NA_real_))
  fit <- tryCatch(rma(yi = g_vec[ok], vi = v_vec[ok],
                      method = "REML", test = "knha"),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(g = NA_real_, I2 = NA_real_))
  list(g = as.numeric(fit$beta), I2 = as.numeric(fit$I2))
}

axis_g_one_method <- function(scoring_fn) {
  # scoring_fn(M, rk, gene_set) → named numeric vector of per-sample scores
  sets <- fread(file.path(RES, "figure5_program_sets.tsv"))
  axes <- list(STRUCTURAL = sets[set_name == "STRUCTURAL", program_id],
               METABOLIC  = sets[set_name == "METABOLIC",  program_id])

  # Bulk magnitude per program (healthy GSE112339)
  bulk_mag <- function(le, bd) {
    s <- scoring_fn(mats[["GSE112339"]]$M, mats[["GSE112339"]]$rk, le)
    bd_d <- acai_donors(s, mats[["GSE112339"]]$meta, bd,
                        sample_col = "sample_id",
                        donor_col = "donor",
                        chamber_col = "chamber")
    abs(mean(bd_d$aCAI, na.rm = TRUE))
  }

  per_donor <- list()
  for (cohort in COHORTS) {
    md <- as.data.table(mats[[cohort]]$meta)
    keep_donor <- md[, .(rhythm = unique(rhythm)[1]), by = donor_id]
    keep_donor <- keep_donor[!is.na(rhythm)]
    rec <- list()
    for (i in seq_len(nrow(prog))) {
      pid <- prog$program_id[i]; bd <- prog$baseline_direction[i]
      le <- prog$le_genes[[i]]
      bm <- bulk_mag(le, bd)
      if (is.na(bm) || bm == 0) next
      sc <- scoring_fn(mats[[cohort]]$M, mats[[cohort]]$rk, le)
      bdat <- acai_donors(sc, mats[[cohort]]$meta, bd,
                          sample_col = "gsm", donor_col = "donor_id",
                          chamber_col = "chamber",
                          rhythm_col = "rhythm")
      bdat[, program_id := pid]
      bdat[, aCAI_rel := aCAI / bm]
      rec[[pid]] <- bdat[, .(donor_id, rhythm, program_id, aCAI_rel)]
    }
    cohort_long <- rbindlist(rec, use.names = TRUE)
    cohort_long[, cohort := cohort]
    per_donor[[cohort]] <- cohort_long
  }
  long <- rbindlist(per_donor, use.names = TRUE)

  # Axis composite: equal-weight mean of aCAI_rel within axis support set
  rows <- list()
  for (axis_lab in names(axes)) {
    ids <- axes[[axis_lab]]
    sub <- long[program_id %in% ids,
                .(composite = mean(aCAI_rel, na.rm = TRUE)),
                by = .(cohort, donor_id, rhythm)]
    gs <- numeric(length(COHORTS)); vs <- numeric(length(COHORTS))
    for (k in seq_along(COHORTS)) {
      .coh <- COHORTS[k]
      ss <- sub[cohort == .coh]
      g <- hedges_g(ss$composite[ss$rhythm == "AF"],
                    ss$composite[ss$rhythm == "SR"])
      gs[k] <- g$g; vs[k] <- g$var
    }
    pool <- reml_pool(gs, vs)
    rows[[axis_lab]] <- data.table(axis = axis_lab,
                                   pooled_g = pool$g, I2 = pool$I2)
  }
  rbindlist(rows, use.names = TRUE)
}

# wrap singscore vs meanZ as scoring functions
ss_fn  <- function(M, rk, le) singscore_one(rk, le)
mz_fn  <- function(M, rk, le) mean_z_one(M, le)

cat("\n=== A4: axis-level pooled g, singscore vs mean-z ===\n")
g_ss <- axis_g_one_method(ss_fn);  g_ss[, method := "singscore"]
g_mz <- axis_g_one_method(mz_fn);  g_mz[, method := "mean_z"]
g_compare <- rbind(g_ss, g_mz, use.names = TRUE)
fwrite(g_compare, file.path(RES, "methodology_audit_axis_g.tsv"), sep = "\t")
print(g_compare)

cat("\nH7 DONE.\n")
