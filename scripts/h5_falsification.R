#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# H5 — Specificity / Falsification (Figure S8).
#
# Three components, all driven by the framework already used for Figs 2/4/6:
#   A  size-matched random gene-set permutations  (1000 / program)
#   B  axis-membership permutation                (1000)
#   C  negative-control pathways                  (6 pre-registered)
#
# Spec: FIGURE_S8_FALSIFICATION_SPEC.md
#
# Component selectable via CLI:   Rscript scripts/h5_falsification.R [A|B|C|all]
# Defaults: B then C then A   (cheapest first).
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(metafor); library(singscore)
})

args <- commandArgs(trailingOnly = TRUE)
do_run <- if (length(args)) args else c("B", "C", "A")
do_run <- toupper(do_run); if ("ALL" %in% do_run) do_run <- c("A", "B", "C")

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
# Shared utilities
# ---------------------------------------------------------------------------
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
  if (sum(ok) < 2) return(NA_real_)
  fit <- tryCatch(rma(yi = g_vec[ok], vi = v_vec[ok],
                      method = "REML", test = "knha"),
                  error = function(e) NULL)
  if (is.null(fit)) NA_real_ else as.numeric(fit$beta)
}

# ---------------------------------------------------------------------------
# Component B — axis-membership permutation
# ---------------------------------------------------------------------------
component_B <- function(B = 1000L, seed = 20260419L) {
  cat("\n=== Component B: axis-membership permutation ===\n")
  set.seed(seed)

  long <- fread(file.path(P1, "program_scores_long.tsv"))
  sets <- fread(file.path(ROOT, "results", "figure5_program_sets.tsv"))
  prog_ids <- sets[set_name %in% c("STRUCTURAL", "METABOLIC"), program_id]
  stopifnot(length(prog_ids) == 6L)

  # donor x program aCAI_rel matrix per cohort
  axis_pool_g <- function(cohort_donor, ids_struct, ids_metab) {
    out <- list()
    for (.coh in COHORTS) {
      sub <- cohort_donor[cohort == .coh]
      sub_w <- dcast(sub[program_id %in% prog_ids],
                     donor_id + rhythm ~ program_id, value.var = "aCAI_rel")
      if (nrow(sub_w) == 0) next
      m_struct <- rowMeans(as.matrix(sub_w[, ..ids_struct]),
                           na.rm = TRUE)
      m_metab  <- rowMeans(as.matrix(sub_w[, ..ids_metab]),
                           na.rm = TRUE)
      gA <- hedges_g(m_struct[sub_w$rhythm == "AF"],
                     m_struct[sub_w$rhythm == "SR"])
      gB <- hedges_g(m_metab[sub_w$rhythm == "AF"],
                     m_metab[sub_w$rhythm == "SR"])
      out[[.coh]] <- data.table(cohort = .coh,
                                gS = gA$g, vS = gA$var,
                                gM = gB$g, vM = gB$var)
    }
    cc <- rbindlist(out, use.names = TRUE)
    list(g_struct = reml_pool(cc$gS, cc$vS),
         g_metab  = reml_pool(cc$gM, cc$vM))
  }

  cohort_donor <- long[program_id %in% prog_ids,
                       .(cohort, donor_id, rhythm, program_id, aCAI_rel)]

  # real partition
  real_S <- sets[set_name == "STRUCTURAL", program_id]
  real_M <- sets[set_name == "METABOLIC",  program_id]
  real <- axis_pool_g(cohort_donor, real_S, real_M)
  cat(sprintf("  real:   g_STRUCT = %+.3f   g_METAB = %+.3f\n",
              real$g_struct, real$g_metab))

  # permutation
  gS_perm <- numeric(B); gM_perm <- numeric(B); gDiff_perm <- numeric(B)
  for (b in seq_len(B)) {
    s <- sample(prog_ids, 3L)
    m <- setdiff(prog_ids, s)
    res <- axis_pool_g(cohort_donor, s, m)
    gS_perm[b] <- res$g_struct; gM_perm[b] <- res$g_metab
    gDiff_perm[b] <- abs(res$g_struct - res$g_metab)
  }

  real_diff <- abs(real$g_struct - real$g_metab)
  p_S    <- (1 + sum(abs(gS_perm)    >= abs(real$g_struct), na.rm = TRUE)) / (1 + B)
  p_M    <- (1 + sum(abs(gM_perm)    >= abs(real$g_metab),  na.rm = TRUE)) / (1 + B)
  p_diff <- (1 + sum(gDiff_perm >= real_diff, na.rm = TRUE)) / (1 + B)
  cat(sprintf("  perm |g_STRUCT|: q025/q975 = %.2f / %.2f   p_emp = %.4f\n",
              quantile(abs(gS_perm), 0.025, na.rm = TRUE),
              quantile(abs(gS_perm), 0.975, na.rm = TRUE), p_S))
  cat(sprintf("  perm |g_METAB|:  q025/q975 = %.2f / %.2f   p_emp = %.4f\n",
              quantile(abs(gM_perm), 0.025, na.rm = TRUE),
              quantile(abs(gM_perm), 0.975, na.rm = TRUE), p_M))
  cat(sprintf("  perm |g_S - g_M|:q025/q975 = %.2f / %.2f   p_emp = %.4f\n",
              quantile(gDiff_perm, 0.025, na.rm = TRUE),
              quantile(gDiff_perm, 0.975, na.rm = TRUE), p_diff))

  out_summary <- data.table(
    axis = c("STRUCTURAL", "METABOLIC", "axis_separation"),
    real_g = c(real$g_struct, real$g_metab, real_diff),
    perm_q025 = c(quantile(abs(gS_perm), 0.025, na.rm = TRUE),
                  quantile(abs(gM_perm), 0.025, na.rm = TRUE),
                  quantile(gDiff_perm,    0.025, na.rm = TRUE)),
    perm_q975 = c(quantile(abs(gS_perm), 0.975, na.rm = TRUE),
                  quantile(abs(gM_perm), 0.975, na.rm = TRUE),
                  quantile(gDiff_perm,    0.975, na.rm = TRUE)),
    p_emp_two_sided = c(p_S, p_M, p_diff)
  )
  fwrite(out_summary, file.path(RESULTS, "figure_s8b_axis_perm.tsv"),
         sep = "\t")
  fwrite(data.table(b = seq_len(B), gS = gS_perm, gM = gM_perm,
                    gDiff = gDiff_perm),
         file.path(RESULTS, "figure_s8b_axis_perm_distribution.tsv.gz"),
         sep = "\t", compress = "gzip")
  cat("wrote results/figure_s8b_axis_perm.tsv\n")
}

# ---------------------------------------------------------------------------
# Helpers for Components A + C — score arbitrary gene sets
# ---------------------------------------------------------------------------
load_cohort_matrices <- function() {
  cat("\nloading cohort expression matrices ...\n")
  out <- list()
  # paired cohorts
  for (cohort in COHORTS) {
    expr_path <- file.path(P1, cohort, "expr_symbol.tsv")
    meta_path <- file.path(P1, cohort, "sample_metadata.tsv")
    expr <- fread(expr_path, sep = "\t", header = TRUE)
    sym <- expr[[1]]
    M <- as.matrix(expr[, -1, with = FALSE])
    rownames(M) <- sym
    if (anyDuplicated(sym)) {
      v <- apply(M, 1, var, na.rm = TRUE)
      keep <- tapply(seq_along(sym), sym, function(idx) idx[which.max(v[idx])])
      M <- M[unlist(keep, use.names = FALSE), , drop = FALSE]
    }
    meta <- fread(meta_path, sep = "\t", header = TRUE)
    out[[cohort]] <- list(M = M, meta = meta)
    cat(sprintf("  %s: %d genes x %d samples\n", cohort, nrow(M), ncol(M)))
  }
  # GSE112339 healthy bulk reference
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
  out[["GSE112339"]] <- list(M = Mb, meta = bulk_meta)
  cat(sprintf("  GSE112339 healthy bulk: %d genes x %d samples\n",
              nrow(Mb), ncol(Mb)))
  out
}

# Pre-rank matrices once (singscore::rankGenes) for fast re-scoring
prepare_ranks <- function(mats) {
  cat("ranking matrices ...\n")
  out <- list()
  for (k in names(mats)) {
    out[[k]] <- list(rk = rankGenes(mats[[k]]$M),
                     meta = mats[[k]]$meta)
  }
  out
}

score_one_set <- function(gene_set, ranks_obj) {
  hit <- intersect(gene_set, rownames(ranks_obj$rk))
  if (length(hit) < 3) return(rep(NA_real_, ncol(ranks_obj$rk)))
  s <- simpleScore(ranks_obj$rk, upSet = hit,
                   knownDirection = TRUE, centerScore = TRUE)
  setNames(s$TotalScore, rownames(s))
}

# Compute aCAI per donor in a paired cohort given per-sample scores
acai_per_donor <- function(scores, meta, baseline_dir) {
  meta <- as.data.table(meta)
  if (!all(c("donor_id", "chamber", "rhythm", "gsm") %in% names(meta)))
    return(data.table(donor_id = character(0), rhythm = character(0),
                      aCAI = numeric(0)))
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

# Compute |aCAI_bulk(p)| from the GSE112339 healthy bulk reference
acai_bulk_magnitude <- function(scores_bulk, meta_bulk) {
  meta_bulk <- as.data.table(meta_bulk)
  meta_bulk[, score := scores_bulk[match(sample_id, names(scores_bulk))]]
  by_d <- meta_bulk[!is.na(score), .(
    LA = mean(score[chamber == "LA"], na.rm = TRUE),
    RA = mean(score[chamber == "RA"], na.rm = TRUE)
  ), by = donor]
  abs(mean(by_d$LA - by_d$RA, na.rm = TRUE))
}

# Full pipeline: gene set + baseline_dir → pooled g across cohorts
score_set_pooled_g <- function(gene_set, baseline_dir, ranks_all) {
  s_bulk <- score_one_set(gene_set, ranks_all[["GSE112339"]])
  bulk_mag <- acai_bulk_magnitude(s_bulk, ranks_all[["GSE112339"]]$meta)
  if (is.na(bulk_mag) || bulk_mag == 0) return(NA_real_)
  gs <- numeric(length(COHORTS)); vs <- numeric(length(COHORTS))
  for (k in seq_along(COHORTS)) {
    cohort <- COHORTS[k]
    s <- score_one_set(gene_set, ranks_all[[cohort]])
    by_donor <- acai_per_donor(s, ranks_all[[cohort]]$meta, baseline_dir)
    if (nrow(by_donor) < 3) { gs[k] <- NA_real_; vs[k] <- NA_real_; next }
    by_donor[, aCAI_rel := aCAI / bulk_mag]
    g <- hedges_g(by_donor$aCAI_rel[by_donor$rhythm == "AF"],
                  by_donor$aCAI_rel[by_donor$rhythm == "SR"])
    gs[k] <- g$g; vs[k] <- g$var
  }
  reml_pool(gs, vs)
}

# ---------------------------------------------------------------------------
# Component A — size-matched random gene-set permutations
# ---------------------------------------------------------------------------
component_A <- function(B = 1000L, seed = 20260418L, mats = NULL) {
  cat("\n=== Component A: size-matched random gene-set permutations ===\n")
  if (is.null(mats)) mats <- load_cohort_matrices()
  ranks_all <- prepare_ranks(mats)

  # Universe = intersection of gene symbols across all 4 matrices
  univ <- Reduce(intersect, lapply(ranks_all, function(x) rownames(x$rk)))
  cat(sprintf("expressed gene universe: %d symbols\n", length(univ)))

  prog <- read.delim(file.path(ROOT, "normal_chamber_programs.tsv"),
                     sep = "\t", stringsAsFactors = FALSE)
  for (col in c("program_id", "leading_edge_genes"))
    if (col %in% colnames(prog)) prog[[col]] <- gsub('^"|"$', '', prog[[col]])
  sets <- fread(file.path(ROOT, "results", "figure5_program_sets.tsv"))
  panel <- sets[set_name %in% c("STRUCTURAL", "METABOLIC")]
  panel <- merge(panel,
                 prog[, c("program_id", "direction", "leading_edge_genes")],
                 by = "program_id", sort = FALSE)
  panel[, le_genes := lapply(leading_edge_genes,
                             function(s) unique(trimws(strsplit(s, ",", fixed = TRUE)[[1]])))]
  panel[, baseline_dir := ifelse(direction == "LA-biased", 1L, -1L)]

  results <- list()
  for (i in seq_len(nrow(panel))) {
    pid <- panel$program_id[i]; plab <- panel$program_label[i]
    le <- panel$le_genes[[i]]; bdir <- panel$baseline_dir[i]
    n_LE <- length(le)
    cat(sprintf("\n[%d/%d] %s  (n_LE = %d)\n",
                i, nrow(panel), plab, n_LE))

    real_g <- score_set_pooled_g(le, bdir, ranks_all)
    cat(sprintf("  real pooled g = %+.3f\n", real_g))

    set.seed(seed + i)
    perm_g <- numeric(B)
    t0 <- Sys.time()
    for (b in seq_len(B)) {
      rs <- sample(univ, n_LE, replace = FALSE)
      perm_g[b] <- score_set_pooled_g(rs, bdir, ranks_all)
      if (b %% 100 == 0) cat(sprintf("    b=%d  (%.1fs)\n", b,
                                     as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
    p_emp <- (1 + sum(abs(perm_g) >= abs(real_g), na.rm = TRUE)) / (1 + B)
    results[[pid]] <- data.table(
      program_id = pid, program_label = plab,
      n_LE = n_LE, real_g = real_g,
      perm_q025 = quantile(perm_g, 0.025, na.rm = TRUE),
      perm_q975 = quantile(perm_g, 0.975, na.rm = TRUE),
      perm_med  = median(perm_g, na.rm = TRUE),
      p_emp_two_sided = p_emp,
      n_perm = B, n_valid = sum(!is.na(perm_g))
    )
    fwrite(data.table(b = seq_len(B), g = perm_g)[!is.na(g)],
           file.path(RESULTS,
                     sprintf("figure_s8a_perm_dist_%s.tsv.gz", pid)),
           sep = "\t", compress = "gzip")
  }
  out <- rbindlist(results, use.names = TRUE)
  fwrite(out, file.path(RESULTS, "figure_s8a_size_matched.tsv"), sep = "\t")
  cat("\n--- Component A summary ---\n"); print(out)
}

# ---------------------------------------------------------------------------
# Component C — negative-control pathways
# ---------------------------------------------------------------------------
load_gmt <- function(path) {
  L <- readLines(path)
  out <- list()
  for (l in L) {
    parts <- strsplit(l, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 3) next
    name <- parts[1]
    out[[name]] <- parts[3:length(parts)]
  }
  out
}

component_C <- function(mats = NULL) {
  cat("\n=== Component C: negative-control pathways ===\n")
  if (is.null(mats)) mats <- load_cohort_matrices()
  ranks_all <- prepare_ranks(mats)

  hall <- load_gmt(file.path(ROOT, "h.all.v2026.1.Hs.symbols.gmt"))
  reac <- load_gmt(file.path(ROOT, "c2.cp.reactome.v2026.1.Hs.symbols.gmt"))

  controls <- list(
    list(id = "HALLMARK_INFLAMMATORY_RESPONSE",      cat = "axis-specificity",
         genes = hall[["HALLMARK_INFLAMMATORY_RESPONSE"]]),
    list(id = "HALLMARK_KRAS_SIGNALING_UP",          cat = "axis-specificity",
         genes = hall[["HALLMARK_KRAS_SIGNALING_UP"]]),
    list(id = "HALLMARK_E2F_TARGETS",                cat = "axis-specificity",
         genes = hall[["HALLMARK_E2F_TARGETS"]]),
    list(id = "REACTOME_OLFACTORY_SIGNALING_PATHWAY", cat = "unrelated",
         genes = reac[["REACTOME_OLFACTORY_SIGNALING_PATHWAY"]]),
    list(id = "REACTOME_KERATINIZATION",              cat = "unrelated",
         genes = reac[["REACTOME_KERATINIZATION"]]),
    list(id = "REACTOME_SENSORY_PERCEPTION_OF_TASTE", cat = "unrelated",
         genes = reac[["REACTOME_SENSORY_PERCEPTION_OF_TASTE"]])
  )
  for (k in seq_along(controls)) {
    g <- controls[[k]]$genes
    if (is.null(g)) {
      cat(sprintf("  [WARN] %s not in GMT, skipping\n", controls[[k]]$id))
      next
    }
    cat(sprintf("  %s  n=%d\n", controls[[k]]$id, length(g)))
  }

  # Healthy-determined baseline direction (sign of mean LA-RA in GSE112339)
  # then a single pooled-g per pathway (not "best of") + I² + LAA ρ.
  healthy_dir <- function(gene_set) {
    s <- score_one_set(gene_set, ranks_all[["GSE112339"]])
    md <- as.data.table(ranks_all[["GSE112339"]]$meta)
    md[, score := s[match(sample_id, names(s))]]
    by_d <- md[!is.na(score), .(
      LA = mean(score[chamber == "LA"], na.rm = TRUE),
      RA = mean(score[chamber == "RA"], na.rm = TRUE)
    ), by = donor]
    if (!nrow(by_d)) return(list(dir = NA_integer_, mag = NA_real_))
    md_diff <- mean(by_d$LA - by_d$RA, na.rm = TRUE)
    list(dir = if (is.na(md_diff)) NA_integer_ else
                  ifelse(md_diff >= 0, +1L, -1L),
         mag = abs(md_diff))
  }

  # extend score_set_pooled_g to also return per-cohort g/var + pooled I²
  full_axis <- function(gene_set, bdir, ranks_all) {
    s_bulk <- score_one_set(gene_set, ranks_all[["GSE112339"]])
    bulk_mag <- acai_bulk_magnitude(s_bulk, ranks_all[["GSE112339"]]$meta)
    if (is.na(bulk_mag) || bulk_mag == 0)
      return(list(g = NA_real_, I2 = NA_real_,
                  per_cohort = data.table()))
    rows <- list()
    for (k in seq_along(COHORTS)) {
      cohort <- COHORTS[k]
      sc <- score_one_set(gene_set, ranks_all[[cohort]])
      bd <- acai_per_donor(sc, ranks_all[[cohort]]$meta, bdir)
      if (nrow(bd) < 3) next
      bd[, aCAI_rel := aCAI / bulk_mag]
      gg <- hedges_g(bd$aCAI_rel[bd$rhythm == "AF"],
                     bd$aCAI_rel[bd$rhythm == "SR"])
      rows[[cohort]] <- data.table(cohort = cohort, g = gg$g, var = gg$var)
    }
    cc <- rbindlist(rows, use.names = TRUE)
    if (nrow(cc) < 2) return(list(g = NA_real_, I2 = NA_real_,
                                  per_cohort = cc))
    fit <- tryCatch(rma(yi = cc$g, vi = cc$var, method = "REML",
                        test = "knha"),
                    error = function(e) NULL)
    if (is.null(fit)) return(list(g = NA_real_, I2 = NA_real_,
                                  per_cohort = cc))
    list(g = as.numeric(fit$beta), I2 = as.numeric(fit$I2),
         per_cohort = cc)
  }

  # LAA duration ρ from Fig 4 framework (use existing program_scores if present)
  # Otherwise compute fresh: score on GSE41177 LAA, Spearman with duration.
  af_dur_rho <- function(gene_set, ranks_all) {
    # use Project_P4/GSE41177 expr if available
    laa <- file.path(ROOT, "Project_P4", "GSE41177", "expr_symbol_LAA.tsv")
    meta_p <- file.path(ROOT, "Project_P4", "GSE41177", "sample_metadata.tsv")
    if (!file.exists(laa) || !file.exists(meta_p)) return(NA_real_)
    expr <- fread(laa, sep = "\t", header = TRUE)
    sym <- expr[[1]]; M <- as.matrix(expr[, -1, with = FALSE])
    rownames(M) <- sym
    if (anyDuplicated(sym)) {
      v <- apply(M, 1, var, na.rm = TRUE)
      keep <- tapply(seq_along(sym), sym, function(idx) idx[which.max(v[idx])])
      M <- M[unlist(keep, use.names = FALSE), , drop = FALSE]
    }
    rk <- rankGenes(M)
    hit <- intersect(gene_set, rownames(rk))
    if (length(hit) < 3) return(NA_real_)
    s <- simpleScore(rk, upSet = hit, knownDirection = TRUE,
                     centerScore = TRUE)
    sc <- s$TotalScore; names(sc) <- rownames(s)
    md <- fread(meta_p)
    md <- md[tissue == "LAA" & rhythm == "AF" & !is.na(duration_months)
             & keep_for_primary == TRUE]
    md[, score := sc[match(gsm, names(sc))]]
    md <- md[!is.na(score)]
    if (nrow(md) < 6) return(NA_real_)
    suppressWarnings(cor(md$score, md$duration_months, method = "spearman"))
  }

  results <- list()
  for (k in seq_along(controls)) {
    cc <- controls[[k]]
    if (is.null(cc$genes)) next
    hd <- healthy_dir(cc$genes)
    fa <- full_axis(cc$genes, ifelse(is.na(hd$dir), 1L, hd$dir), ranks_all)
    rho <- af_dur_rho(cc$genes, ranks_all)
    cat(sprintf("  %s  bdir=%+d  pooled g = %+.3f  I² = %.0f%%  ρ_LAA = %+.2f\n",
                cc$id, hd$dir, fa$g, fa$I2 * 100, rho))
    results[[cc$id]] <- data.table(
      pathway = cc$id, category = cc$cat, n_genes = length(cc$genes),
      healthy_dir = hd$dir, healthy_mag = hd$mag,
      pooled_g = fa$g, I2 = fa$I2,
      duration_rho_LAA = rho
    )
  }
  out <- rbindlist(results, use.names = TRUE)
  fwrite(out, file.path(RESULTS, "figure_s8c_negative_controls.tsv"),
         sep = "\t")
  cat("\n--- Component C summary ---\n"); print(out)
}

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
mats_cache <- NULL
if ("B" %in% do_run) component_B()
if ("C" %in% do_run) {
  if (is.null(mats_cache)) mats_cache <- load_cohort_matrices()
  component_C(mats_cache)
}
if ("A" %in% do_run) {
  if (is.null(mats_cache)) mats_cache <- load_cohort_matrices()
  component_A(mats = mats_cache)
}
cat("\nH5 DONE.\n")
