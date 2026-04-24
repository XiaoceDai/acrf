#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# P3 — Unidirectional singscore per sample, per program + Fig 1 bulk aCAI
#       reference.
#
# For each cohort and each of the 10 stable chamber programs defined in
# `normal_chamber_programs.tsv`, compute one singscore per sample using the
# program's Figure 1 leading-edge gene set.  Rank-based, within-sample
# normalisation is what makes scores comparable across platforms (voom
# logCPM vs RMA log2) without cross-cohort merging.
#
# Also computes Figure 1's healthy chamber bias on GSE112339 using the same
# singscore pipeline, which all downstream cohort aCAI will be normalised
# against (`aCAI_rel = aCAI / |aCAI_bulk|`).  This removes the LE-set-size
# bias that plagues absolute singscore magnitudes (OXPHOS with 120 genes
# has much smaller absolute aCAI than BMP Signalling with 6 genes).
#
# Spec: FIGURE2_WORKFLOW.md §3.1 / §3.1b / §3.2.
#
# Inputs
# ------
#   normal_chamber_programs.tsv                        Fig 1 prior
#   bulk/logcpm_voom.tsv + bulk/sample_meta.tsv        Fig 1 bulk source (GSE112339)
#   Project_P1/<cohort>/expr_symbol.tsv                from P2
#   Project_P1/<cohort>/sample_metadata.tsv            from P1
#
# Outputs
# -------
#   Project_P1/<cohort>/program_scores.tsv             wide: gsm x 10 programs (+ meta)
#   Project_P1/<cohort>/program_coverage.tsv           per-program LE coverage in this cohort
#   Project_P1/bulk_reference_aCAI.tsv                 per-program bulk aCAI normaliser
#   Project_P1/program_scores_long.tsv                 union, long format, incl. aCAI_rel
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(singscore)
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
  cwd <- getwd()
  if (basename(cwd) == "scripts") return(dirname(cwd))
  cwd
}
ROOT <- find_root()
P1   <- file.path(ROOT, "Project_P1")
cat(sprintf("ROOT = %s\n", ROOT))


# --------------------------------------------------------------------
# load Fig 1 prior
# --------------------------------------------------------------------

load_programs <- function() {
  df <- read.delim(file.path(ROOT, "normal_chamber_programs.tsv"),
                   sep = "\t", stringsAsFactors = FALSE)
  # Strip residual quoting left over from the P01 export.
  for (col in c("program_id", "program", "direction", "leading_edge_genes"))
    if (col %in% colnames(df)) df[[col]] <- gsub('^"|"$', '', df[[col]])

  # Lock the 10-program order to the same signed-NES order used by Fig 1C:
  # LA-biased (positive NES) top to bottom by descending NES, then RA-biased
  # (negative NES) top to bottom by descending NES (i.e. least negative first
  # is not what 1C does — 1C sorts by signed NES descending overall, so the
  # most negative sits at the bottom).
  df <- df[order(-df$NES), , drop = FALSE]
  rownames(df) <- NULL
  stopifnot(nrow(df) == 10)

  # Split the pipe/comma-separated LE gene string into a character vector.
  df$le_genes <- strsplit(df$leading_edge_genes, ",", fixed = TRUE)
  df$le_genes <- lapply(df$le_genes, function(x) unique(trimws(x)))
  df$n_LE_total <- vapply(df$le_genes, length, integer(1))

  # Short label map: hand-curated to match FIGURE2_WORKFLOW.md §5.
  short_label_map <- c(
    "HALLMARK__HALLMARK_OXIDATIVE_PHOSPHORYLATION"                             = "OXPHOS",
    "REACTOME__REACTOME_RESPIRATORY_ELECTRON_TRANSPORT"                        = "Respiratory ETC",
    "REACTOME__REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT" = "Aerobic Resp. & ETC",
    "REACTOME__REACTOME_MITOCHONDRIAL_TRANSLATION"                             = "Mito. Translation",
    "REACTOME__REACTOME_TRANSLATION"                                           = "Cytoplasmic Translation",
    "REACTOME__REACTOME_HEPARAN_SULFATE_HEPARIN_HS_GAG_METABOLISM"             = "HS-GAG Metabolism",
    "REACTOME__REACTOME_SIGNALING_BY_BMP"                                      = "BMP Signalling",
    "REACTOME__REACTOME_HS_GAG_BIOSYNTHESIS"                                   = "HS-GAG Biosynthesis",
    "REACTOME__REACTOME_DISEASES_ASSOCIATED_WITH_GLYCOSAMINOGLYCAN_METABOLISM" = "ECM proteoglycans",
    "REACTOME__REACTOME_RUNX2_REGULATES_BONE_DEVELOPMENT"                      = "RUNX2 Bone Dev."
  )
  df$program_label <- unname(short_label_map[df$program_id])
  stopifnot(!any(is.na(df$program_label)))
  df$program_order      <- seq_len(nrow(df))
  df$baseline_direction <- ifelse(df$direction == "LA-biased", 1L, -1L)
  df
}


# --------------------------------------------------------------------
# Figure 1 bulk aCAI reference (GSE112339, 8 paired donors)
# --------------------------------------------------------------------

score_bulk_reference <- function(programs) {
  cat("\n=== Fig 1 bulk reference (GSE112339) ===\n")
  expr_path <- file.path(ROOT, "bulk", "logcpm_voom.tsv")
  meta_path <- file.path(ROOT, "bulk", "sample_meta.tsv")
  raw <- read.delim(expr_path, check.names = FALSE, stringsAsFactors = FALSE)
  # columns: gene_key, gene_symbol, <sample1>, <sample2>, ...
  stopifnot(all(c("gene_key", "gene_symbol") %in% colnames(raw)))
  sym <- raw$gene_symbol
  # If multiple rows share a symbol (shouldn't in a voom matrix but be safe),
  # pick the row with the largest variance across samples.
  samp_cols <- setdiff(colnames(raw), c("gene_key", "gene_symbol"))
  mat_full  <- as.matrix(raw[, samp_cols, drop = FALSE])
  rownames(mat_full) <- sym
  if (anyDuplicated(sym)) {
    v <- apply(mat_full, 1, var, na.rm = TRUE)
    keep <- tapply(seq_along(sym), sym, function(idx) idx[which.max(v[idx])])
    mat_full <- mat_full[unlist(keep, use.names = FALSE), , drop = FALSE]
  }
  m <- mat_full

  meta <- read.delim(meta_path, stringsAsFactors = FALSE)
  # strip surrounding quotes carried over from the export
  for (col in c("sample_id", "chamber", "donor"))
    if (col %in% colnames(meta)) meta[[col]] <- gsub('^"|"$', '', meta[[col]])
  meta <- meta[match(colnames(m), meta$sample_id), , drop = FALSE]
  stopifnot(all(meta$sample_id == colnames(m)))
  cat(sprintf("  bulk expr: %d genes x %d samples, %d donors\n",
              nrow(m), ncol(m), length(unique(meta$donor))))

  ranked <- rankGenes(m)

  out <- data.frame(
    program_id         = programs$program_id,
    program_label      = programs$program_label,
    program_order      = programs$program_order,
    baseline_direction = programs$baseline_direction,
    n_LE_total         = programs$n_LE_total,
    n_LE_in_bulk       = NA_integer_,
    le_coverage_bulk   = NA_real_,
    n_donors_bulk      = NA_integer_,
    mean_CAI_bulk      = NA_real_,
    mean_aCAI_bulk     = NA_real_,
    sd_aCAI_bulk       = NA_real_,
    stringsAsFactors   = FALSE
  )

  for (i in seq_len(nrow(programs))) {
    le   <- programs$le_genes[[i]]
    have <- intersect(le, rownames(m))
    out$n_LE_in_bulk[i]    <- length(have)
    out$le_coverage_bulk[i] <- length(have) / length(le)
    if (length(have) < 3) next

    res <- simpleScore(ranked, upSet = have, centerScore = TRUE)
    sc  <- res$TotalScore
    names(sc) <- rownames(res)

    cai <- tapply(seq_along(sc), meta$donor, function(idx) {
      la <- sc[meta$sample_id[idx][meta$chamber[idx] == "LA"]]
      ra <- sc[meta$sample_id[idx][meta$chamber[idx] == "RA"]]
      if (length(la) == 1 && length(ra) == 1) la - ra else NA_real_
    })
    cai <- cai[!is.na(cai)]
    aCAI <- programs$baseline_direction[i] * cai
    out$n_donors_bulk[i]   <- length(aCAI)
    out$mean_CAI_bulk[i]   <- mean(cai)
    out$mean_aCAI_bulk[i]  <- mean(aCAI)
    out$sd_aCAI_bulk[i]    <- sd(aCAI)
  }
  out$abs_aCAI_bulk <- abs(out$mean_aCAI_bulk)

  ref_path <- file.path(P1, "bulk_reference_aCAI.tsv")
  write.table(out, ref_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  wrote %s\n", ref_path))
  for (i in seq_len(nrow(out))) {
    cat(sprintf("    %-25s dir=%+d  aCAI_bulk=%+.4f  (n_LE=%d/%d, n_donors=%d)\n",
                out$program_label[i], out$baseline_direction[i],
                out$mean_aCAI_bulk[i],
                out$n_LE_in_bulk[i], out$n_LE_total[i], out$n_donors_bulk[i]))
  }
  out
}


# --------------------------------------------------------------------
# singscore per cohort
# --------------------------------------------------------------------

score_cohort <- function(cohort, programs, bulk_ref) {
  cat(sprintf("\n=== %s ===\n", cohort))
  expr_path <- file.path(P1, cohort, "expr_symbol.tsv")
  meta_path <- file.path(P1, cohort, "sample_metadata.tsv")
  m <- as.matrix(read.delim(expr_path, row.names = 1, check.names = FALSE))
  meta <- read.delim(meta_path, stringsAsFactors = FALSE)
  meta <- meta[match(colnames(m), meta$gsm), , drop = FALSE]
  stopifnot(all(meta$gsm == colnames(m)))
  cat(sprintf("  expr: %d genes x %d samples\n", nrow(m), ncol(m)))

  # rankGenes handles ties robustly; we feed the whole matrix once.
  ranked <- rankGenes(m)

  score_mat <- matrix(
    NA_real_, nrow = ncol(m), ncol = nrow(programs),
    dimnames = list(colnames(m), programs$program_id)
  )
  coverage <- data.frame(
    program_id       = programs$program_id,
    program_label    = programs$program_label,
    program_order    = programs$program_order,
    baseline_direction = programs$baseline_direction,
    n_LE_total       = programs$n_LE_total,
    n_LE_in_matrix   = NA_integer_,
    le_coverage      = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(programs))) {
    le   <- programs$le_genes[[i]]
    have <- intersect(le, rownames(m))
    coverage$n_LE_in_matrix[i] <- length(have)
    coverage$le_coverage[i]    <- length(have) / length(le)

    if (length(have) < 3) {
      # singscore needs at least a handful of genes; record NA rather than
      # a meaningless score.
      warning(sprintf("  %s: only %d LE genes found in %s expression, skipping",
                      programs$program_id[i], length(have), cohort))
      next
    }
    # Unidirectional singscore (all LE genes are one-sided by construction),
    # default centerScore = TRUE so score is in ~[-0.5, 0.5].
    res <- simpleScore(ranked, upSet = have, centerScore = TRUE)
    score_mat[, i] <- res$TotalScore
  }

  # wide table: gsm + donor + chamber + rhythm + 10 score columns (by program_label)
  wide <- data.frame(
    gsm       = meta$gsm,
    cohort    = meta$cohort,
    donor_id  = meta$donor_id,
    chamber   = meta$chamber,
    rhythm    = meta$rhythm,
    is_paired = meta$is_paired,
    stringsAsFactors = FALSE
  )
  score_df <- as.data.frame(score_mat, stringsAsFactors = FALSE)
  colnames(score_df) <- programs$program_label
  wide <- cbind(wide, score_df)

  out_wide <- file.path(P1, cohort, "program_scores.tsv")
  write.table(wide, out_wide, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  wrote %s\n", out_wide))

  out_cov <- file.path(P1, cohort, "program_coverage.tsv")
  write.table(coverage, out_cov, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  wrote %s\n", out_cov))

  # build long table for union: one row per (donor, program) with CAI + aCAI + aCAI_rel.
  # score_mat is sample × program; collapse to donor × program within this cohort.
  long_list <- list()
  for (i in seq_len(nrow(programs))) {
    pid   <- programs$program_id[i]
    pdir  <- programs$baseline_direction[i]
    ab    <- bulk_ref$abs_aCAI_bulk[match(pid, bulk_ref$program_id)]
    for (d in unique(meta$donor_id)) {
      rows <- which(meta$donor_id == d)
      if (length(rows) != 2) next
      la_row <- rows[meta$chamber[rows] == "LA"]
      ra_row <- rows[meta$chamber[rows] == "RA"]
      if (length(la_row) != 1 || length(ra_row) != 1) next
      sc_LA <- score_mat[la_row, i]
      sc_RA <- score_mat[ra_row, i]
      cai   <- sc_LA - sc_RA
      aCAI  <- pdir * cai
      aCAI_rel <- if (!is.na(ab) && ab > 0) aCAI / ab else NA_real_

      long_list[[length(long_list) + 1]] <- data.frame(
        cohort             = cohort,
        donor_id           = d,
        rhythm             = meta$rhythm[rows[1]],
        gsm_LA             = meta$gsm[la_row],
        gsm_RA             = meta$gsm[ra_row],
        program_id         = pid,
        program_label      = programs$program_label[i],
        program_order      = programs$program_order[i],
        baseline_direction = pdir,
        score_LA           = sc_LA,
        score_RA           = sc_RA,
        CAI                = cai,
        aCAI               = aCAI,
        abs_aCAI_bulk      = ab,
        aCAI_rel           = aCAI_rel,
        stringsAsFactors   = FALSE
      )
    }
  }
  do.call(rbind, long_list)
}


# --------------------------------------------------------------------
# main
# --------------------------------------------------------------------

programs <- load_programs()
cat(sprintf("loaded %d programs (%d LA-biased, %d RA-biased)\n",
            nrow(programs), sum(programs$baseline_direction == 1),
            sum(programs$baseline_direction == -1)))
cat("program order (Fig 1C order = Fig 2 row order):\n")
for (i in seq_len(nrow(programs))) {
  cat(sprintf("  %2d  %-25s  %+d  NES=%+.2f  n_LE=%d\n",
              programs$program_order[i], programs$program_label[i],
              programs$baseline_direction[i], programs$NES[i],
              programs$n_LE_total[i]))
}

# Bulk reference first, because cohort scoring normalises against it.
bulk_ref <- score_bulk_reference(programs)

cohorts <- c("GSE128188", "GSE79768", "GSE115574")
all_long <- do.call(rbind, lapply(
  cohorts, function(c) score_cohort(c, programs = programs, bulk_ref = bulk_ref)
))

out_long <- file.path(P1, "program_scores_long.tsv")
write.table(all_long, out_long, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\nunion long table: %d rows -> %s\n", nrow(all_long),
            sub(paste0("^", ROOT, "/"), "", out_long)))

cat("\nDONE.\n")
