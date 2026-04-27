#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# R3 - Program singscore on two pseudobulk panels.
#
# Two panels, one aggregator each:
#
#   GSE255612 (from R2)
#     expr = mean log-normalized expression per gene per (donor, compartment).
#     Already on a log scale + rank-preserving within sample;
#     rankGenes on raw matrix, then singscore.
#
#   Fig 1 pseudobulk (reuse: pseudobulk/<ct>/counts.tsv)
#     expr = raw integer UMI counts per gene per (donor, chamber, ct).
#     Normalise via edgeR::cpm(... log = TRUE, prior.count = 3) to match
#     Fig 2's voom-like log-CPM handling, then rankGenes + singscore.
#
# Program panel (frozen leading-edge unions): union of
#     (strict consensus != NA) U (115574 strict in {Amp, Att, Rev}) U
#     (any discovery trend_class != Preserved)
# with selection_reason recorded per program.
#
# Outputs:
#   results/figure3_program_panel.tsv
#   Project_P3/GSE255612/program_scores.tsv           (donor x ct x program)
#   Project_P3/fig1_healthy_program_scores.tsv        (donor x chamber x ct x program)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(edgeR)
  library(singscore)
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
cat(sprintf("ROOT = %s\n", ROOT))


# --------------------------------------------------------------------
# 1. Load the 10 Fig 1 stable programs (shared LE gene sets)
# --------------------------------------------------------------------

programs <- fread(file.path(ROOT, "normal_chamber_programs.tsv"))
for (col in c("program_id", "program", "direction", "leading_edge_genes"))
  if (col %in% colnames(programs))
    programs[[col]] <- gsub('^"|"$', '', programs[[col]])

programs[, le_genes := lapply(strsplit(leading_edge_genes, ",", fixed = TRUE),
                              function(x) unique(trimws(x)))]
programs[, n_LE_total := vapply(le_genes, length, integer(1))]

# sort by signed NES descending (Fig 1C / Fig 2 row order)
setorder(programs, -NES)
programs[, program_order := seq_len(.N)]

short_label_map <- c(
  "HALLMARK__HALLMARK_OXIDATIVE_PHOSPHORYLATION"                              = "OXPHOS",
  "REACTOME__REACTOME_RESPIRATORY_ELECTRON_TRANSPORT"                         = "Respiratory ETC",
  "REACTOME__REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT" = "Aerobic Resp. & ETC",
  "REACTOME__REACTOME_MITOCHONDRIAL_TRANSLATION"                              = "Mito. Translation",
  "REACTOME__REACTOME_TRANSLATION"                                            = "Cytoplasmic Translation",
  "REACTOME__REACTOME_HEPARAN_SULFATE_HEPARIN_HS_GAG_METABOLISM"              = "HS-GAG Metabolism",
  "REACTOME__REACTOME_SIGNALING_BY_BMP"                                       = "BMP Signalling",
  "REACTOME__REACTOME_HS_GAG_BIOSYNTHESIS"                                    = "HS-GAG Biosynthesis",
  "REACTOME__REACTOME_DISEASES_ASSOCIATED_WITH_GLYCOSAMINOGLYCAN_METABOLISM"  = "ECM proteoglycans",
  "REACTOME__REACTOME_RUNX2_REGULATES_BONE_DEVELOPMENT"                       = "RUNX2 Bone Dev."
)
programs[, program_label := unname(short_label_map[program_id])]
stopifnot(!any(is.na(programs$program_label)))
programs[, baseline_direction := ifelse(direction == "LA-biased", 1L, -1L)]


# --------------------------------------------------------------------
# 2. Build Figure 3 program panel (selection_reason)
# --------------------------------------------------------------------

rewiring <- fread(file.path(ROOT, "results", "af_chamber_rewiring_programs.tsv"))

panel <- merge(
  programs[, .(program_id, program_label, program_order, baseline_direction, NES)],
  rewiring[, .(program_id,
               consensus_class,
               rewiring_class_GSE128188, rewiring_class_GSE79768,
               rewiring_class_GSE115574,
               trend_class_GSE128188, trend_class_GSE79768, trend_class_GSE115574)],
  by = "program_id", all.x = TRUE, sort = FALSE
)

# pandas writes pd.NA as empty string; fread reads it as "" not NA, so
# use nzchar() to test "is non-empty" (robust against both conventions).
is_set <- function(x) !is.na(x) & nzchar(x)
panel[, reason_strict := is_set(consensus_class)]
panel[, reason_115574 := is_set(rewiring_class_GSE115574) &
                         rewiring_class_GSE115574 %in% c("Amplified", "Attenuated", "Reversed")]
panel[, reason_trend  := (is_set(trend_class_GSE128188) &
                          trend_class_GSE128188 %in% c("Amplified", "Attenuated", "Reversed")) |
                         (is_set(trend_class_GSE79768) &
                          trend_class_GSE79768  %in% c("Amplified", "Attenuated", "Reversed"))]

panel[, selection_reason := mapply(function(s, r, t) {
  reasons <- c()
  if (isTRUE(s)) reasons <- c(reasons, "strict_consensus")
  if (isTRUE(r)) reasons <- c(reasons, "115574_strict")
  if (isTRUE(t)) reasons <- c(reasons, "discovery_trend")
  paste(reasons, collapse = ",")
}, reason_strict, reason_115574, reason_trend)]

panel <- panel[reason_strict | reason_115574 | reason_trend]
setorder(panel, program_order)

panel_out <- panel[, .(program_id, program_label, program_order, baseline_direction,
                       consensus_class,
                       trend_class_GSE128188, trend_class_GSE79768,
                       rewiring_class_GSE115574,
                       selection_reason)]
fwrite(panel_out, file.path(ROOT, "results", "figure3_program_panel.tsv"), sep = "\t")
cat(sprintf("panel: %d programs selected\n", nrow(panel)))
cat("  reasons count: strict=%d, 115574=%d, trend=%d\n" |>
      sprintf(sum(panel$reason_strict),
              sum(panel$reason_115574),
              sum(panel$reason_trend)))
for (i in seq_len(nrow(panel))) {
  cat(sprintf("    %2d  %-28s  %+d  [%s]\n",
              panel$program_order[i], panel$program_label[i],
              panel$baseline_direction[i], panel$selection_reason[i]))
}

# restrict programs to panel (Figure 3 scoring is only on the panel)
programs_panel <- programs[program_id %in% panel$program_id]
programs_panel <- programs_panel[order(program_order)]


# --------------------------------------------------------------------
# 3. Scoring on GSE255612 pseudobulks
# --------------------------------------------------------------------

score_one_matrix <- function(expr_mat, programs_panel, sample_meta,
                             id_cols, extra_cols = list()) {
  # expr_mat: genes x samples with rownames=symbols, colnames=sample_ids
  stopifnot(all(sample_meta[[id_cols[1]]] %in% colnames(expr_mat)))
  expr_mat <- expr_mat[, sample_meta[[id_cols[1]]], drop = FALSE]
  ranked <- rankGenes(expr_mat)

  out_list <- list()
  for (i in seq_len(nrow(programs_panel))) {
    le <- programs_panel$le_genes[[i]]
    have <- intersect(le, rownames(expr_mat))
    if (length(have) < 3) {
      next
    }
    res <- simpleScore(ranked, upSet = have, centerScore = TRUE)
    df <- data.table(
      sample_id          = rownames(res),
      program_id         = programs_panel$program_id[i],
      program_label      = programs_panel$program_label[i],
      program_order      = programs_panel$program_order[i],
      baseline_direction = programs_panel$baseline_direction[i],
      n_LE_total         = programs_panel$n_LE_total[i],
      n_LE_matched       = length(have),
      score              = res$TotalScore
    )
    out_list[[i]] <- df
  }
  scores <- rbindlist(out_list)

  # attach sample-level metadata
  sample_meta <- as.data.table(sample_meta)
  meta_keys <- c("sample_id", unlist(extra_cols), id_cols[-1])
  sample_meta[, sample_id := get(id_cols[1])]
  scores <- merge(scores, unique(sample_meta[, ..meta_keys]),
                  by = "sample_id", sort = FALSE)
  scores
}


gse255612_groups <- list(
  atrialCM                      = list(label = "atrialCM",                      supp = FALSE),
  fibroblast                    = list(label = "fibroblast",                    supp = FALSE),
  endothelial_or_vascular       = list(label = "endothelial_or_vascular",       supp = FALSE),
  immune                        = list(label = "immune",                        supp = FALSE),
  supp_pericytes                = list(label = "supp_pericytes",                supp = TRUE),
  supp_vascular_smooth_muscle   = list(label = "supp_vascular_smooth_muscle",   supp = TRUE)
)

all_gse_scores <- list()
for (g in names(gse255612_groups)) {
  gdir <- file.path(ROOT, "Project_P3", "GSE255612", "pseudobulk", g)
  expr <- read.delim(file.path(gdir, "expr.tsv"),
                     check.names = FALSE, row.names = 1)
  expr <- as.matrix(expr)
  samps <- fread(file.path(gdir, "samples.tsv"))

  # column name in expr is donor_id (per r2 choice); align to samples$donor_id
  expr <- expr[, samps$donor_id, drop = FALSE]

  # singscore expects no ties-of-inf or NAs
  stopifnot(all(is.finite(expr)))

  ranked <- rankGenes(expr)
  for (i in seq_len(nrow(programs_panel))) {
    le <- programs_panel$le_genes[[i]]
    have <- intersect(le, rownames(expr))
    if (length(have) < 3) next
    res <- simpleScore(ranked, upSet = have, centerScore = TRUE)
    df <- data.table(
      cohort             = "GSE255612",
      group              = g,
      is_supp            = gse255612_groups[[g]]$supp,
      donor_id           = samps$donor_id,
      biosample_id       = samps$biosample_id,
      rhythm             = samps$rhythm,
      sex                = samps$sex,
      n_cells            = samps$n_cells,
      program_id         = programs_panel$program_id[i],
      program_label      = programs_panel$program_label[i],
      program_order      = programs_panel$program_order[i],
      baseline_direction = programs_panel$baseline_direction[i],
      n_LE_total         = programs_panel$n_LE_total[i],
      n_LE_matched       = length(have),
      score              = res$TotalScore
    )
    all_gse_scores[[paste(g, i, sep = ":")]] <- df
  }
  cat(sprintf("  scored GSE255612/%s (%d donors, %d programs)\n",
              g, nrow(samps), nrow(programs_panel)))
}
gse_scores <- rbindlist(all_gse_scores)
fwrite(gse_scores,
       file.path(ROOT, "Project_P3", "GSE255612", "program_scores.tsv"),
       sep = "\t")
cat(sprintf("wrote Project_P3/GSE255612/program_scores.tsv  (%d rows)\n",
            nrow(gse_scores)))


# --------------------------------------------------------------------
# 4. Scoring on Fig 1 pseudobulks (reuse)
# --------------------------------------------------------------------

# Fig 1 pseudobulk layout uses the long cell-type name (atrial_cardiomyocyte),
# which we map to the Fig 3 compartment4 key.
fig1_groups <- list(
  atrial_cardiomyocyte    = "atrialCM",
  fibroblast              = "fibroblast",
  endothelial_or_vascular = "endothelial_or_vascular",
  immune                  = "immune"
)

all_fig1_scores <- list()
for (g in names(fig1_groups)) {
  gdir <- file.path(ROOT, "pseudobulk", g)
  counts <- read.delim(file.path(gdir, "counts.tsv"),
                       check.names = FALSE, row.names = 1)
  counts <- as.matrix(counts)
  samps <- fread(file.path(gdir, "samples.tsv"))

  # Normalise raw counts to log-CPM.  edgeR requires non-zero library sizes.
  stopifnot(all(colSums(counts) > 0))
  logcpm <- edgeR::cpm(counts, log = TRUE, prior.count = 3)

  ranked <- rankGenes(logcpm)
  for (i in seq_len(nrow(programs_panel))) {
    le <- programs_panel$le_genes[[i]]
    have <- intersect(le, rownames(logcpm))
    if (length(have) < 3) next
    res <- simpleScore(ranked, upSet = have, centerScore = TRUE)

    # Match sample_id in res (=colname of logcpm) to samps$sample_id
    m <- match(rownames(res), samps$sample_id)
    stopifnot(!anyNA(m))
    df <- data.table(
      cohort             = "GSE112339",    # Fig 1 healthy single-cell source
      compartment4       = fig1_groups[[g]],
      sample_id          = rownames(res),
      donor_id           = samps$donor[m],
      chamber            = samps$chamber[m],
      n_cells            = samps$n_cells[m],
      program_id         = programs_panel$program_id[i],
      program_label      = programs_panel$program_label[i],
      program_order      = programs_panel$program_order[i],
      baseline_direction = programs_panel$baseline_direction[i],
      n_LE_total         = programs_panel$n_LE_total[i],
      n_LE_matched       = length(have),
      score              = res$TotalScore
    )
    all_fig1_scores[[paste(g, i, sep = ":")]] <- df
  }
  cat(sprintf("  scored Fig1 pseudobulk %s -> %s (%d samples, %d programs)\n",
              g, fig1_groups[[g]], nrow(samps), nrow(programs_panel)))
}
fig1_scores <- rbindlist(all_fig1_scores)
fwrite(fig1_scores,
       file.path(ROOT, "Project_P3", "fig1_healthy_program_scores.tsv"),
       sep = "\t")
cat(sprintf("wrote Project_P3/fig1_healthy_program_scores.tsv  (%d rows)\n",
            nrow(fig1_scores)))

cat("\nR3 DONE.\n")
