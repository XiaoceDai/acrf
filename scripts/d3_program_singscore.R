#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# D3 - Panel-program singscore on GSE41177 LAA + LAPV.
#
# Uses the 10-program panel inherited from Figure 3
# (`results/figure3_program_panel.tsv`) - same LE gene sets as Fig 2/3.
# Score = unidirectional singscore on RMA log2 expression (from D2).
#
# Output:
#   Project_P4/GSE41177/program_scores.tsv
#
# Columns: gsm, cohort, donor_id, rhythm, tissue, age, sex,
#          duration_raw, duration_months, duration_sensitivity_A,
#          duration_tier_B, keep_for_primary,
#          program_id, program_label, program_order, baseline_direction,
#          n_LE_total, n_LE_matched, score
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
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
IN_DIR <- file.path(ROOT, "Project_P4", "GSE41177")

# --------------------------------------------------------------------
# 1. Load panel + LE gene sets
# --------------------------------------------------------------------
programs <- fread(file.path(ROOT, "normal_chamber_programs.tsv"))
for (col in c("program_id", "program", "direction", "leading_edge_genes"))
  if (col %in% colnames(programs))
    programs[[col]] <- gsub('^"|"$', '', programs[[col]])
programs[, le_genes := lapply(strsplit(leading_edge_genes, ",", fixed = TRUE),
                              function(x) unique(trimws(x)))]
programs[, n_LE_total := vapply(le_genes, length, integer(1))]

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
setorder(programs, -NES)
programs[, program_order := seq_len(.N)]

# --------------------------------------------------------------------
# 2. Score each tissue matrix
# --------------------------------------------------------------------

meta <- fread(file.path(IN_DIR, "sample_metadata.tsv"))

score_tissue <- function(tissue_name) {
  expr_path <- file.path(IN_DIR, sprintf("expr_symbol_%s.tsv", tissue_name))
  expr <- as.matrix(read.delim(expr_path, check.names = FALSE, row.names = 1))
  cat(sprintf("\n=== %s : %d genes x %d samples ===\n",
              tissue_name, nrow(expr), ncol(expr)))

  ranked <- rankGenes(expr)
  rows <- list()
  for (i in seq_len(nrow(programs))) {
    le   <- programs$le_genes[[i]]
    have <- intersect(le, rownames(expr))
    if (length(have) < 3) next
    res <- simpleScore(ranked, upSet = have, centerScore = TRUE)
    # samples are the columns of expr; match to meta
    for (gsm_id in rownames(res)) {
      meta_row <- meta[meta$gsm == gsm_id]
      stopifnot(nrow(meta_row) == 1)
      rows[[length(rows) + 1]] <- data.table(
        gsm                    = gsm_id,
        cohort                 = meta_row$cohort,
        donor_id               = meta_row$donor_id,
        rhythm                 = meta_row$rhythm,
        tissue                 = meta_row$tissue,
        age                    = meta_row$age,
        sex                    = meta_row$sex,
        duration_raw           = meta_row$duration_raw,
        duration_months        = meta_row$duration_months,
        duration_sensitivity_A = meta_row$duration_sensitivity_A,
        duration_tier_B        = meta_row$duration_tier_B,
        keep_for_primary       = meta_row$keep_for_primary,
        program_id             = programs$program_id[i],
        program_label          = programs$program_label[i],
        program_order          = programs$program_order[i],
        baseline_direction     = programs$baseline_direction[i],
        n_LE_total             = programs$n_LE_total[i],
        n_LE_matched           = length(have),
        score                  = res[gsm_id, "TotalScore"]
      )
    }
  }
  rbindlist(rows)
}

all_scores <- rbind(score_tissue("LAA"), score_tissue("LAPV"))
out_path <- file.path(IN_DIR, "program_scores.tsv")
fwrite(all_scores, out_path, sep = "\t")
cat(sprintf("\nwrote %s  (%d rows = %d samples x %d programs)\n",
            sub(paste0("^", ROOT, "/"), "", out_path),
            nrow(all_scores),
            length(unique(all_scores$gsm)),
            length(unique(all_scores$program_id))))

# --------------------------------------------------------------------
# 3. Stdout preview: LE coverage + mean score per (rhythm, tissue)
# --------------------------------------------------------------------

cat("\n=== LE coverage check ===\n")
cov <- all_scores[, .(n_LE_matched = first(n_LE_matched),
                      n_LE_total   = first(n_LE_total),
                      le_coverage  = first(n_LE_matched)/first(n_LE_total)),
                  by = .(tissue, program_id, program_label)]
print(as.data.frame(cov[tissue == "LAA"]), row.names = FALSE)

cat("\nD3 DONE.\n")
