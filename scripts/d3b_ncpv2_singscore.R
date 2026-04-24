#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# D3b - NCP_v2 program singscore on GSE41177 LAA + LAPV.
#
# Refreshes Fig 4 data against the 30-program NCP_v2 panel using the
# leading-edge gene sets re-derived in R3b
# (results/fig3_ncpv2_program_le.tsv).
#
# Output:
#   Project_P4/GSE41177/program_scores_v2.tsv
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
ROOT   <- find_root()
IN_DIR <- file.path(ROOT, "Project_P4", "GSE41177")

# --------------------------------------------------------------------
# 1. Load NCP_v2 panel + LE genes (R3b output) + axis families (Fig 1e)
# --------------------------------------------------------------------
panel <- fread(file.path(ROOT, "results", "fig3_ncpv2_program_le.tsv"))
panel[, le_genes := strsplit(leading_edge_genes, ",", fixed = TRUE)]
panel[, le_genes := lapply(le_genes, function(x) unique(trimws(x)))]
panel[, n_LE_total := vapply(le_genes, length, integer(1))]
setorder(panel, program_order)

# pretty short labels (mirror R3b)
.short_map <- c(
  "HALLMARK__HALLMARK_OXIDATIVE_PHOSPHORYLATION"                              = "OXPHOS",
  "REACTOME__REACTOME_RESPIRATORY_ELECTRON_TRANSPORT"                         = "ETC",
  "REACTOME__REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT" = "Aerobic Resp & ETC",
  "REACTOME__REACTOME_MITOCHONDRIAL_TRANSLATION"                              = "Mito translation",
  "REACTOME__REACTOME_TRANSLATION"                                            = "Translation (cyto)",
  "REACTOME__REACTOME_HEPARAN_SULFATE_HEPARIN_HS_GAG_METABOLISM"              = "HS-GAG metabolism",
  "REACTOME__REACTOME_SIGNALING_BY_BMP"                                       = "BMP signalling",
  "REACTOME__REACTOME_HS_GAG_BIOSYNTHESIS"                                    = "HS-GAG biosynthesis",
  "REACTOME__REACTOME_DISEASES_ASSOCIATED_WITH_GLYCOSAMINOGLYCAN_METABOLISM"  = "GAG metabolism",
  "REACTOME__REACTOME_RUNX2_REGULATES_BONE_DEVELOPMENT"                       = "RUNX2 / bone dev"
)
.short <- function(s) {
  out <- unname(.short_map[s])
  fb  <- gsub("^[A-Z]+__[A-Z]+_+", "", s)
  fb  <- tolower(gsub("_+", " ", fb))
  ifelse(is.na(out), substr(fb, 1, 30), out)
}
panel[, program_label := .short(program_id)]


# --------------------------------------------------------------------
# 2. Score each tissue matrix (LAA + LAPV)
# --------------------------------------------------------------------
meta <- fread(file.path(IN_DIR, "sample_metadata.tsv"))

score_tissue <- function(tissue_name) {
  expr_path <- file.path(IN_DIR, sprintf("expr_symbol_%s.tsv", tissue_name))
  expr <- as.matrix(read.delim(expr_path, check.names = FALSE, row.names = 1))
  cat(sprintf("\n=== %s : %d genes x %d samples ===\n",
              tissue_name, nrow(expr), ncol(expr)))

  ranked <- rankGenes(expr)
  rows <- list()
  for (i in seq_len(nrow(panel))) {
    le   <- panel$le_genes[[i]]
    have <- intersect(le, rownames(expr))
    if (length(have) < 3) next
    res <- simpleScore(ranked, upSet = have, centerScore = TRUE)
    for (gsm_id in rownames(res)) {
      meta_row <- meta[meta$gsm == gsm_id]
      if (nrow(meta_row) != 1) next
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
        program_id             = panel$program_id[i],
        program_label          = panel$program_label[i],
        program_order          = panel$program_order[i],
        family_id              = panel$family_id[i],
        baseline_direction     = panel$baseline_direction[i],
        n_LE_total             = panel$n_LE_total[i],
        n_LE_matched           = length(have),
        score                  = res[gsm_id, "TotalScore"]
      )
    }
  }
  rbindlist(rows)
}

all_scores <- rbind(score_tissue("LAA"), score_tissue("LAPV"))
out_path   <- file.path(IN_DIR, "program_scores_v2.tsv")
fwrite(all_scores, out_path, sep = "\t")
cat(sprintf("\nwrote %s  (%d rows = %d samples x %d programs)\n",
            sub(paste0("^", ROOT, "/"), "", out_path),
            nrow(all_scores),
            length(unique(all_scores$gsm)),
            length(unique(all_scores$program_id))))

# LE coverage QC
cat("\n=== LE coverage check (LAA) ===\n")
cov <- all_scores[, .(n_LE_matched = first(n_LE_matched),
                       n_LE_total   = first(n_LE_total),
                       le_coverage  = round(first(n_LE_matched) /
                                              first(n_LE_total), 3)),
                   by = .(tissue, program_id, program_label, family_id)]
print(as.data.frame(cov[tissue == "LAA"]), row.names = FALSE)

cat("\nD3b DONE.\n")
