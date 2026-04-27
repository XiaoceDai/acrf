#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# H1 — Donor-level axis composite from aCAI_rel.
#
# Reads `Project_P1/program_scores_long.tsv` (per donor x program aCAI_rel
# already aligned + normalised against Fig 1 healthy bulk baseline) and
# computes:
#
#   STRUCTURAL_axis_d = mean over S of aCAI_rel[d,p]
#   METABOLIC_axis_d  = mean over M of aCAI_rel[d,p]
#
# Axis composite: equal weight, no cohort z-scoring.
# Donor-level NA if fewer than 2/3 of the axis programs are present.
#
# Outputs:
#   results/figure6_axis_coverage.tsv        rows = program, cols = cohort
#   results/figure6_axis_composite_long.tsv  donor x axis composite (long)
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
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
SETS <- file.path(ROOT, "results", "figure5_program_sets.tsv")
OUT <- file.path(ROOT, "results")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

MIN_FRAC <- 2 / 3

cat(sprintf("ROOT = %s\n", ROOT))

long <- fread(LONG, showProgress = FALSE)
sets <- fread(SETS, showProgress = FALSE)
struct_ids <- sets[set_name == "STRUCTURAL", program_id]
metab_ids  <- sets[set_name == "METABOLIC",  program_id]

cat(sprintf("STRUCTURAL_SUPPORT_IDS: %d\n", length(struct_ids)))
cat(sprintf("METABOLIC_SUPPORT_IDS : %d\n", length(metab_ids)))

# ---------------------------------------------------------------------
# Coverage table — for each (cohort, program), # of donors with aCAI_rel
# ---------------------------------------------------------------------
support_long <- long[program_id %in% c(struct_ids, metab_ids)]
cov_dt <- support_long[, .(n_donors = sum(!is.na(aCAI_rel))),
                       by = .(cohort, program_id, program_label)]
cov_wide <- dcast(cov_dt, program_id + program_label ~ cohort,
                  value.var = "n_donors", fill = 0L)
fwrite(cov_wide, file.path(OUT, "figure6_axis_coverage.tsv"), sep = "\t")
cat("\n--- Coverage (n_donors with aCAI_rel) ---\n")
print(cov_wide)

# ---------------------------------------------------------------------
# Donor-level composites (per cohort)
# ---------------------------------------------------------------------
compute_axis <- function(sub, axis_ids, axis_name) {
  sub_a <- sub[program_id %in% axis_ids]
  if (nrow(sub_a) == 0) return(NULL)
  sub_a <- sub_a[!is.na(aCAI_rel)]
  by_donor <- sub_a[, .(n_present = .N,
                        composite = mean(aCAI_rel, na.rm = TRUE)),
                    by = .(cohort, donor_id, rhythm)]
  needed_min <- ceiling(length(axis_ids) * MIN_FRAC)
  by_donor[, axis := axis_name]
  by_donor[, total_in_axis := length(axis_ids)]
  by_donor[, eligible := n_present >= needed_min]
  by_donor[eligible == FALSE, composite := NA_real_]
  by_donor
}

struct_comp <- compute_axis(support_long, struct_ids, "STRUCTURAL")
metab_comp  <- compute_axis(support_long, metab_ids,  "METABOLIC")

comp <- rbind(struct_comp, metab_comp)
setorder(comp, axis, cohort, rhythm, donor_id)
fwrite(comp, file.path(OUT, "figure6_axis_composite_long.tsv"), sep = "\t")

cat("\n--- Composite per (axis, cohort, rhythm) ---\n")
summ <- comp[!is.na(composite),
             .(n = .N,
               n_AF = sum(rhythm == "AF"),
               n_SR = sum(rhythm == "SR"),
               mean_AF = mean(composite[rhythm == "AF"], na.rm = TRUE),
               sd_AF   = sd(composite[rhythm == "AF"],   na.rm = TRUE),
               mean_SR = mean(composite[rhythm == "SR"], na.rm = TRUE),
               sd_SR   = sd(composite[rhythm == "SR"],   na.rm = TRUE)),
             by = .(axis, cohort)]
setorder(summ, axis, cohort)
print(summ)

# Drop notice if any (cohort, axis) cell was excluded by 2/3 gate
drop <- comp[eligible == FALSE, .N, by = .(axis, cohort)]
if (nrow(drop)) {
  cat("\n--- Dropped donors by 2/3 coverage gate ---\n")
  print(drop)
}

cat("\nH1 DONE.\n")
