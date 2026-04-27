#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# R1 - Parse GSE255612 (Hill et al. 2024, NComms) cell metadata and emit
# the three derived files Figure 3 requires:
#
#   Project_P3/GSE255612/derived/
#   ├── gse255612_celltype_mapping.tsv     mapping dictionary (frozen)
#   ├── gse255612_obs_4compartment.tsv     per-nucleus obs trimmed to main 4 compartments
#   └── gse255612_qc_meta.tsv              donor-level QC summary
#
# All three files are meant to be the ONLY handoff from metadata to r2
# (r2 never re-reads the raw MetaData.txt.gz).
#
# Frozen decisions (pseudobulk compartments + QC flags):
#   Main 4 compartments: atrialCM / fibroblast / endothelial_or_vascular / immune
#   Pericytes, Vascular_Smooth_Muscle: keep_for_main=FALSE, flagged_for_supp=TRUE
#   Mast_Cells: immune, keep_for_main=TRUE
#   Everything else (Adipocytes / Neuronal / Epicardial_Cells / Megakaryocytes):
#       keep_for_main=FALSE, flagged_for_supp=FALSE
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
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
ROOT    <- find_root()
IN_DIR  <- file.path(ROOT, "Project_P3", "GSE255612")
OUT_DIR <- file.path(IN_DIR, "derived")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("ROOT    = %s\n", ROOT))
cat(sprintf("IN_DIR  = %s\n", IN_DIR))
cat(sprintf("OUT_DIR = %s\n", OUT_DIR))

# --------------------------------------------------------------------
# 1. Frozen cell-type mapping
# --------------------------------------------------------------------

mapping <- tribble(
  ~original_cell_type,        ~compartment4,              ~keep_for_main, ~flagged_for_supp, ~note,
  "Cardiomyocytes",           "atrialCM",                  TRUE,          FALSE,             "Atrial CM",
  "Fibroblasts",              "fibroblast",                TRUE,          FALSE,             "Interstitial FB",
  "Endothelial_1",            "endothelial_or_vascular",   TRUE,          FALSE,             "Vascular endothelial (subset 1)",
  "Endothelial_2",            "endothelial_or_vascular",   TRUE,          FALSE,             "Vascular endothelial (subset 2)",
  "Lymphatic_Endothelial",    "endothelial_or_vascular",   TRUE,          FALSE,             "Lymphatic endothelial; Fig 1 RULES include here",
  "Macrophages",              "immune",                    TRUE,          FALSE,             "Myeloid (dominant immune pop)",
  "Lymphocytes",              "immune",                    TRUE,          FALSE,             "T/B-lineage lymphoid",
  "B_Cells",                  "immune",                    TRUE,          FALSE,             "Lymphoid",
  "Mast_Cells",               "immune",                    TRUE,          FALSE,             "Fig 1 RULES explicitly include mast in immune",
  "Pericytes",                NA_character_,               FALSE,         TRUE,              "Mural cell; excluded from main to match Fig 1 definition; kept for Supp",
  "Vascular_Smooth_Muscle",   NA_character_,               FALSE,         TRUE,              "Mural cell; excluded from main; kept for Supp",
  "Adipocytes",               NA_character_,               FALSE,         FALSE,             "Not in Fig 1 scope",
  "Neuronal",                 NA_character_,               FALSE,         FALSE,             "Not in Fig 1 scope",
  "Epicardial_Cells",         NA_character_,               FALSE,         FALSE,             "Not in Fig 1 scope",
  "Megakaryocytes",           NA_character_,               FALSE,         FALSE,             "n=58, too sparse + not in Fig 1 scope"
)
stopifnot(!anyDuplicated(mapping$original_cell_type))


# --------------------------------------------------------------------
# 2. Load the GSE255612 per-cell MetaData
# --------------------------------------------------------------------
# Broad/SCP format: first data row is "NAME\t..." header we want,
# second data row is "TYPE\tgroup/numeric/...".
meta_path <- file.path(IN_DIR, "GSE255612_AF_snRNA_MetaData.txt.gz")
stopifnot(file.exists(meta_path))

cat("reading MetaData ...\n")
# pick only the columns we need for r2; the full metadata (~20 cols) is not required downstream
obs_full <- read_tsv(meta_path, col_types = cols(.default = col_character()))
# Drop the "TYPE" row (first data row after header).
obs <- obs_full[obs_full$NAME != "TYPE", , drop = FALSE]
cat(sprintf("  raw rows: %d  (after TYPE row drop)\n", nrow(obs)))

# Numeric coercion for the columns we actually use downstream.
num_cols <- c("n_umi", "n_genes", "cellranger_percent_mito",
              "doublet_score", "exon_prop", "entropy")
for (c in num_cols) if (c %in% colnames(obs)) obs[[c]] <- as.numeric(obs[[c]])

# Sanity: one barcode per row.
stopifnot(anyDuplicated(obs$NAME) == 0)


# --------------------------------------------------------------------
# 3. Apply mapping; emit obs_4compartment and celltype_mapping outputs
# --------------------------------------------------------------------
obs <- obs %>%
  rename(cell_id = NAME, cell_type_original = cell_type) %>%
  left_join(mapping, by = c("cell_type_original" = "original_cell_type"))

# Sanity: every cell_type_original covered by the mapping.
uncovered <- obs %>% filter(is.na(compartment4) & !is.na(cell_type_original) &
                            !cell_type_original %in% mapping$original_cell_type) %>%
                     distinct(cell_type_original)
if (nrow(uncovered) > 0) {
  stop("Uncovered cell_type(s) found in MetaData: ",
       paste(uncovered$cell_type_original, collapse = ", "))
}

# Harmonise rhythm column: "case" -> "AF", "control" -> "Ctrl".
obs <- obs %>%
  mutate(rhythm = case_when(
    af == "case"    ~ "AF",
    af == "control" ~ "Ctrl",
    TRUE            ~ NA_character_
  ))
stopifnot(all(obs$rhythm %in% c("AF", "Ctrl")))

# obs_4compartment: keep MAIN + flagged_for_supp rows only; drop excluded
# non-scope cells (Adipocytes / Neuronal / Epicardial / Megakaryocyte) so
# r2 never has to see them.  Downstream main panel filters on keep_for_main.
obs_out <- obs %>%
  filter(keep_for_main | flagged_for_supp) %>%
  select(cell_id, biosample_id, donor_id, rhythm, sex,
         cell_type_original, compartment4,
         keep_for_main, flagged_for_supp,
         n_umi, n_genes, cellranger_percent_mito,
         doublet_score)

out_obs_path <- file.path(OUT_DIR, "gse255612_obs_4compartment.tsv")
write_tsv(obs_out, out_obs_path, na = "")
cat(sprintf("wrote %s  (%d nuclei kept: %d main, %d supp)\n",
            sub(paste0("^", ROOT, "/"), "", out_obs_path),
            nrow(obs_out), sum(obs_out$keep_for_main), sum(obs_out$flagged_for_supp)))

# celltype_mapping output: mapping table + live counts from this dataset
# so r2 / Supp scripts can index without re-reading MetaData.
type_counts <- obs %>%
  count(cell_type_original, name = "n_nuclei_total") %>%
  arrange(desc(n_nuclei_total))
mapping_out <- mapping %>%
  left_join(type_counts, by = c("original_cell_type" = "cell_type_original")) %>%
  mutate(n_nuclei_total = tidyr::replace_na(n_nuclei_total, 0L)) %>%
  arrange(desc(n_nuclei_total))

out_map_path <- file.path(OUT_DIR, "gse255612_celltype_mapping.tsv")
write_tsv(mapping_out, out_map_path, na = "")
cat(sprintf("wrote %s\n", sub(paste0("^", ROOT, "/"), "", out_map_path)))


# --------------------------------------------------------------------
# 4. Donor-level QC summary
# --------------------------------------------------------------------
donors <- obs %>%
  distinct(donor_id, biosample_id, rhythm, sex)

# per donor x compartment4 cell counts (main 4 compartments only + Pericyte/VSMC as supp)
donor_ct_wide <- obs %>%
  filter(keep_for_main | flagged_for_supp) %>%
  mutate(ct_key = ifelse(is.na(compartment4), cell_type_original, compartment4)) %>%
  count(donor_id, ct_key) %>%
  pivot_wider(names_from = ct_key, values_from = n, values_fill = 0L)

# total + percent_mito summary per donor
donor_qc <- obs %>%
  group_by(donor_id) %>%
  summarise(
    n_nuclei_total   = n(),
    n_main           = sum(keep_for_main),
    n_supp           = sum(flagged_for_supp),
    n_excluded       = sum(!(keep_for_main | flagged_for_supp)),
    median_n_umi     = median(n_umi, na.rm = TRUE),
    median_n_genes   = median(n_genes, na.rm = TRUE),
    median_pct_mito  = median(cellranger_percent_mito, na.rm = TRUE),
    median_doublet   = median(doublet_score, na.rm = TRUE),
    .groups = "drop"
  )

donor_meta <- donors %>%
  left_join(donor_qc, by = "donor_id") %>%
  left_join(donor_ct_wide, by = "donor_id") %>%
  arrange(rhythm, donor_id)

out_qc_path <- file.path(OUT_DIR, "gse255612_qc_meta.tsv")
write_tsv(donor_meta, out_qc_path, na = "")
cat(sprintf("wrote %s\n", sub(paste0("^", ROOT, "/"), "", out_qc_path)))


# --------------------------------------------------------------------
# 5. Stdout summary
# --------------------------------------------------------------------
cat("\n===== summary =====\n")
cat(sprintf("Total donors: %d  (AF=%d, Ctrl=%d)\n",
            nrow(donors), sum(donors$rhythm == "AF"), sum(donors$rhythm == "Ctrl")))

# rough main-compartment donor x cell distribution
compartment_stats <- obs %>%
  filter(keep_for_main) %>%
  group_by(compartment4, rhythm) %>%
  summarise(
    n_donors = n_distinct(donor_id),
    n_nuclei = n(),
    median_per_donor = as.integer(median(table(donor_id))),
    min_per_donor = min(table(donor_id)),
    .groups = "drop"
  )
cat("\n--- per-compartment donor coverage (main only) ---\n")
print(as.data.frame(compartment_stats), row.names = FALSE)

# donors at risk of being filtered by §3.6 gate (n_cells < 30 per compartment)
at_risk <- obs %>%
  filter(keep_for_main) %>%
  group_by(donor_id, rhythm, compartment4) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n < 30)
cat(sprintf("\n--- donors with compartment n_cells < 30 (would be dropped by spec §3.6): %d entries ---\n",
            nrow(at_risk)))
if (nrow(at_risk) > 0) {
  print(as.data.frame(at_risk), row.names = FALSE)
}

cat("\nR1 DONE.\n")
