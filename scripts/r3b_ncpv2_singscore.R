#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# R3b - NCP_v2 program singscore on GSE255612 LA snRNA pseudobulk.
#
# Refreshes Fig 3 data against the NCP_v2 panel (30 programs, 11 axis families)
# defined in results/fig1c_ncp_panel_v2.tsv.
#
# Pipeline:
#   1. Load NCP_v2 30 pathways (program IDs + direction).
#   2. Re-derive leading-edge genes per pathway by running fgsea on the
#      bulk LA-vs-RA DE rank (bulk_chamber_DE_full.tsv) restricted to the
#      30 pathways.  This gives a self-consistent LE per program in the
#      same DE rank that defined NCP_v2.
#   3. Score each (donor, compartment) sample of the GSE255612 LA snRNA
#      pseudobulk with singscore on those LE genes.
#
# Outputs:
#   results/fig3_ncpv2_program_le.tsv             (program_id × LE genes)
#   Project_P3/GSE255612/program_scores_v2.tsv     (donor × ct × program)
#
# Hard rule: NCP_v2 panel is fixed - we never re-select pathways here.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(fgsea)
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
# 1. Load NCP_v2 30 pathways
# --------------------------------------------------------------------

panel <- fread(file.path(ROOT, "results", "fig1c_ncp_panel_v2.tsv"))
cat(sprintf("loaded NCP_v2 panel: %d programs\n", nrow(panel)))

# attach axis-family ids from fig1e
fam <- fread(file.path(ROOT, "results", "fig1e_baseline_axis_families.tsv"))
panel <- merge(panel,
               fam[, .(pathway, family_id)],
               by = "pathway", all.x = TRUE, sort = FALSE)
setorder(panel, program_order)
panel[, baseline_direction := ifelse(direction == "LA-biased", 1L, -1L)]


# --------------------------------------------------------------------
# 2. Build leading-edge genes per NCP_v2 pathway
#    via fgsea on the bulk LA-vs-RA DE rank.
# --------------------------------------------------------------------

de <- fread(file.path(ROOT, "bulk", "Project_P01", "bulk_chamber_DE_full.tsv"))
# rank: signed t-statistic (LA > 0, RA < 0).  Use moderated_t / equiv.
if ("moderated_t" %in% colnames(de)) {
  de[, t_stat := moderated_t]
}
stopifnot("t_stat" %in% colnames(de))
de <- de[!is.na(gene_symbol) & nzchar(gene_symbol) & is.finite(t_stat)]
# collapse duplicate gene symbols by absolute t (keep extreme)
de[, abs_t := abs(t_stat)]
setorder(de, gene_symbol, -abs_t)
de <- de[!duplicated(gene_symbol)]
rank_vec <- de$t_stat
names(rank_vec) <- de$gene_symbol
cat(sprintf("DE rank: %d unique genes\n", length(rank_vec)))


# load gene sets (Hallmark + Reactome) from the same gmt files used originally
gmt_dir  <- file.path(ROOT, "bulk", "Project_P01", "gene_sets")
gmt_h    <- gmtPathways(file.path(gmt_dir, "h.all.v2026.1.Hs.symbols.gmt"))
gmt_r    <- gmtPathways(file.path(gmt_dir, "c2.cp.reactome.v2026.1.Hs.symbols.gmt"))
all_sets <- c(gmt_h, gmt_r)

# the panel's pathway IDs use the prefix "HALLMARK__..." or "REACTOME__..."
# whereas gmt names are just the bare set names (e.g. HALLMARK_OXIDATIVE_PHOSPHORYLATION).
strip_prefix <- function(s) sub("^[A-Z]+__", "", s)
panel[, gmt_key := strip_prefix(pathway)]

missing <- panel[!gmt_key %in% names(all_sets)]
if (nrow(missing)) {
  cat("WARN: pathways not found in GMT (will be skipped):\n")
  print(missing[, .(pathway, gmt_key)])
}
panel <- panel[gmt_key %in% names(all_sets)]

panel_sets <- all_sets[panel$gmt_key]
names(panel_sets) <- panel$pathway

cat(sprintf("running fgsea on %d NCP_v2 pathways...\n", length(panel_sets)))
set.seed(20260418)
fg <- fgsea(pathways = panel_sets, stats = rank_vec, nPermSimple = 50000,
             scoreType = "std")
cat(sprintf("fgsea done: %d pathways scored\n", nrow(fg)))

# extract leading-edge per pathway (signed by panel direction)
fg <- as.data.table(fg)
fg[, le_genes_str := vapply(leadingEdge, paste, collapse = ",",
                              FUN.VALUE = character(1))]
fg[, n_LE_total := vapply(leadingEdge, length, integer(1))]

panel <- merge(panel,
               fg[, .(pathway, NES_v2 = NES, padj_v2 = padj,
                       n_LE_total, le_genes_str)],
               by = "pathway", all.x = TRUE, sort = FALSE)
setorder(panel, program_order)

# Write LE table
le_out <- panel[, .(program_id = pathway,
                       program_order, family_id, direction,
                       baseline_direction,
                       n_LE_total, leading_edge_genes = le_genes_str)]
fwrite(le_out, file.path(ROOT, "results", "fig3_ncpv2_program_le.tsv"),
        sep = "\t")
cat(sprintf("wrote results/fig3_ncpv2_program_le.tsv (%d rows)\n",
            nrow(le_out)))


# --------------------------------------------------------------------
# 3. Score each (donor x compartment) sample on GSE255612 LA snRNA
# --------------------------------------------------------------------

panel[, le_genes := strsplit(le_genes_str, ",", fixed = TRUE)]

# short label map for plotting (reuse Fig 1 / Fig 2 abbreviation table)
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
  fb <- gsub("^[A-Z]+__[A-Z]+_+", "", s)
  fb <- tolower(gsub("_+", " ", fb))
  ifelse(is.na(out), substr(fb, 1, 30), out)
}
panel[, program_label := .short(pathway)]


# Use the 13-compartment pseudobulk_full (matches Fig 1f attribution).
# Compartments without enough donors per AF/Ctrl arm will be filtered
# downstream by R4b (MIN_DONORS = 4); we still score everything here.
PB_BASE <- file.path(ROOT, "Project_P3", "GSE255612", "pseudobulk_full")
gse_groups <- list()
for (g in list.files(PB_BASE)) {
  gdir <- file.path(PB_BASE, g)
  if (file.exists(file.path(gdir, "expr.tsv")) &&
      file.exists(file.path(gdir, "samples.tsv"))) {
    gse_groups[[g]] <- list(supp = FALSE)
  }
}
cat(sprintf("found %d compartments under pseudobulk_full\n",
            length(gse_groups)))

all_scores <- list()
for (g in names(gse_groups)) {
  gdir <- file.path(PB_BASE, g)
  if (!file.exists(file.path(gdir, "expr.tsv"))) {
    cat(sprintf("  SKIP %s (no pseudobulk available)\n", g))
    next
  }
  expr <- read.delim(file.path(gdir, "expr.tsv"),
                     check.names = FALSE, row.names = 1)
  expr <- as.matrix(expr)
  samps <- fread(file.path(gdir, "samples.tsv"))
  expr <- expr[, samps$donor_id, drop = FALSE]
  stopifnot(all(is.finite(expr)))

  ranked <- rankGenes(expr)
  for (i in seq_len(nrow(panel))) {
    le <- panel$le_genes[[i]]
    have <- intersect(le, rownames(expr))
    if (length(have) < 3) next
    res <- simpleScore(ranked, upSet = have, centerScore = TRUE)
    df <- data.table(
      cohort             = "GSE255612",
      compartment        = g,
      is_supp            = gse_groups[[g]]$supp,
      donor_id           = samps$donor_id,
      biosample_id       = if ("biosample_id" %in% names(samps))
                              samps$biosample_id else NA_character_,
      rhythm             = samps$rhythm,
      sex                = if ("sex" %in% names(samps)) samps$sex
                              else NA_character_,
      n_cells            = samps$n_cells,
      program_id         = panel$pathway[i],
      program_label      = panel$program_label[i],
      program_order      = panel$program_order[i],
      family_id          = panel$family_id[i],
      baseline_direction = panel$baseline_direction[i],
      n_LE_total         = panel$n_LE_total[i],
      n_LE_matched       = length(have),
      score              = res$TotalScore
    )
    all_scores[[paste(g, i, sep = ":")]] <- df
  }
  cat(sprintf("  scored GSE255612/%s (%d donors, %d programs)\n",
              g, nrow(samps), nrow(panel)))
}
scores <- rbindlist(all_scores)
fwrite(scores, file.path(ROOT, "Project_P3", "GSE255612",
                          "program_scores_v2.tsv"),
        sep = "\t")
cat(sprintf("wrote Project_P3/GSE255612/program_scores_v2.tsv (%d rows)\n",
            nrow(scores)))

cat("\nR3b DONE.\n")
