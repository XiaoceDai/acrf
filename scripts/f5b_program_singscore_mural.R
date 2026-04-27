#!/usr/bin/env Rscript
# F5b — TMM/logCPM + singscore on the GSE238242 mural-supp pseudobulk for the
#       6 frozen support programs.  Mural compartment (PC + SMC) is reported
#       as a sensitivity layer (mural carve-out), not part of the main
#       4-compartment analysis.
#
# Output:
#   Project_P5/GSE238242/results/figure5_program_scores_mural.tsv
#   Project_P5/GSE238242/results/figure5_program_delta_mural.tsv
suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(edgeR); library(singscore)
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
PB <- file.path(ROOT, "Project_P5", "GSE238242", "pseudobulk")
OUT <- file.path(ROOT, "Project_P5", "GSE238242", "results")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Frozen 6-program panel
panel_in <- fread(file.path(ROOT, "results", "figure5_program_sets.tsv"))
panel <- panel_in[set_name %in% c("STRUCTURAL", "METABOLIC")]
prog_in <- read.delim(file.path(ROOT, "normal_chamber_programs.tsv"),
                     sep="\t", stringsAsFactors=FALSE)
for (col in c("program_id","program","direction","leading_edge_genes")) {
  if (col %in% colnames(prog_in)) prog_in[[col]] <- gsub('^"|"$','',prog_in[[col]])
}
idx <- match(panel$program_id, prog_in$program_id)
panel[, baseline_direction := ifelse(prog_in$direction[idx] == "LA-biased", 1L, -1L)]
panel[, le_str := prog_in$leading_edge_genes[idx]]
panel[, le_genes := lapply(le_str,
                          function(s) unique(trimws(strsplit(s,",",fixed=TRUE)[[1]])))]

M <- as.matrix(readRDS(file.path(PB, "gse238242_pseudobulk_counts_supp_mural.rds")))
meta <- fread(file.path(PB, "gse238242_pseudobulk_meta_supp_mural.tsv"))
M <- M[, meta$pb_id]
cat(sprintf("loaded mural pseudobulk: %d genes x %d pb\n", nrow(M), ncol(M)))

# TMM + logCPM (one compartment, so all together)
y <- DGEList(counts = M)
keep <- filterByExpr(y, group = meta$rhythm, min.count = 5)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method = "TMM")
lcpm <- cpm(y, log=TRUE, prior.count=1)
cat(sprintf("genes after filter = %d\n", nrow(lcpm)))

# singscore each program
rk <- rankGenes(lcpm)
score_program <- function(le, pid, plab) {
  hit <- intersect(le, rownames(lcpm))
  if (length(hit) < 3) {
    return(data.table(program_id=pid, program_label=plab, n_LE=length(hit),
                      pb_id=colnames(lcpm), score=NA_real_))
  }
  s <- simpleScore(rk, upSet=hit, knownDirection=TRUE, centerScore=TRUE)
  data.table(program_id=pid, program_label=plab, n_LE=length(hit),
             pb_id=rownames(s), score=s$TotalScore)
}

scores <- list()
for (i in seq_len(nrow(panel))) {
  scores[[i]] <- score_program(panel$le_genes[[i]],
                              panel$program_id[i], panel$program_label[i])
}
sc <- rbindlist(scores)
sc <- merge(sc, meta[, .(pb_id, donor, rhythm, n_cells)],
            by = "pb_id", sort = FALSE)
fwrite(sc, file.path(OUT, "figure5_program_scores_mural.tsv"), sep = "\t")
cat(sprintf("wrote %s (%d rows)\n",
            "figure5_program_scores_mural.tsv", nrow(sc)))

# Δ summary
delta <- sc[, .(
  n_AF = sum(rhythm == "AF"), n_SR = sum(rhythm == "SR"),
  mean_AF = mean(score[rhythm == "AF"], na.rm = TRUE),
  mean_SR = mean(score[rhythm == "SR"], na.rm = TRUE)
), by = .(program_id, program_label)]
delta[, delta_score := mean_AF - mean_SR]
panel_axis <- panel[, .(program_id, axis = set_name)]
delta <- merge(delta, panel_axis, by = "program_id", sort = FALSE)

fwrite(delta, file.path(OUT, "figure5_program_delta_mural.tsv"), sep = "\t")
cat("\n--- Δscore (AF − SR) on mural compartment (PC + SMC) ---\n")
print(delta[, .(axis, program_label,
                d = round(delta_score, 3),
                mean_AF = round(mean_AF, 3), mean_SR = round(mean_SR, 3))])
