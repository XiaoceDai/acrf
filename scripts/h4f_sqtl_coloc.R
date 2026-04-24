#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# H4f — Atrial cis-sQTL × AF-GWAS coloc on the 10 anchors + PITX2 calibrator.
#
# Spec lock: FIGURE6_WORKFLOW.md §10.14 (sQTL RNA-layer extension).
#
# Source (locked):
#   sQTL : eQTL Catalogue dataset QTD000255
#          (GTEx v8 heart atrial appendage, leafcutter quantification, n=372)
#   AF GWAS : GCST90624411 (Miyazawa 2023 multi-ancestry, harmonised tabix)
#
# Pipeline per anchor:
#   1. Tabix-slice sQTL within ±1 Mb of the gene.
#   2. Filter rows to those whose gene_id matches the anchor.
#      (Each anchor maps to MANY leafcutter intron clusters /
#       molecular_trait_ids — we coloc each cluster separately.)
#   3. Per-cluster: variant-flip-aware merge with AF GWAS slice → coloc.abf
#      (only run if ≥ 1 SNP has nominal p < 1e-4 in eQTL window — same
#      prefilter convention as h4d).
#   4. Aggregate per-anchor: report PP4_max, PP3_at_max,
#      lead_junction_id, n_clusters_tested.
#
# Output:
#   results/figure6_mr_anchors_v3_sqtl.tsv
#       anchor, family, role, source, tissue,
#       n_clusters_tested, n_clusters_skipped_low_signal,
#       n_clusters_eligible, lead_junction_id,
#       lead_eqtl_p, lead_PP4, lead_PP3,
#       PP4_max, PP3_at_max, n_overlap_at_max
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Rsamtools); library(data.table); library(coloc)
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
ROOT  <- find_root()
CACHE <- file.path(ROOT, "Project_P6", "cache")
OUT_DIR <- file.path(ROOT, "results")
dir.create(CACHE,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GWAS_URL <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000/GCST90624411/harmonised/GCST90624411.h.tsv.gz"
# eQTL Catalogue uses .cc.tsv.gz for sQTL (leafcutter) tabix-indexed sumstats.
SQTL_URL <- "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000255/QTD000255.cc.tsv.gz"

# .cc.tsv.gz schema (verified 2026-04-23 from QTD000255 file head):
#   molecular_trait_id, chromosome, position, ref, alt, variant,
#   ma_samples, maf, pvalue, beta, se, type, ac, an, r2,
#   molecular_trait_object_id, gene_id, median_tpm, rsid
# (same fields as .all.tsv.gz minus the trailing nlog10p column)
SQTL_NAMES <- c("molecular_trait_id","chromosome","position","ref","alt",
                "variant","ma_samples","maf","pvalue","beta","se","type",
                "ac","an","r2","molecular_trait_object_id","gene_id",
                "median_tpm","rsid")
GWAS_NAMES <- c("chromosome","base_pair_location","effect_allele",
                "other_allele","beta","standard_error",
                "effect_allele_frequency","p_value","rsid","rs_id","n",
                "hm_coordinate_conversion","hm_code","variant_id")

# 11 anchors + calibrator — coordinates locked, identical to h4 / h4d.
ANCHORS <- list(
  list(symbol = "PRELP",  family = "F10",        role = "matrix",
       gene_id = "ENSG00000188783", chrom = 1L,
       start = 203475828L, end = 203491352L),
  list(symbol = "OMD",    family = "F10",        role = "matrix",
       gene_id = "ENSG00000127083", chrom = 9L,
       start = 92414245L,  end = 92424461L),
  list(symbol = "IDH2",   family = "Metabolic",  role = "metabolic",
       gene_id = "ENSG00000182054", chrom = 15L,
       start = 90083045L,  end = 90102504L),
  list(symbol = "NDUFA8", family = "F8",         role = "metabolic",
       gene_id = "ENSG00000119421", chrom = 9L,
       start = 122144058L, end = 122159819L),
  list(symbol = "NTRK2",  family = "F15",        role = "receptor",
       gene_id = "ENSG00000148053", chrom = 9L,
       start = 84668522L,  end = 85029699L),
  list(symbol = "NTF3",   family = "F15",        role = "ligand",
       gene_id = "ENSG00000185652", chrom = 12L,
       start = 5404888L,   end = 5495026L),
  list(symbol = "FRS2",   family = "F15",        role = "adaptor",
       gene_id = "ENSG00000166225", chrom = 12L,
       start = 69464533L,  end = 69513650L),
  list(symbol = "BMP10",  family = "F11",        role = "ligand",
       gene_id = "ENSG00000163217", chrom = 2L,
       start = 68869054L,  end = 68886842L),
  list(symbol = "SMAD9",  family = "F11",        role = "R-SMAD",
       gene_id = "ENSG00000120693", chrom = 13L,
       start = 36879326L,  end = 36925087L),
  list(symbol = "SMAD6",  family = "F11",        role = "i-SMAD",
       gene_id = "ENSG00000137834", chrom = 15L,
       start = 66703041L,  end = 66782082L),
  list(symbol = "PITX2",  family = "calibrator", role = "calibrator",
       gene_id = "ENSG00000164093", chrom = 4L,
       start = 110610001L, end = 110633149L)
)

WINDOW_BP        <- 1e6
EQTL_PREFILTER_P <- 1e-4
MIN_OVERLAP      <- 50L

# ---------------------------------------------------------------------------
# Tabix slice with on-disk cache (identical convention to h4 / h4d).
# ---------------------------------------------------------------------------
slice_tabix <- function(url, chrom, start, end, source_tag, retries = 3) {
  key <- sprintf("%s_chr%s_%s_%s.tsv.gz", source_tag, chrom, start, end)
  cache <- file.path(CACHE, key)
  if (file.exists(cache) && file.size(cache) > 0) {
    cat(sprintf("  cached: %s (%.1f MB)\n",
                basename(cache), file.size(cache) / 1024^2))
    return(fread(cache, sep = "\t", header = FALSE, showProgress = FALSE))
  }
  rng <- GRanges(as.character(chrom), IRanges(start = start, end = end))
  tf  <- TabixFile(url, index = paste0(url, ".tbi"))
  cat(sprintf("  tabix slice: %s @ %s:%d-%d\n",
              source_tag, as.character(chrom), start, end))
  for (k in seq_len(retries)) {
    txt <- tryCatch(scanTabix(tf, param = rng)[[1]],
                    error = function(e) {
                      cat(sprintf("    [retry %d] %s\n", k, conditionMessage(e)))
                      NULL
                    })
    if (!is.null(txt)) break
    Sys.sleep(2 * k)
  }
  if (is.null(txt) || length(txt) == 0L) {
    cat("    [empty slice]\n")
    return(data.table())
  }
  con <- gzfile(cache, "wb"); writeLines(txt, con); close(con)
  fread(cache, sep = "\t", header = FALSE, showProgress = FALSE)
}

# ---------------------------------------------------------------------------
# Per-cluster coloc helper
# ---------------------------------------------------------------------------
run_cluster_coloc <- function(e_sub, g, cluster_id) {
  e_sub <- e_sub[!is.na(pvalue) & !is.na(beta) & !is.na(se) & se > 0
                  & maf > 0 & maf < 0.5]
  if (nrow(e_sub) < MIN_OVERLAP)
    return(list(skip = TRUE, reason = sprintf("eqtl rows %d", nrow(e_sub))))

  lead_p <- min(e_sub$pvalue, na.rm = TRUE)
  if (lead_p >= EQTL_PREFILTER_P)
    return(list(skip = TRUE, reason = sprintf("lead p %.2e >= 1e-4", lead_p)))

  e_sub[, key1 := paste(chromosome, position, ref, alt, sep = "_")]
  e_sub[, key2 := paste(chromosome, position, alt, ref, sep = "_")]

  joined_a <- merge(e_sub, g, by.x = "key1", by.y = "key",
                    suffixes = c(".e", ".g"))
  joined_a[, allele_flip := FALSE]
  joined_b <- merge(e_sub, g, by.x = "key2", by.y = "key",
                    suffixes = c(".e", ".g"))
  if (nrow(joined_b)) {
    joined_b[, allele_flip := TRUE]
    joined_b[, beta.g := -beta.g]
  }
  setnames(joined_a, "key1", "key"); setnames(joined_b, "key2", "key")
  joined <- rbindlist(list(joined_a, joined_b),
                      use.names = TRUE, fill = TRUE)
  joined <- unique(joined, by = "key")
  joined <- joined[!is.na(maf) & maf > 0 & maf < 0.5]
  joined[, n_eqtl := an %/% 2L]

  if (nrow(joined) < MIN_OVERLAP)
    return(list(skip = TRUE,
                reason = sprintf("overlap %d < %d", nrow(joined), MIN_OVERLAP)))

  ds_eqtl <- list(snp = joined$key, beta = joined$beta.e,
                  varbeta = joined$se^2, MAF = joined$maf,
                  N = max(joined$n_eqtl, na.rm = TRUE),
                  type = "quant", sdY = 1)
  ds_gwas <- list(snp = joined$key, beta = joined$beta.g,
                  varbeta = joined$standard_error^2, MAF = joined$maf,
                  N = as.integer(median(joined$n, na.rm = TRUE)),
                  s = 0.05, type = "cc")

  res <- tryCatch(
    coloc.abf(dataset1 = ds_eqtl, dataset2 = ds_gwas,
              p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
    error = function(e) list(error = conditionMessage(e))
  )
  if (!is.null(res$error))
    return(list(skip = TRUE,
                reason = paste("coloc error:", res$error)))

  pp <- as.numeric(res$summary[c("PP.H0.abf","PP.H1.abf","PP.H2.abf",
                                  "PP.H3.abf","PP.H4.abf")])
  list(
    skip = FALSE,
    cluster_id = cluster_id,
    lead_p     = lead_p,
    n_overlap  = nrow(joined),
    PP0 = pp[1], PP1 = pp[2], PP2 = pp[3],
    PP3 = pp[4], PP4 = pp[5]
  )
}

# ---------------------------------------------------------------------------
# Per-anchor processing
# ---------------------------------------------------------------------------
process_anchor <- function(a) {
  cat(sprintf("\n=== %s [%s · %s] chr%s:%d-%d ===\n",
              a$symbol, a$family, a$role,
              a$chrom, a$start, a$end))
  win_lo <- max(1L, a$start - WINDOW_BP)
  win_hi <- a$end + WINDOW_BP

  # ---- sQTL slice ----
  e <- slice_tabix(SQTL_URL, a$chrom, win_lo, win_hi,
                   sprintf("sqtl_%s", a$symbol))
  if (nrow(e) == 0L) {
    return(data.table(anchor = a$symbol, family = a$family,
                       role = a$role, source = "sqtl_v3_extension",
                       tissue = "atrial_appendage",
                       n_clusters_tested = 0L,
                       n_clusters_skipped_low_signal = 0L,
                       n_clusters_eligible = 0L,
                       lead_junction_id = NA_character_,
                       lead_eqtl_p = NA_real_,
                       lead_PP4 = NA_real_, lead_PP3 = NA_real_,
                       PP4_max = NA_real_, PP3_at_max = NA_real_,
                       n_overlap_at_max = NA_integer_,
                       reason = "sQTL slice empty"))
  }
  setnames(e, SQTL_NAMES[1:ncol(e)])
  # Filter to clusters whose gene_id maps to this anchor.
  e <- e[grepl(paste0("^", a$gene_id, "(\\.|$)"), gene_id)]
  cat(sprintf("  sQTL rows for %s: %d (across %d clusters)\n",
              a$symbol, nrow(e),
              if (nrow(e)) length(unique(e$molecular_trait_id)) else 0L))
  if (nrow(e) == 0L) {
    return(data.table(anchor = a$symbol, family = a$family,
                       role = a$role, source = "sqtl_v3_extension",
                       tissue = "atrial_appendage",
                       n_clusters_tested = 0L,
                       n_clusters_skipped_low_signal = 0L,
                       n_clusters_eligible = 0L,
                       lead_junction_id = NA_character_,
                       lead_eqtl_p = NA_real_,
                       lead_PP4 = NA_real_, lead_PP3 = NA_real_,
                       PP4_max = NA_real_, PP3_at_max = NA_real_,
                       n_overlap_at_max = NA_integer_,
                       reason = "no sQTL rows for gene_id"))
  }
  e[, position := as.integer(position)]
  e[, pvalue := as.numeric(pvalue)]
  e[, beta := as.numeric(beta)]; e[, se := as.numeric(se)]
  e[, maf := as.numeric(maf)]; e[, an := as.integer(an)]

  # ---- GWAS slice (cached from h4 / h4d) ----
  g <- slice_tabix(GWAS_URL, a$chrom, win_lo, win_hi,
                   sprintf("gwas_%s", a$symbol))
  setnames(g, GWAS_NAMES[1:ncol(g)])
  g[, base_pair_location := as.integer(base_pair_location)]
  g[, beta := as.numeric(beta)]
  g[, standard_error := as.numeric(standard_error)]
  g[, p_value := as.numeric(p_value)]
  g[, n := as.integer(n)]
  g <- g[!is.na(beta) & !is.na(standard_error) & standard_error > 0
         & !is.na(p_value)]
  g[, key := paste(chromosome, base_pair_location, effect_allele,
                    other_allele, sep = "_")]
  cat(sprintf("  GWAS rows: %d   min p: %.2e\n",
              nrow(g),
              if (nrow(g)) min(g$p_value, na.rm = TRUE) else NA_real_))

  # ---- per-cluster coloc ----
  clusters <- unique(e$molecular_trait_id)
  per_cluster <- list()
  n_skip <- 0L; n_elig <- 0L
  for (cl in clusters) {
    e_sub <- e[molecular_trait_id == cl]
    res <- run_cluster_coloc(e_sub, g, cl)
    if (isTRUE(res$skip)) {
      n_skip <- n_skip + 1L
    } else {
      n_elig <- n_elig + 1L
      per_cluster[[cl]] <- res
    }
  }
  cat(sprintf("  clusters tested: %d   eligible (passed prefilter+overlap): %d\n",
              length(clusters), n_elig))

  if (n_elig == 0L) {
    return(data.table(anchor = a$symbol, family = a$family,
                       role = a$role, source = "sqtl_v3_extension",
                       tissue = "atrial_appendage",
                       n_clusters_tested = length(clusters),
                       n_clusters_skipped_low_signal = n_skip,
                       n_clusters_eligible = 0L,
                       lead_junction_id = NA_character_,
                       lead_eqtl_p = NA_real_,
                       lead_PP4 = NA_real_, lead_PP3 = NA_real_,
                       PP4_max = NA_real_, PP3_at_max = NA_real_,
                       n_overlap_at_max = NA_integer_,
                       reason = sprintf("0 / %d clusters passed prefilter",
                                         length(clusters))))
  }

  # Aggregate: pick cluster with strongest lead p ("lead junction") AND
  # cluster with max PP4 ("max-PP4 junction").  These can differ.
  pcdt <- rbindlist(lapply(per_cluster, function(x) {
    data.table(cluster_id = x$cluster_id, lead_p = x$lead_p,
               n_overlap = x$n_overlap,
               PP0 = x$PP0, PP1 = x$PP1, PP2 = x$PP2,
               PP3 = x$PP3, PP4 = x$PP4)
  }))
  i_lead <- which.min(pcdt$lead_p)
  i_max  <- which.max(pcdt$PP4)
  cat(sprintf("  lead-p cluster: %s (p=%.2e, PP4=%.3f)\n",
              pcdt$cluster_id[i_lead], pcdt$lead_p[i_lead],
              pcdt$PP4[i_lead]))
  cat(sprintf("  max-PP4 cluster: %s (PP4=%.3f, PP3=%.3f)\n",
              pcdt$cluster_id[i_max], pcdt$PP4[i_max],
              pcdt$PP3[i_max]))

  data.table(
    anchor = a$symbol, family = a$family, role = a$role,
    source = "sqtl_v3_extension", tissue = "atrial_appendage",
    n_clusters_tested = length(clusters),
    n_clusters_skipped_low_signal = n_skip,
    n_clusters_eligible = n_elig,
    lead_junction_id  = pcdt$cluster_id[i_max],   # the junction we report
    lead_eqtl_p       = pcdt$lead_p[i_max],
    lead_PP4          = pcdt$PP4[i_lead],
    lead_PP3          = pcdt$PP3[i_lead],
    PP4_max           = pcdt$PP4[i_max],
    PP3_at_max        = pcdt$PP3[i_max],
    n_overlap_at_max  = pcdt$n_overlap[i_max],
    reason            = sprintf("max PP4 = %.3f over %d eligible clusters",
                                 pcdt$PP4[i_max], n_elig)
  )
}

# ---------------------------------------------------------------------------
# Run all anchors
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
do_subset <- if (length(args)) args else NULL
results <- list()
for (a in ANCHORS) {
  if (!is.null(do_subset) && !(a$symbol %in% do_subset)) next
  results[[a$symbol]] <- process_anchor(a)
}

out <- rbindlist(results, use.names = TRUE, fill = TRUE)
OUT_TSV <- file.path(OUT_DIR, "figure6_mr_anchors_v3_sqtl.tsv")
fwrite(out, OUT_TSV, sep = "\t")

cat("\n=== Final sQTL summary ===\n")
print(out[, .(anchor, family, n_clusters_tested,
              n_clusters_eligible, lead_eqtl_p,
              PP4_max, PP3_at_max, lead_junction_id)])
cat(sprintf("\nWrote %s\n", OUT_TSV))
cat("H4f DONE.\n")
