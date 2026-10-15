#!/usr/bin/env Rscript
# Recheck coloc.abf for the 12 atrial anchors from the existing cache.
# This is a new run in the environment recorded by sessionInfo().
# It does not download summary statistics and does not overwrite results/.

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

ROOT <- "/home/bio/Project_P02"
CACHE <- file.path(ROOT, "Project_P6", "cache")
OUT <- file.path(ROOT, "revision/acrf_provenance_check_20261015/outputs")
LOG <- file.path(ROOT, "revision/acrf_provenance_check_20261015/logs/coloc_recheck.log")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Declared before the comparison, not adjusted after seeing deltas.
PP_ABS_TOL <- 1e-6
N_EXACT <- TRUE

EQTL_NAMES <- c("molecular_trait_id", "chromosome", "position", "ref", "alt",
                "variant", "ma_samples", "maf", "pvalue", "beta", "se",
                "type", "ac", "an", "r2", "molecular_trait_object_id",
                "gene_id", "median_tpm", "rsid", "nlog10p")
GWAS_NAMES <- c("chromosome", "base_pair_location", "effect_allele",
                "other_allele", "beta", "standard_error",
                "effect_allele_frequency", "p_value", "rsid", "rs_id",
                "n", "hm_coordinate_conversion", "hm_code", "variant_id")

anchors <- list(
  list(symbol="PRELP",  gene_id="ENSG00000188783", chrom=1L,  start=203475828L, end=203491352L),
  list(symbol="OMD",    gene_id="ENSG00000127083", chrom=9L,  start=92414245L,  end=92424461L),
  list(symbol="IDH2",   gene_id="ENSG00000182054", chrom=15L, start=90083045L,  end=90102504L),
  list(symbol="NDUFA8", gene_id="ENSG00000119421", chrom=9L,  start=122144058L, end=122159819L),
  list(symbol="NTRK2",  gene_id="ENSG00000148053", chrom=9L,  start=84668522L,  end=85029699L),
  list(symbol="NTF3",   gene_id="ENSG00000185652", chrom=12L, start=5404888L,   end=5495026L),
  list(symbol="FRS2",   gene_id="ENSG00000166225", chrom=12L, start=69464533L,  end=69513650L),
  list(symbol="BMP10",  gene_id="ENSG00000163217", chrom=2L,  start=68869054L,  end=68886842L),
  list(symbol="SMAD9",  gene_id="ENSG00000120693", chrom=13L, start=36879326L,  end=36925087L),
  list(symbol="SMAD6",  gene_id="ENSG00000137834", chrom=15L, start=66703041L,  end=66782082L),
  list(symbol="PITX2",  gene_id="ENSG00000164093", chrom=4L,  start=110610001L, end=110633149L),
  list(symbol="KCNN3",  gene_id="ENSG00000143603", chrom=1L,  start=154837802L, end=154961963L)
)

log_con <- file(LOG, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
cat("command: Rscript revision/acrf_provenance_check_20261015/coloc_recheck.R\n")
cat("started: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
cat("coloc version: ", as.character(packageVersion("coloc")), "\n", sep = "")
cat("coloc path: ", find.package("coloc"), "\n", sep = "")
cat("PP absolute tolerance predeclared: ", PP_ABS_TOL, "\n", sep = "")
cat("priors p1=p2=1e-4 p12=1e-5; eqtl type=quant sdY=1; gwas type=cc s=0.05\n")
cat("inputs: existing Project_P6/cache only; no download\n\n")

one <- function(a) {
  win_lo <- max(1L, a$start - 1000000L)
  win_hi <- a$end + 1000000L
  eqtl_path <- file.path(CACHE, sprintf("eqtl_%s_chr%s_%s_%s.tsv.gz", a$symbol, a$chrom, win_lo, win_hi))
  gwas_path <- file.path(CACHE, sprintf("gwas_%s_chr%s_%s_%s.tsv.gz", a$symbol, a$chrom, win_lo, win_hi))
  base <- data.table(anchor = a$symbol, eqtl_cache = eqtl_path, gwas_cache = gwas_path)
  if (!file.exists(eqtl_path) || !file.exists(gwas_path)) {
    base[, `:=`(status = "input_missing", warning = "cache file missing")]
    return(base)
  }
  e <- fread(eqtl_path, sep = "\t", header = FALSE, showProgress = FALSE)
  g <- fread(gwas_path, sep = "\t", header = FALSE, showProgress = FALSE)
  setnames(e, EQTL_NAMES[1:ncol(e)])
  setnames(g, GWAS_NAMES[1:ncol(g)])
  e <- e[grepl(paste0("^", a$gene_id, "(\\.|$)"), gene_id) |
         grepl(paste0("^", a$gene_id, "(\\.|$)"), molecular_trait_id)]
  e[, position := as.integer(position)]
  e[, pvalue := as.numeric(pvalue)]
  e[, beta := as.numeric(beta)]
  e[, se := as.numeric(se)]
  e[, maf := as.numeric(maf)]
  e[, an := as.integer(an)]
  e <- e[!is.na(pvalue) & !is.na(beta) & !is.na(se) & se > 0 & maf > 0 & maf < 0.5]
  g[, base_pair_location := as.integer(base_pair_location)]
  g[, beta := as.numeric(beta)]
  g[, standard_error := as.numeric(standard_error)]
  g[, p_value := as.numeric(p_value)]
  g[, n := as.integer(n)]
  g <- g[!is.na(beta) & !is.na(standard_error) & standard_error > 0 & !is.na(p_value)]
  e[, key1 := paste(chromosome, position, ref, alt, sep = "_")]
  e[, key2 := paste(chromosome, position, alt, ref, sep = "_")]
  g[, key := paste(chromosome, base_pair_location, effect_allele, other_allele, sep = "_")]
  joined_a <- merge(e, g, by.x = "key1", by.y = "key", suffixes = c(".e", ".g"))
  joined_a[, allele_flip := FALSE]
  joined_b <- merge(e, g, by.x = "key2", by.y = "key", suffixes = c(".e", ".g"))
  if (nrow(joined_b)) {
    joined_b[, allele_flip := TRUE]
    joined_b[, beta.g := -beta.g]
  }
  setnames(joined_a, "key1", "key")
  setnames(joined_b, "key2", "key")
  joined <- rbindlist(list(joined_a, joined_b), use.names = TRUE, fill = TRUE)
  joined <- unique(joined, by = "key")
  joined <- joined[!is.na(maf) & maf > 0 & maf < 0.5]
  joined[, sample_n := an %/% 2L]
  warn <- character()
  if (nrow(joined) < 50L) {
    base[, `:=`(
      status = "below_overlap_minimum",
      n_eqtl = nrow(e), n_gwas = nrow(g), n_overlap = nrow(joined),
      lead_eqtl_p = min(e$pvalue), gwas_min_p = if (nrow(g)) min(g$p_value) else NA_real_,
      warning = "overlap < 50; coloc not called; this is not a no-colocalization result"
    )]
    return(base)
  }
  ds_eqtl <- list(snp = joined$key, beta = joined$beta.e, varbeta = joined$se^2,
                  MAF = joined$maf, N = max(joined$sample_n, na.rm = TRUE),
                  type = "quant", sdY = 1)
  ds_gwas <- list(snp = joined$key, beta = joined$beta.g, varbeta = joined$standard_error^2,
                  MAF = joined$maf, N = as.integer(median(joined$n, na.rm = TRUE)),
                  s = 0.05, type = "cc")
  res <- tryCatch(
    withCallingHandlers(
      coloc.abf(dataset1 = ds_eqtl, dataset2 = ds_gwas, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
      warning = function(w) {
        warn <<- c(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(err) list(error = conditionMessage(err))
  )
  if (!is.null(res$error)) {
    base[, `:=`(
      status = "coloc_error", n_eqtl = nrow(e), n_gwas = nrow(g), n_overlap = nrow(joined),
      warning = res$error
    )]
    return(base)
  }
  pp <- as.numeric(res$summary[c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")])
  base[, `:=`(
    status = "rechecked",
    n_eqtl = nrow(e), n_gwas = nrow(g), n_overlap = nrow(joined),
    n_allele_flip = sum(joined$allele_flip),
    eqtl_N = max(joined$sample_n, na.rm = TRUE),
    gwas_N = as.integer(median(joined$n, na.rm = TRUE)),
    lead_eqtl_p = min(e$pvalue), gwas_min_p = min(g$p_value),
    coloc_PP0 = pp[1], coloc_PP1 = pp[2], coloc_PP2 = pp[3],
    coloc_PP3 = pp[4], coloc_PP4 = pp[5],
    PP4_ge_0.70 = pp[5] >= 0.70,
    warning = paste(unique(warn), collapse = " | ")
  )]
  cat(sprintf("%s overlap=%d PP4=%.6g PP3=%.6g\n", a$symbol, nrow(joined), pp[5], pp[4]))
  base
}

out <- rbindlist(lapply(anchors, one), fill = TRUE)
fwrite(out, file.path(OUT, "coloc_recheck_raw.tsv"), sep = "\t")
cat("exit_status: completed\n")
sink(type = "message")
sink(type = "output")
close(log_con)

si <- file(file.path(OUT, "coloc_sessionInfo.txt"), open = "wt")
sink(si)
cat("coloc_version\t", as.character(packageVersion("coloc")), "\n", sep = "")
cat("coloc_path\t", find.package("coloc"), "\n", sep = "")
cat("checked_at\t", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
cat("historical_version\tNOT_RECOVERABLE\n")
cat("this_run\tnew recheck of cached inputs; not a recovered historical session\n\n")
print(sessionInfo())
sink()
close(si)
cat("wrote coloc_recheck_raw.tsv\n")
