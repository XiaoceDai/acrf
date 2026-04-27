#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# H4b — Sensitivity coloc on the same 4 anchors using a blood-tissue
#       eQTL panel (eQTL Catalogue dataset QTD000356, GTEx whole blood,
#       n = 670).  Reported as a sensitivity column in fig 6C; verdicts
#       remain those of the tissue-matched primary (atrial appendage).
#
# Spec note: blood sensitivity prefers eQTLGen-style blood panel when available;
# eQTLGen is not in the eQTL Catalogue tabix infrastructure, so GTEx whole
# blood (also tabix-indexed and accessible) is used as the cross-tissue
# sensitivity proxy.  This substitution is documented in the figure caption.
#
# Pipeline mirrors h4_mr_coloc_anchors.R but writes to
#   results/figure6_mr_anchors_blood.tsv
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
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GWAS_URL <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000/GCST90624411/harmonised/GCST90624411.h.tsv.gz"
EQTL_URL <- "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000356/QTD000356.all.tsv.gz"

EQTL_NAMES <- c("molecular_trait_id","chromosome","position","ref","alt",
                "variant","ma_samples","maf","pvalue","beta","se","type",
                "ac","an","r2","molecular_trait_object_id","gene_id",
                "median_tpm","rsid","nlog10p")
GWAS_NAMES <- c("chromosome","base_pair_location","effect_allele",
                "other_allele","beta","standard_error",
                "effect_allele_frequency","p_value","rsid","rs_id","n",
                "hm_coordinate_conversion","hm_code","variant_id")

ANCHORS <- list(
  list(symbol="PRELP",  axis="matrix",    gene_id="ENSG00000188783",
       chrom=1L,  start=203475828L, end=203491352L),
  list(symbol="OMD",    axis="matrix",    gene_id="ENSG00000127083",
       chrom=9L,  start=92414245L,  end=92424461L),
  list(symbol="IDH2",   axis="metabolic", gene_id="ENSG00000182054",
       chrom=15L, start=90083045L,  end=90102504L),
  list(symbol="NDUFA8", axis="metabolic", gene_id="ENSG00000119421",
       chrom=9L,  start=122144058L, end=122159819L)
)

WINDOW_BP    <- 1e6
EQTL_GATE_P  <- 5e-8
PP4_HI <- 0.7; PP4_LO <- 0.3

slice_tabix <- function(url, chrom, start, end, source_tag, retries=3) {
  key <- sprintf("%s_chr%s_%s_%s.tsv.gz", source_tag, chrom, start, end)
  cache <- file.path(CACHE, key)
  if (file.exists(cache) && file.size(cache) > 0) {
    cat(sprintf("  cached: %s\n", basename(cache)))
    return(fread(cache, sep="\t", header=FALSE, showProgress=FALSE))
  }
  rng <- GRanges(as.character(chrom), IRanges(start=start, end=end))
  tf <- TabixFile(url, index=paste0(url, ".tbi"))
  cat(sprintf("  tabix slice: %s @ %s:%d-%d\n",
              source_tag, as.character(chrom), start, end))
  for (k in seq_len(retries)) {
    txt <- tryCatch(scanTabix(tf, param=rng)[[1]],
                    error=function(e) { cat("    [retry]", conditionMessage(e), "\n"); NULL })
    if (!is.null(txt)) break
    Sys.sleep(2*k)
  }
  if (is.null(txt) || length(txt)==0L) return(data.table())
  con <- gzfile(cache, "wb"); writeLines(txt, con); close(con)
  fread(cache, sep="\t", header=FALSE, showProgress=FALSE)
}

process_anchor <- function(a) {
  cat(sprintf("\n=== %s (%s axis, BLOOD sensitivity) chr%s:%d-%d ===\n",
              a$symbol, a$axis, a$chrom, a$start, a$end))
  win_lo <- max(1L, a$start - WINDOW_BP); win_hi <- a$end + WINDOW_BP

  e <- slice_tabix(EQTL_URL, a$chrom, win_lo, win_hi,
                   sprintf("eqtl_blood_%s", a$symbol))
  if (nrow(e) == 0L) {
    return(data.table(anchor=a$symbol, axis=a$axis, n_eqtl=0L,
                      lead_eqtl_p=NA_real_, gate_passed=FALSE,
                      coloc_PP4=NA_real_, blood_verdict="not supported",
                      reason="empty"))
  }
  setnames(e, EQTL_NAMES[1:ncol(e)])
  e <- e[grepl(paste0("^", a$gene_id, "(\\.|$)"), gene_id) |
         grepl(paste0("^", a$gene_id, "(\\.|$)"), molecular_trait_id)]
  cat(sprintf("  eQTL rows for %s in BLOOD: %d\n", a$symbol, nrow(e)))
  if (nrow(e) == 0L) {
    return(data.table(anchor=a$symbol, axis=a$axis, n_eqtl=0L,
                      lead_eqtl_p=NA_real_, gate_passed=FALSE,
                      coloc_PP4=NA_real_, blood_verdict="not supported",
                      reason="no rows for gene"))
  }
  e[, position:=as.integer(position)]; e[, pvalue:=as.numeric(pvalue)]
  e[, beta:=as.numeric(beta)]; e[, se:=as.numeric(se)]
  e[, maf:=as.numeric(maf)]; e[, an:=as.integer(an)]
  e <- e[!is.na(pvalue) & !is.na(beta) & !is.na(se) & se>0 & maf>0 & maf<0.5]
  lead_eqtl_p <- min(e$pvalue, na.rm=TRUE)
  gate_passed <- isTRUE(lead_eqtl_p < EQTL_GATE_P)
  cat(sprintf("  BLOOD eQTL lead p: %.2e   gate (<5e-8): %s\n",
              lead_eqtl_p, gate_passed))

  g <- slice_tabix(GWAS_URL, a$chrom, win_lo, win_hi,
                   sprintf("gwas_%s", a$symbol))
  setnames(g, GWAS_NAMES[1:ncol(g)])
  g[, base_pair_location:=as.integer(base_pair_location)]
  g[, beta:=as.numeric(beta)]; g[, standard_error:=as.numeric(standard_error)]
  g[, p_value:=as.numeric(p_value)]; g[, n:=as.integer(n)]
  g <- g[!is.na(beta) & !is.na(standard_error) & standard_error>0]

  e[, key1:=paste(chromosome, position, ref, alt, sep="_")]
  e[, key2:=paste(chromosome, position, alt, ref, sep="_")]
  g[, key:=paste(chromosome, base_pair_location, effect_allele, other_allele, sep="_")]
  ja <- merge(e, g, by.x="key1", by.y="key", suffixes=c(".e",".g"))
  ja[, allele_flip:=FALSE]
  jb <- merge(e, g, by.x="key2", by.y="key", suffixes=c(".e",".g"))
  if (nrow(jb)) { jb[, allele_flip:=TRUE]; jb[, beta.g:=-beta.g] }
  setnames(ja, "key1", "key"); setnames(jb, "key2", "key")
  joined <- unique(rbindlist(list(ja, jb), use.names=TRUE, fill=TRUE), by="key")
  joined <- joined[!is.na(maf) & maf>0 & maf<0.5]
  joined[, n_eqtl:=an %/% 2L]
  cat(sprintf("  overlap variants: %d\n", nrow(joined)))

  if (nrow(joined) < 50L) {
    return(data.table(anchor=a$symbol, axis=a$axis, n_eqtl=nrow(e),
                      lead_eqtl_p=lead_eqtl_p, gate_passed=gate_passed,
                      coloc_PP4=NA_real_,
                      blood_verdict=if (gate_passed) "inconclusive" else "not supported",
                      reason="overlap<50"))
  }
  ds_e <- list(snp=joined$key, beta=joined$beta.e, varbeta=joined$se^2,
               MAF=joined$maf, N=max(joined$n_eqtl), type="quant", sdY=1)
  ds_g <- list(snp=joined$key, beta=joined$beta.g, varbeta=joined$standard_error^2,
               MAF=joined$maf, N=as.integer(median(joined$n)), s=0.05, type="cc")
  res <- tryCatch(coloc.abf(ds_e, ds_g, p1=1e-4, p2=1e-4, p12=1e-5),
                  error=function(e) list(error=conditionMessage(e)))
  if (!is.null(res$error)) {
    return(data.table(anchor=a$symbol, axis=a$axis, n_eqtl=nrow(e),
                      lead_eqtl_p=lead_eqtl_p, gate_passed=gate_passed,
                      coloc_PP4=NA_real_,
                      blood_verdict=if (gate_passed) "inconclusive" else "not supported",
                      reason=paste("coloc error:", res$error)))
  }
  pp4 <- as.numeric(res$summary["PP.H4.abf"])
  cat(sprintf("  BLOOD coloc PP4 = %.3f\n", pp4))
  v <- if (!gate_passed) "not supported" else
       if (pp4 >= PP4_HI) "supported" else
       if (pp4 >= PP4_LO) "inconclusive" else "not supported"
  data.table(anchor=a$symbol, axis=a$axis, n_eqtl=nrow(e),
             lead_eqtl_p=lead_eqtl_p, gate_passed=gate_passed,
             coloc_PP4=pp4, blood_verdict=v,
             reason=sprintf("gate=%s PP4=%.2f", gate_passed, pp4))
}

results <- list()
for (a in ANCHORS) results[[a$symbol]] <- process_anchor(a)
out <- rbindlist(results, use.names=TRUE, fill=TRUE)

OUT_TSV <- file.path(OUT_DIR, "figure6_mr_anchors_blood.tsv")
fwrite(out, OUT_TSV, sep="\t")
cat("\n=== BLOOD sensitivity summary ===\n"); print(out)
cat("\nH4b DONE.\n")
