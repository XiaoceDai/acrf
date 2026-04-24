#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# D1 - GSE41177 sample identity card.
#
# Parses series_matrix (4 characteristics rows: organ, gender, age, af_duration)
# plus sample_title to produce the authoritative per-sample metadata for
# Figure 4.  Title pattern is "{Sinus rhythm|AF}_{LAA|LAPV}_{donor_number}".
#
# Outputs (FIGURE4_WORKFLOW.md §5):
#   Project_P4/GSE41177/sample_metadata.tsv
#
# Columns (locked schema):
#   gsm
#   cohort              = "GSE41177"
#   sample_title        raw title
#   donor_id            N{donor_number}  (e.g. N7, N14, N15)
#   rhythm              AF / SR
#   tissue              LAA / LAPV
#   organ_raw           free text from series_matrix (audit)
#   sex                 male / female / NA
#   age                 integer years
#   duration_raw        raw "XM" / "0M" / ">1M"
#   duration_months     integer months for numeric codings; NA for ">1M" / SR
#   duration_sensitivity_A    integer months; ">1M" coded as 1; SR still NA
#   duration_tier_B     SR / short_AF (<48M) / long_AF (>=48M) / NA (ambiguous)
#   cel_filename        exact filename inside GSE41177_RAW.tar
#   keep_for_primary    bool: rhythm==AF AND duration_months NOT NA AND tissue==LAA
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
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
IN_DIR  <- file.path(ROOT, "Project_P4", "GSE41177")
META_DIR <- file.path(ROOT, "Project_P4", "_meta")
SM_PATH <- file.path(IN_DIR, "GSE41177_series_matrix.txt.gz")
if (!file.exists(SM_PATH)) {
  SM_PATH <- file.path(META_DIR, "GSE41177_series_matrix.txt.gz")
}
stopifnot(file.exists(SM_PATH))
cat(sprintf("series_matrix: %s\n", SM_PATH))

# ---------- read raw lines ----------
con <- gzfile(SM_PATH, "rt")
lines <- readLines(con); close(con)

split_line <- function(prefix) {
  ln <- grep(paste0("^!", prefix, "\\b"), lines, value = TRUE)
  if (length(ln) == 0) return(NULL)
  if (length(ln) > 1 && prefix == "Sample_characteristics_ch1") return(ln)
  if (length(ln) > 1) ln <- ln[1]
  parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
  gsub("^\"|\"$", "", parts[-1])
}

gsms   <- split_line("Sample_geo_accession")
titles <- split_line("Sample_title")
stopifnot(length(gsms) == length(titles))
n <- length(gsms)
cat(sprintf("n_samples: %d\n", n))

char_rows <- split_line("Sample_characteristics_ch1")
# each char_rows entry is a full line; split each into values
parse_char_row <- function(line) {
  parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
  gsub("^\"|\"$", "", parts[-1])
}
char_mat <- do.call(rbind, lapply(char_rows, parse_char_row))
# derive field name from first value of each row (before ":")
field_names <- sub(":.*", "", char_mat[, 1])
field_names <- gsub(" ", "_", field_names)
rownames(char_mat) <- field_names
cat(sprintf("characteristics rows: %s\n",
            paste(field_names, collapse = ", ")))

# strip "key:" prefix on each cell
strip_key <- function(v, key) {
  v <- sub(paste0("^", key, ":\\s*"), "", v)
  v
}
organ_raw    <- strip_key(char_mat["organ", ], "organ")
gender       <- strip_key(char_mat["gender", ], "gender")
age_raw      <- strip_key(char_mat["age", ], "age")
duration_raw <- strip_key(char_mat["af_duration", ], "af duration")

# ---------- derive fields ----------
# title pattern
pat <- regmatches(titles,
  regexec("^(Sinus rhythm|AF)_(LAA|LAPV)_([A-Za-z0-9]+)$", titles))
parsed <- do.call(rbind, lapply(pat, function(x) {
  if (length(x) == 0) stop("title regex failed: ", x)
  data.frame(rhythm_txt = x[2], tissue = x[3], donor_num = x[4],
             stringsAsFactors = FALSE)
}))
parsed$rhythm <- ifelse(parsed$rhythm_txt == "AF", "AF", "SR")
parsed$donor_id <- paste0("N", parsed$donor_num)

age <- suppressWarnings(as.integer(sub("Y$", "", age_raw)))

parse_duration <- function(x) {
  x <- trimws(x)
  if (x == "0M") return(0L)
  # strict "NM" integer months
  m <- regmatches(x, regexpr("^([0-9]+)M$", x))
  if (length(m) == 1 && m == x) {
    return(as.integer(sub("M$", "", x)))
  }
  # ">1M" or other non-numeric -> NA
  NA_integer_
}
duration_months <- vapply(duration_raw, parse_duration, integer(1))

parse_duration_sensA <- function(x) {
  x <- trimws(x)
  if (x == ">1M") return(1L)
  m <- regmatches(x, regexpr("^([0-9]+)M$", x))
  if (length(m) == 1 && m == x) return(as.integer(sub("M$", "", x)))
  NA_integer_
}
duration_sensA <- vapply(duration_raw, parse_duration_sensA, integer(1))

# 3-tier ordinal (sensitivity B)
tier <- rep(NA_character_, n)
tier[parsed$rhythm == "SR"] <- "SR"
tier[parsed$rhythm == "AF" & duration_sensA < 48 & !is.na(duration_sensA)] <- "short_AF"
tier[parsed$rhythm == "AF" & duration_sensA >= 48 & !is.na(duration_sensA)] <- "long_AF"

# ---------- CEL filename lookup ----------
cel_list_path <- file.path(IN_DIR, "GSE41177_RAW.tar")
cel_files <- character(0)
if (file.exists(cel_list_path)) {
  cel_files <- system(paste0("tar -tf ", shQuote(cel_list_path),
                             " | grep -E '\\.CEL\\.gz$'"),
                      intern = TRUE)
}
# filename begins with GSM...; match by GSM prefix
cel_filename <- rep(NA_character_, n)
for (i in seq_along(gsms)) {
  hit <- grep(paste0("^", gsms[i], "_.*\\.CEL\\.gz$"), cel_files, value = TRUE)
  if (length(hit) == 1) cel_filename[i] <- hit
}
cat(sprintf("CEL filename matched: %d / %d\n",
            sum(!is.na(cel_filename)), n))

# ---------- assemble ----------
out <- tibble(
  gsm                     = gsms,
  cohort                  = "GSE41177",
  sample_title            = titles,
  donor_id                = parsed$donor_id,
  rhythm                  = parsed$rhythm,
  tissue                  = parsed$tissue,
  organ_raw               = organ_raw,
  sex                     = gender,
  age                     = age,
  duration_raw            = duration_raw,
  duration_months         = duration_months,
  duration_sensitivity_A  = duration_sensA,
  duration_tier_B         = tier,
  cel_filename            = cel_filename,
)
out$keep_for_primary <- (out$rhythm == "AF"
                         & !is.na(out$duration_months)
                         & out$tissue == "LAA")

out_path <- file.path(IN_DIR, "sample_metadata.tsv")
write_tsv(out, out_path, na = "")
cat(sprintf("wrote %s  (%d rows)\n",
            sub(paste0("^", ROOT, "/"), "", out_path), nrow(out)))

# ---------- stdout summary ----------
cat("\n===== sample count summary =====\n")
print(as.data.frame(
  out %>% count(rhythm, tissue)
), row.names = FALSE)

cat("\n===== AF duration distribution (numeric) =====\n")
af_laa_numeric <- out %>%
  filter(rhythm == "AF", tissue == "LAA", !is.na(duration_months))
cat(sprintf("LAA  AF with numeric duration: %d / %d\n",
            nrow(af_laa_numeric),
            sum(out$rhythm == "AF" & out$tissue == "LAA")))
if (nrow(af_laa_numeric) > 0) {
  cat(sprintf("  range: %d-%d months  median=%d\n",
              min(af_laa_numeric$duration_months),
              max(af_laa_numeric$duration_months),
              as.integer(median(af_laa_numeric$duration_months))))
}

cat("\n===== samples with duration_raw = '>1M' (dropped from primary) =====\n")
gt1m <- out %>% filter(duration_raw == ">1M") %>%
  select(gsm, donor_id, tissue, rhythm, duration_raw, duration_sensitivity_A, duration_tier_B)
print(as.data.frame(gt1m), row.names = FALSE)

cat("\n===== keep_for_primary tally =====\n")
print(as.data.frame(out %>% count(keep_for_primary)), row.names = FALSE)
