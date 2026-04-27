#!/usr/bin/env bash
# G1 — Download GSE261363 (NanoString GeoMx DSP, EAT + atrial appendage) supplementary files.
# GSE261363 spatial download (GeoMx ROI / inset cohort).
#
# Files (HTTPS):
#   GSE261363_All_Data_WTA_1_7_22.xlsx  (~8.7M, normalized WTA matrix)
#   GSE261363_RAW.tar                   (~2.4M, per-ROI .dcc.gz files)
#   GSE261363_series_matrix.txt.gz      (~tens KB, sample/ROI metadata)
#
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${ROOT}/Project_P5/GSE261363"
SUPP="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE261nnn/GSE261363/suppl"
MAT="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE261nnn/GSE261363/matrix"
mkdir -p "${OUT}"

fetch() {
  local url="$1" dest="$2"
  if [[ -f "${dest}" ]]; then echo "exists: ${dest}"; return 0; fi
  echo "downloading → ${dest}"
  curl -fsSL --connect-timeout 30 --retry 3 --retry-delay 2 -o "${dest}.partial" "${url}"
  mv "${dest}.partial" "${dest}"
}

fetch "${SUPP}/filelist.txt" "${OUT}/filelist.txt"
fetch "${SUPP}/GSE261363_All_Data_WTA_1_7_22.xlsx" "${OUT}/GSE261363_All_Data_WTA_1_7_22.xlsx"
fetch "${SUPP}/GSE261363_RAW.tar" "${OUT}/GSE261363_RAW.tar"
fetch "${MAT}/GSE261363_series_matrix.txt.gz" "${OUT}/GSE261363_series_matrix.txt.gz"

echo "done. Next: bash scripts/g2_gse261363_extract.sh"
