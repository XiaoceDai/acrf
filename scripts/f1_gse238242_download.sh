#!/usr/bin/env bash
# F1 — Download GSE238242 (LAA snRNA + snATAC) supplementary files from NCBI GEO.
# GSE238242 multiome download (RNA-side support / regulatory layer).
#
# Primary files (HTTPS):
#   GSE238242_snAF.metadata.tsv.gz  (~73K)
#   GSE238242_RAW.tar               (~179M)  — per-sample RNA/ATAC count matrices
#
# Usage:
#   bash scripts/f1_gse238242_download.sh
# Optional:
#   SKIP_RAW=1  bash scripts/f1_gse238242_download.sh   # metadata only
#
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${ROOT}/Project_P5/GSE238242"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE238nnn/GSE238242/suppl"
mkdir -p "${OUT}"

fetch() {
  local url="$1" dest="$2"
  if [[ -f "${dest}" ]]; then
    echo "exists: ${dest}"
    return 0
  fi
  echo "downloading → ${dest}"
  curl -fsSL --connect-timeout 30 --retry 3 --retry-delay 2 -o "${dest}.partial" "${url}"
  mv "${dest}.partial" "${dest}"
}

fetch "${BASE}/GSE238242_snAF.metadata.tsv.gz" "${OUT}/GSE238242_snAF.metadata.tsv.gz"
fetch "${BASE}/filelist.txt" "${OUT}/filelist.txt"

if [[ "${SKIP_RAW:-0}" == "1" ]]; then
  echo "SKIP_RAW=1 — not downloading GSE238242_RAW.tar"
  exit 0
fi

fetch "${BASE}/GSE238242_RAW.tar" "${OUT}/GSE238242_RAW.tar"

echo "done. Next: bash scripts/f2_gse238242_extract.sh"
