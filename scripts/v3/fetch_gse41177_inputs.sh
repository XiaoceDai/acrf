#!/usr/bin/env bash
# Fetch authoritative GSE41177 inputs for Figure 4 duration analysis.
#
# Downloads:
#   GSE41177_series_matrix.txt.gz  (metadata for D1)
#   GSE41177_RAW.tar               (CEL files for D2)
#
# Usage:
#   bash scripts/v3/fetch_gse41177_inputs.sh
# Optional:
#   SKIP_RAW=1  bash scripts/v3/fetch_gse41177_inputs.sh
#
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="${ROOT}/Project_P4/GSE41177"
MATRIX_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41177/matrix"
SUPP_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41177/suppl"
mkdir -p "${OUT}"

fetch() {
  local url="$1" dest="$2"
  if [[ -f "${dest}" ]]; then
    echo "exists: ${dest}"
    return 0
  fi
  echo "downloading → ${dest}"
  curl -fsSL --connect-timeout 60 --retry 3 --retry-delay 5 -o "${dest}.partial" "${url}"
  mv "${dest}.partial" "${dest}"
}

fetch "${MATRIX_BASE}/GSE41177_series_matrix.txt.gz" "${OUT}/GSE41177_series_matrix.txt.gz"

if [[ "${SKIP_RAW:-0}" == "1" ]]; then
  echo "SKIP_RAW=1 — not downloading GSE41177_RAW.tar"
  exit 0
fi

fetch "${SUPP_BASE}/GSE41177_RAW.tar" "${OUT}/GSE41177_RAW.tar"

sha256sum "${OUT}/GSE41177_series_matrix.txt.gz" "${OUT}/GSE41177_RAW.tar" \
  > "${OUT}/fetch_manifest.sha256"
echo "Wrote ${OUT}/fetch_manifest.sha256"
