#!/usr/bin/env bash
# G2 — Extract GSE261363_RAW.tar (per-ROI .dcc.gz files).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIR="${ROOT}/Project_P5/GSE261363"
TAR="${DIR}/GSE261363_RAW.tar"
DEST="${DIR}/dcc"
[[ -f "${TAR}" ]] || { echo "Run scripts/g1_gse261363_download.sh first"; exit 1; }
mkdir -p "${DEST}"
echo "extracting ${TAR} → ${DEST}"
tar -xf "${TAR}" -C "${DEST}"
echo "n_files: $(ls -1 "${DEST}" | wc -l)"
echo "Next: Rscript scripts/g3_gse261363_qc.R"
