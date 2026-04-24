#!/usr/bin/env bash
# F2 — Extract GSE238242_RAW.tar into Project_P5/GSE238242/extracted/
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIR="${ROOT}/Project_P5/GSE238242"
TAR="${DIR}/GSE238242_RAW.tar"
DEST="${DIR}/extracted"

if [[ ! -f "${TAR}" ]]; then
  echo "Missing ${TAR}. Run: bash scripts/f1_gse238242_download.sh"
  exit 1
fi

mkdir -p "${DEST}"
echo "extracting ${TAR} → ${DEST}"
tar -xf "${TAR}" -C "${DEST}"
echo "done. Files:"
ls -1 "${DEST}" | head -20
echo " ..."
echo "Next: Rscript scripts/f3_gse238242_qc.R"
