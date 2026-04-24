#!/usr/bin/env bash
# ----------------------------------------------------------------------
# Background downloader for the 4 LDSC reference archives needed by
# the sLDSC pipeline (FIGURE6_WORKFLOW.md §10.15, planned).
#
# Network reality (alkesgroup.broadinstitute.org observed ~1-10 KB/s):
#   * w_hm3.snplist            ~6 MB    (HM3 SNP list for munge)
#   * weights_hm3_no_MHC.tgz   ~150 MB
#   * baselineLD_v2.2_lds.tgz  ~675 MB
#   * 1000G_Phase3_plinkfiles  ~2 GB    (raw genotypes for custom annot)
#
# Total ~3 GB at 5 KB/s = ~7 days.  At 1 KB/s = 35+ days.
# Whichever — well beyond a single shell session.  This script is
# launched detached (nohup) so it survives the parent shell exit;
# user polls progress via `ls -la ldsc_data/` and reads
# `ldsc_data/_download_status.log`.
#
# Each wget call uses --continue (resume partial downloads) and
# unlimited retries so transient drops don't kill the run.
# ----------------------------------------------------------------------
set -u
ROOT="$(cd "$(dirname "$0")/.."; pwd)"
DEST="${ROOT}/ldsc_data"
LOG="${DEST}/_download_status.log"
mkdir -p "${DEST}"

URLS=(
  "https://alkesgroup.broadinstitute.org/LDSCORE/w_hm3.snplist"
  "https://alkesgroup.broadinstitute.org/LDSCORE/1000G_Phase3_weights_hm3_no_MHC.tgz"
  "https://alkesgroup.broadinstitute.org/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz"
  "https://alkesgroup.broadinstitute.org/LDSCORE/1000G_Phase3_plinkfiles.tgz"
)

echo "=== LDSC download started: $(date -Is)" >> "${LOG}"
echo "PID: $$" >> "${LOG}"

for url in "${URLS[@]}"; do
  fname="$(basename "${url}")"
  dest="${DEST}/${fname}"
  echo "" >> "${LOG}"
  echo "--- $(date -Is)  ${fname} ---" >> "${LOG}"
  # --continue: resume partial; --tries=0 = infinite retries;
  # --waitretry: backoff between retries; --no-verbose: less spam.
  wget \
    --continue \
    --tries=0 \
    --waitretry=30 \
    --timeout=120 \
    --no-verbose \
    --append-output="${LOG}" \
    --output-document="${dest}" \
    "${url}"
  echo "  status: $?  size: $(stat -c '%s' "${dest}" 2>/dev/null || echo NA)" >> "${LOG}"
done

echo "" >> "${LOG}"
echo "=== LDSC download finished: $(date -Is) ===" >> "${LOG}"
echo "Final inventory:" >> "${LOG}"
ls -la "${DEST}" >> "${LOG}"
