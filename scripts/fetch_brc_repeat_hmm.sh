#!/bin/bash
# Fetch the BRC-repeat Pfam HMM (PF00634) into data/external/BRC_repeat.hmm.
# Run once on a login node (needs internet access):
#
#   bash scripts/fetch_brc_repeat_hmm.sh
#
# Tries the lightweight InterPro single-entry API first; if that doesn't
# return a genuine HMMER profile (verified by checking the file actually
# starts with "HMMER3", not just that the HTTP request succeeded), falls
# back to downloading the full Pfam-A HMM library and extracting just this
# one family - slower, but Pfam's own documented, stable distribution
# channel rather than a single-entry API endpoint whose exact output format
# isn't easy to verify in advance.

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PROJECT_DIR}"

PFAM_ACC="PF00634"
OUT="data/external/BRC_repeat.hmm"
mkdir -p data/external

echo "Trying the InterPro single-entry API..."
curl -sL -o "${OUT}.tmp" "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/${PFAM_ACC}?annotation=hmm" || true

if [[ -s "${OUT}.tmp" ]] && head -1 "${OUT}.tmp" | grep -q "^HMMER3"; then
  mv "${OUT}.tmp" "${OUT}"
  echo "Got a valid HMMER3 profile from the InterPro API -> ${OUT}"
else
  echo "That didn't return a usable HMM file (empty, unreachable, or not"
  echo "HMMER3 format) - falling back to the full Pfam-A library."
  rm -f "${OUT}.tmp"

  FETCH_DIR="${PROJECT_DIR}/.tmp/pfam_fetch"
  mkdir -p "${FETCH_DIR}"

  if [[ ! -f "${FETCH_DIR}/Pfam-A.hmm" ]]; then
    echo "Downloading the full Pfam-A HMM library (one-time, a few hundred MB)..."
    curl -L -o "${FETCH_DIR}/Pfam-A.hmm.gz" \
      "http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    gunzip "${FETCH_DIR}/Pfam-A.hmm.gz"
  fi

  python3 - "${FETCH_DIR}/Pfam-A.hmm" "${PFAM_ACC}" "${OUT}" << 'PYEOF'
import sys
hmm_path, target, out_path = sys.argv[1:4]
with open(hmm_path) as f:
    content = f.read()
match = None
for entry in content.split("//\n"):
    for line in entry.splitlines():
        if line.startswith("ACC") and line.split()[1].split(".")[0] == target:
            match = entry + "//\n"
            break
    if match:
        break
if match is None:
    sys.exit(f"Could not find accession {target} in {hmm_path}")
with open(out_path, "w") as f:
    f.write(match)
print(f"Extracted {target} -> {out_path}")
PYEOF

  echo "Done. The full Pfam-A library is cached at ${FETCH_DIR} - delete it"
  echo "if you don't want to keep it around: rm -rf ${FETCH_DIR}"
fi

echo
echo "First lines of ${OUT}:"
head -4 "${OUT}"
