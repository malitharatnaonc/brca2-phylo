#!/bin/bash
# No IQ-TREE module on Nemo (confirmed via probe_cluster_environment.sh).
# IQ-TREE ships precompiled static Linux binaries that need no compilation
# and no root access - just download, untar, and point at it directly.
# Run on a login node (needs internet access), once:
#
#   bash scripts/install_iqtree_static.sh

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VERSION="${IQTREE_VERSION:-2.4.0}"
INSTALL_DIR="${PROJECT_DIR}/tools/iqtree-${VERSION}"
URL="https://github.com/iqtree/iqtree2/releases/download/v${VERSION}/iqtree-${VERSION}-Linux-intel.tar.gz"

mkdir -p "${PROJECT_DIR}/tools"
cd "${PROJECT_DIR}/tools"

if [[ -x "${INSTALL_DIR}/bin/iqtree2" ]]; then
  echo "Already installed: ${INSTALL_DIR}/bin/iqtree2"
else
  echo "Downloading IQ-TREE ${VERSION}..."
  curl -L -o "iqtree-${VERSION}.tar.gz" "${URL}"
  tar xzf "iqtree-${VERSION}.tar.gz"
  rm "iqtree-${VERSION}.tar.gz"
fi

"${INSTALL_DIR}/bin/iqtree2" --version

echo
echo "Installed. Set this in slurm_jobs/03_iqtree.sbatch (or export before"
echo "submitting):"
echo "  IQTREE_BIN=${INSTALL_DIR}/bin/iqtree2"
