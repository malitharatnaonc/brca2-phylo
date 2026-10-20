#!/bin/bash
# Build the Python environment for any pipeline steps that run ON the
# cluster (the BRC-repeat hmmsearch wrapper, and a Python interpreter for
# scripts invoked from inside SLURM jobs). Run on a login node, not a job:
#
#   bash scripts/setup_cluster_venv.sh 2>&1 | tee cluster_venv_setup.log
#
# Most of this pipeline's Python work (parsing, metrics, plotting) runs on
# the laptop - see environment/environment.yml for that. This venv is only
# for the handful of steps that have to run where the cluster-only tools
# (hmmer, mafft, esmfold) live.
#
# Safe to re-run: reuses an existing venv, pip skips what's installed.

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
VENV_DIR="${VENV_DIR:-${PROJECT_DIR}/venvs/brca2phylo}"
PYTHON_MODULE="${PYTHON_MODULE:-}"   # fill in after probe_cluster_environment.sh, e.g. Python/3.11.3-GCCcore-12.3.0

# pip's cache and temp build dir default under $HOME, which is a 5GB quota on
# Nemo and fills fast once biopython/metapredict start compiling. Redirect
# both into lab space before anything else runs.
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-${PROJECT_DIR}/.cache/pip}"
export TMPDIR="${TMPDIR_OVERRIDE:-${PROJECT_DIR}/.tmp}"
mkdir -p "${PIP_CACHE_DIR}" "${TMPDIR}" "$(dirname "${VENV_DIR}")"

echo "Project dir: ${PROJECT_DIR}"
echo "Venv:        ${VENV_DIR}"
echo "Pip cache:   ${PIP_CACHE_DIR}"
echo "TMPDIR:      ${TMPDIR}"
echo

ml purge
if [[ -n "${PYTHON_MODULE}" ]]; then
  echo "Loading ${PYTHON_MODULE}..."
  ml "${PYTHON_MODULE}"
else
  echo "PYTHON_MODULE not set - using whatever 'python3' resolves to by default."
  echo "Run scripts/probe_cluster_environment.sh first and set PYTHON_MODULE"
  echo "to one of the Python modules it lists if this picks up the wrong one."
fi

PY="$(command -v python3)"
echo "python3: ${PY}  ($(${PY} --version 2>&1))"
echo

if [[ -x "${VENV_DIR}/bin/python" ]]; then
  echo "Reusing existing venv."
else
  echo "Creating venv..."
  "${PY}" -m venv "${VENV_DIR}"
fi

VPY="${VENV_DIR}/bin/python"
"${VPY}" -m pip install --upgrade pip --quiet
"${VPY}" -m pip install pandas biopython requests --quiet

echo
echo "Done. Activate with:"
echo "  source ${VENV_DIR}/bin/activate"
echo "or call scripts directly with:"
echo "  ${VPY} scripts/count_brc_repeats.py ..."
