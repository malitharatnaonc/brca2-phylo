#!/bin/bash
#SBATCH --job-name=brca2-iqtree
#SBATCH --partition=ncpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/iqtree/%x-%j.out
#SBATCH --error=logs/iqtree/%x-%j.err

set -euo pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
cd "${PROJECT_DIR}"

# No IQ-TREE module on Nemo - use the static binary from
# scripts/install_iqtree_static.sh instead.
IQTREE_BIN="${IQTREE_BIN:-${PROJECT_DIR}/tools/iqtree-2.4.0/bin/iqtree2}"
IN="${IN:-results/msa/homologs_aln.fasta}"
PREFIX="${PREFIX:-results/tree/homologs}"
BOOTSTRAP="${BOOTSTRAP:-1000}"

mkdir -p results/tree logs/iqtree

echo "Job ID: ${SLURM_JOB_ID}  Host: $(hostname)  Started: $(date)"

if [[ ! -x "${IQTREE_BIN}" ]]; then
  echo "IQ-TREE binary not found at ${IQTREE_BIN}" >&2
  echo "Run scripts/install_iqtree_static.sh first." >&2
  exit 1
fi

if [[ ! -f "${IN}" ]]; then
  echo "Alignment not found: ${IN}" >&2
  exit 1
fi

"${IQTREE_BIN}" -s "${IN}" -m MFP -bb "${BOOTSTRAP}" -nt "${SLURM_CPUS_PER_TASK}" -pre "${PREFIX}"

echo "Finished: $(date)"
echo "Check efficiency with: seff ${SLURM_JOB_ID}"
