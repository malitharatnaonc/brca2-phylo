#!/bin/bash
#SBATCH --job-name=brca2-iqtree
#SBATCH --partition=ncpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
# Bumped from the original 8h estimate - ModelFinder + 1000 ultrafast
# bootstraps on ~1600+ taxa (vs. the few hundred this was originally sized
# for) is a substantially bigger job. This is still a guess, not a tested
# number: ncpu allows up to 7 days, so there's room to go higher if this
# undershoots. Consider testing with a reduced bootstrap count first
# (BOOTSTRAP=100 sbatch ... or similar) to sanity-check convergence/runtime
# before committing to the full 1000-replicate run.
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
# BOOTSTRAP=0 skips -bb entirely - useful for timing just ModelFinder + the
# ML tree search (the fixed cost, incurred regardless of bootstrap count)
# before committing to the real run. IQ-TREE's ultrafast bootstrap (-bb)
# refuses to run below 1000 replicates at all (hard requirement, not a
# suggestion - confirmed the hard way: "#replicates must be >= 1000"), so
# there's no small-N -bb test available; skipping it entirely is the only
# way to get a partial timing read without doing the full 1000.
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

if [[ "${BOOTSTRAP}" -ne 0 && "${BOOTSTRAP}" -lt 1000 ]]; then
  echo "BOOTSTRAP=${BOOTSTRAP} is invalid - IQ-TREE's -bb requires >= 1000" >&2
  echo "replicates, or use BOOTSTRAP=0 to skip bootstrapping entirely for a" >&2
  echo "timing-only test of ModelFinder + the ML tree search." >&2
  exit 1
fi

if [[ "${BOOTSTRAP}" -eq 0 ]]; then
  echo "BOOTSTRAP=0 - running without bootstrap (timing test only)"
  "${IQTREE_BIN}" -s "${IN}" -m MFP -nt "${SLURM_CPUS_PER_TASK}" -pre "${PREFIX}"
else
  "${IQTREE_BIN}" -s "${IN}" -m MFP -bb "${BOOTSTRAP}" -nt "${SLURM_CPUS_PER_TASK}" -pre "${PREFIX}"
fi

echo "Finished: $(date)"
echo "Check efficiency with: seff ${SLURM_JOB_ID}"
