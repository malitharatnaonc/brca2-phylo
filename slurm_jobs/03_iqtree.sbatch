#!/bin/bash
#SBATCH --job-name=brca2-iqtree
#SBATCH --partition=ncpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
# Bumped from 32G - a BOOTSTRAP=0 test run reported "ModelFinder requires
# 31412 MB RAM", leaving only ~600MB headroom at the old setting. ncpu's
# per-cpu cap (28000mb x 16 cpus here = 448000mb max) means there's no real
# cost to this much margin.
#SBATCH --time=7-00:00:00
# With MSET=LG,WAG,JTT (44 models instead of 1232), timing data shows ~10.6
# min/model (7 models in 74 minutes) - model selection alone should take
# ~8 hours, before the ML tree search or any of the 1000 bootstrap
# replicates. Genuinely uncertain how long the rest takes, so this requests
# ncpu's full 7-day maximum rather than guessing a tighter number and
# risking a timeout that throws away most of a week's progress. If it DOES
# need to be resubmitted, IQ-TREE checkpoints its own progress periodically
# (by default, into a .ckp.gz file alongside the other -pre outputs) -
# resubmitting the exact same command (same PREFIX) should resume from
# there rather than restarting from scratch.
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
# Restricts ModelFinder to a comma-separated list of substitution matrices
# instead of testing all 1232 model combinations (infeasible at this scale -
# see --time comment above). Empty by default, preserving the old
# unrestricted behaviour for smaller datasets; set e.g. MSET=LG,WAG,JTT for
# a large run. A partial ModelFinder test already showed LG+G4/LG+I+G4
# dramatically outperforming plain LG, which is the standard justification
# for using LG as the base matrix in protein phylogenetics generally - not
# a guarantee this is optimal for BRCA2 specifically, just a reasonable,
# well-precedented default given the evidence so far.
MSET="${MSET:-}"
MSET_FLAG=()
if [[ -n "${MSET}" ]]; then
  MSET_FLAG=(-mset "${MSET}")
fi

# MRATE restricts rate-heterogeneity models tested (e.g. MRATE=G,R drops +I
# and +I+R combinations). CMAX caps FreeRate category count (default: 10).
# Both added after real timing data showed diminishing/redundant returns at
# high R-category counts and near-total overlap between +R and +I+R fits -
# an evidence-based cut, not a guess.
MRATE="${MRATE:-}"
MRATE_FLAG=()
if [[ -n "${MRATE}" ]]; then
  MRATE_FLAG=(-mrate "${MRATE}")
fi

CMAX="${CMAX:-}"
CMAX_FLAG=()
if [[ -n "${CMAX}" ]]; then
  CMAX_FLAG=(-cmax "${CMAX}")
fi

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

if [[ -n "${MSET}" ]]; then
  echo "Restricting ModelFinder to: ${MSET}"
fi

if [[ "${BOOTSTRAP}" -eq 0 ]]; then
  echo "BOOTSTRAP=0 - running without bootstrap (timing test only)"
  "${IQTREE_BIN}" -s "${IN}" -m MFP "${MSET_FLAG[@]}" "${MRATE_FLAG[@]}" "${CMAX_FLAG[@]}" -nt "${SLURM_CPUS_PER_TASK}" -pre "${PREFIX}"
else
  "${IQTREE_BIN}" -s "${IN}" -m MFP "${MSET_FLAG[@]}" "${MRATE_FLAG[@]}" "${CMAX_FLAG[@]}" -bb "${BOOTSTRAP}" -nt "${SLURM_CPUS_PER_TASK}" -pre "${PREFIX}"
fi

echo "Finished: $(date)"
echo "Check efficiency with: seff ${SLURM_JOB_ID}"
