#!/bin/bash
#SBATCH --job-name=brca2-colabfold
#SBATCH --partition=ga100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/colabfold/%x-%j.out
#SBATCH --error=logs/colabfold/%x-%j.err

set -euo pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
cd "${PROJECT_DIR}"

COLABFOLD_SETUP="${COLABFOLD_SETUP:-}"   # e.g. 'ml Singularity/3.11.3' - fill in after probing
IN="${IN:-data/processed/colabfold_shortlist.fasta}"
OUT="${OUT:-results/structures/colabfold}"

mkdir -p "${OUT}" logs/colabfold

echo "Job ID: ${SLURM_JOB_ID}  Host: $(hostname)  Started: $(date)"

ml purge
if [[ -n "${COLABFOLD_SETUP}" ]]; then
  eval "${COLABFOLD_SETUP}"
else
  echo "COLABFOLD_SETUP not set - see scripts/probe_cluster_environment.sh." >&2
  exit 1
fi

if [[ ! -s "${IN}" ]]; then
  echo "Shortlist fasta missing or empty: ${IN}" >&2
  echo "Run scripts/build_shortlist.py (after ESMFold) before this stage." >&2
  exit 1
fi

colabfold_batch "${IN}" "${OUT}"

echo "Finished: $(date)"
echo "Check efficiency with: seff ${SLURM_JOB_ID}"
