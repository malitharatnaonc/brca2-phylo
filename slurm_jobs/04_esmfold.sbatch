#!/bin/bash
#SBATCH --job-name=brca2-esmfold
#SBATCH --partition=ga100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/esmfold/%x-%j.out
#SBATCH --error=logs/esmfold/%x-%j.err

# ga100 is open-access (4xA100-80GB/node). gh100/gl40 are faster but need an
# HPC access request first (its-helpdesk@crick.ac.uk) - switch partitions
# above once granted. Don't use `vis`: capped at 1 GPU/user and 16h, meant
# for visualisation sessions, not batch inference.

set -euo pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
cd "${PROJECT_DIR}"

ESMFOLD_SETUP="${ESMFOLD_SETUP:-}"   # e.g. 'ml Singularity/3.11.3' or a conda activate line -
                                       # fill in after scripts/probe_cluster_environment.sh
IN="${IN:-data/processed/needs_prediction.fasta}"
OUT="${OUT:-results/structures/esmfold}"
BATCH_SIZE="${BATCH_SIZE:-8}"

mkdir -p "${OUT}" logs/esmfold

echo "Job ID: ${SLURM_JOB_ID}  Host: $(hostname)  Started: $(date)"

ml purge
if [[ -n "${ESMFOLD_SETUP}" ]]; then
  eval "${ESMFOLD_SETUP}"
else
  echo "ESMFOLD_SETUP not set - see scripts/probe_cluster_environment.sh." >&2
  echo "ESMFold is unlikely to be a plain module; expect to need a" >&2
  echo "Singularity container or a conda/venv environment with torch+esm" >&2
  echo "installed (build it under lab space, same pattern as setup_cluster_venv.sh)." >&2
  exit 1
fi

if [[ ! -s "${IN}" ]]; then
  echo "Input fasta missing or empty: ${IN}" >&2
  echo "(This is expected if AlphaFold DB already covered everything -" >&2
  echo "check data/external/afdb_hits.tsv before treating this as an error.)" >&2
  exit 1
fi

python scripts/run_esmfold_batch.py "${IN}" "${OUT}" --batch-size "${BATCH_SIZE}"

echo "Finished: $(date)"
echo "Check efficiency with: seff ${SLURM_JOB_ID}"
