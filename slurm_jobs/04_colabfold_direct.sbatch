#!/bin/bash
#SBATCH --job-name=brca2-colabfold-direct
#SBATCH --partition=ga100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/colabfold/%x-%j.out
#SBATCH --error=logs/colabfold/%x-%j.err

# DEFAULT structure-prediction path on Nemo: ColabFold is a ready module
# (confirmed via probe_cluster_environment.sh), ESMFold is not. Rather than
# custom-installing torch+esm for a cheap first-pass triage, this runs
# ColabFold directly on everything AlphaFold DB didn't already cover.
# Reasonable as long as that set is a few hundred sequences or fewer, given
# the domain-coverage filtering in config.yaml's `filtering.max_hits`.
# If the set turns out much larger, slurm_jobs/04_esmfold.sbatch +
# 05_colabfold.sbatch (triage then shortlist) is the fallback - see the
# runbook's "structure prediction: which path" note.

set -euo pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
cd "${PROJECT_DIR}"

COLABFOLD_MODULE="${COLABFOLD_MODULE:-ColabFold/1.5.2-foss-2022a-CUDA-11.7.0}"
IN="${IN:-data/processed/needs_prediction.fasta}"
OUT="${OUT:-results/structures/colabfold}"

mkdir -p "${OUT}" logs/colabfold

echo "Job ID: ${SLURM_JOB_ID}  Host: $(hostname)  Started: $(date)"

ml purge
ml "${COLABFOLD_MODULE}"

if [[ ! -s "${IN}" ]]; then
  echo "Input fasta missing or empty: ${IN}" >&2
  echo "(Expected if AlphaFold DB already covered everything - check" >&2
  echo "data/external/afdb_hits.tsv before treating this as an error.)" >&2
  exit 1
fi

N_SEQ=$(grep -c "^>" "${IN}")
echo "Predicting ${N_SEQ} structures with ColabFold"
if [[ "${N_SEQ}" -gt 300 ]]; then
  echo "WARNING: ${N_SEQ} sequences is a lot for direct ColabFold - consider" >&2
  echo "the ESMFold-triage + shortlist path instead (slurm_jobs/04_esmfold.sbatch)" >&2
fi

colabfold_batch "${IN}" "${OUT}"

echo "Finished: $(date)"
echo "Check efficiency with: seff ${SLURM_JOB_ID}"
