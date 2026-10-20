#!/bin/bash
#SBATCH --job-name=brca2-mafft
#SBATCH --partition=ncpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/mafft/%x-%j.out
#SBATCH --error=logs/mafft/%x-%j.err

set -euo pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
cd "${PROJECT_DIR}"

MAFFT_MODULE="${MAFFT_MODULE:-}"
IN="${IN:-data/processed/homologs_filtered.fasta}"
OUT="${OUT:-results/msa/homologs_aln.fasta}"

mkdir -p results/msa logs/mafft

echo "Job ID: ${SLURM_JOB_ID}  Host: $(hostname)  Started: $(date)"

ml purge
if [[ -n "${MAFFT_MODULE}" ]]; then
  ml "${MAFFT_MODULE}"
else
  echo "MAFFT_MODULE not set - see scripts/probe_cluster_environment.sh" >&2
fi

if [[ ! -f "${IN}" ]]; then
  echo "Input fasta not found: ${IN}" >&2
  exit 1
fi

N_SEQ=$(grep -c "^>" "${IN}")
echo "Aligning ${N_SEQ} sequences"

if [[ "${N_SEQ}" -le 200 ]]; then
  mafft --auto --thread "${SLURM_CPUS_PER_TASK}" "${IN}" > "${OUT}"
else
  mafft --localpair --maxiterate 1000 --thread "${SLURM_CPUS_PER_TASK}" "${IN}" > "${OUT}"
fi

echo "Finished: $(date)"
echo "Check efficiency with: seff ${SLURM_JOB_ID}"
