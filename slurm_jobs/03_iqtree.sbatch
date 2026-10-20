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

IQTREE_MODULE="${IQTREE_MODULE:-}"
IN="${IN:-results/msa/homologs_aln.fasta}"
PREFIX="${PREFIX:-results/tree/homologs}"
BOOTSTRAP="${BOOTSTRAP:-1000}"

mkdir -p results/tree logs/iqtree

echo "Job ID: ${SLURM_JOB_ID}  Host: $(hostname)  Started: $(date)"

ml purge
if [[ -n "${IQTREE_MODULE}" ]]; then
  ml "${IQTREE_MODULE}"
else
  echo "IQTREE_MODULE not set - see scripts/probe_cluster_environment.sh" >&2
fi

if [[ ! -f "${IN}" ]]; then
  echo "Alignment not found: ${IN}" >&2
  exit 1
fi

iqtree -s "${IN}" -m MFP -bb "${BOOTSTRAP}" -nt "${SLURM_CPUS_PER_TASK}" -pre "${PREFIX}"

echo "Finished: $(date)"
echo "Check efficiency with: seff ${SLURM_JOB_ID}"
