#!/bin/bash
#SBATCH --job-name=brca2-jackhmmer
#SBATCH --partition=ncpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/jackhmmer/%x-%j.out
#SBATCH --error=logs/jackhmmer/%x-%j.err

# ncpu allows up to 7 days; 5 jackhmmer iterations against a full
# UniRef90/UniProtKB is unpredictable on first run. Test first with
# JACKHMMER_ITER=2 against a small DB slice on `nint` to estimate wall-clock,
# then if the real run looks like it'll fit under 6h, resubmit with
# --time < 06:00:00 for the "short" reservation's preferential scheduling.

set -euo pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
cd "${PROJECT_DIR}"

HMMER_MODULE="${HMMER_MODULE:-}"      # fill in after scripts/probe_cluster_environment.sh
QUERY="${QUERY:-data/raw/human_BRCA2.fasta}"
DB="${DB:-data/external/uniref90.fasta}"
OUTDIR="${OUTDIR:-results/hits}"
JACKHMMER_ITER="${JACKHMMER_ITER:-5}"
JACKHMMER_INCE="${JACKHMMER_INCE:-1e-3}"

mkdir -p "${OUTDIR}" logs/jackhmmer

echo "Job ID: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Project dir: ${PROJECT_DIR}"
echo "Query: ${QUERY}  DB: ${DB}  iterations: ${JACKHMMER_ITER}"
echo "Started: $(date)"

ml purge
if [[ -n "${HMMER_MODULE}" ]]; then
  ml "${HMMER_MODULE}"
else
  echo "HMMER_MODULE not set - run scripts/probe_cluster_environment.sh and" >&2
  echo "set it (e.g. export HMMER_MODULE=HMMER/3.4-...) before relying on this." >&2
fi

if [[ ! -f "${QUERY}" ]]; then
  echo "Query fasta not found: ${QUERY}" >&2
  exit 1
fi
if [[ ! -f "${DB}" ]]; then
  echo "Database not found: ${DB}" >&2
  echo "Check with your HPC team whether an institute-wide copy of" >&2
  echo "UniRef90/UniProtKB already exists before downloading your own." >&2
  exit 1
fi

jackhmmer -N "${JACKHMMER_ITER}" \
    --incE "${JACKHMMER_INCE}" \
    --cpu "${SLURM_CPUS_PER_TASK}" \
    --domtblout "${OUTDIR}/brca2.domtbl" \
    --tblout "${OUTDIR}/brca2_tbl.txt" \
    -A "${OUTDIR}/brca2_align.sto" \
    "${QUERY}" "${DB}"

echo "Finished: $(date)"
echo "Check efficiency with: seff ${SLURM_JOB_ID}"
