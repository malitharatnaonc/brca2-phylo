#!/bin/bash
#SBATCH --job-name=brca2-mafft
#SBATCH --partition=ncpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
# Bumped from the original 4h estimate, which was sized before the real
# dataset turned out to be ~1600+ sequences rather than a few hundred.
# --auto should still handle this reasonably (nowhere near as slow as the
# L-INS-i bug above would have been), but there's no tested estimate at this
# scale yet - watch the first run with `sacct`/`seff` and adjust from there.
#SBATCH --output=logs/mafft/%x-%j.out
#SBATCH --error=logs/mafft/%x-%j.err

set -euo pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
cd "${PROJECT_DIR}"

MAFFT_MODULE="${MAFFT_MODULE:-MAFFT/7.526-GCC-14.2.0-with-extensions}"
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

# This branch was backwards in an earlier version of this script: L-INS-i
# (--localpair --maxiterate 1000) is MAFFT's SLOWEST, most accurate mode -
# practical up to roughly 200 sequences, not beyond. It was being triggered
# for the LARGE case, which at this project's actual scale (~1600+ taxa)
# would likely never finish in any reasonable walltime. --auto lets MAFFT
# itself pick an appropriate algorithm for the input size instead (L-INS-i
# for small inputs, progressively faster methods as N grows), which is the
# right default regardless of N - so just use it unconditionally.
if [[ "${N_SEQ}" -le 200 ]]; then
  mafft --localpair --maxiterate 1000 --thread "${SLURM_CPUS_PER_TASK}" "${IN}" > "${OUT}"
else
  mafft --auto --thread "${SLURM_CPUS_PER_TASK}" "${IN}" > "${OUT}"
fi

echo "Finished: $(date)"
echo "Check efficiency with: seff ${SLURM_JOB_ID}"
