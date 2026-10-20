#!/bin/bash
# Submit MAFFT -> IQ-TREE as a chained pair: IQ-TREE only runs if MAFFT
# exits 0. Checks the input exists and isn't empty before submitting anything,
# so a typo in the filtered-fasta path fails immediately instead of after a
# day in the queue. Run from the project root:
#
#   bash scripts/submit_msa_tree_pipeline.sh

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PROJECT_DIR}"

IN="${IN:-data/processed/homologs_filtered.fasta}"

if [[ ! -s "${IN}" ]]; then
  echo "Filtered homolog fasta not found or empty: ${IN}" >&2
  echo "Run the parse/taxonomy/filter steps on the laptop first (see README)." >&2
  exit 1
fi

N_SEQ=$(grep -c "^>" "${IN}")
echo "Found ${N_SEQ} sequences in ${IN}"

mafft_jobid="$(sbatch --parsable --export=ALL,IN="${IN}" slurm_jobs/02_mafft.sbatch)"
echo "Submitted MAFFT: ${mafft_jobid}"

iqtree_jobid="$(sbatch --parsable --dependency=afterok:"${mafft_jobid}" slurm_jobs/03_iqtree.sbatch)"
echo "Submitted IQ-TREE after MAFFT success: ${iqtree_jobid}"

echo
echo "Monitor with:"
echo "  squeue -u ${USER}"
echo
echo "If MAFFT fails, IQ-TREE will stay PENDING with DependencyNeverSatisfied -"
echo "check logs/mafft/brca2-mafft-${mafft_jobid}.err and resubmit once fixed."
