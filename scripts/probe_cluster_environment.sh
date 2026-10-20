#!/bin/bash
# Check what this pipeline needs is actually available on Nemo, before any
# job submission relies on a guessed module name. Read-only - queries modules
# and storage only, so run it directly on a login node, not as a job:
#
#   bash scripts/probe_cluster_environment.sh 2>&1 | tee cluster_probe.log
#
# Modelled on the colleague repo's probe_mmpbsa_environment.sh pattern:
# establish what's really there before writing a single sbatch script that
# assumes it.

set -uo pipefail

echo "Host: $(hostname)"
echo "Date: $(date)"
echo

echo "=============== 1. home directory usage (5GB limit on Nemo - don't build here)"
du -sh "${HOME}" 2>/dev/null | sed 's/^/  /' || echo "  could not stat \$HOME"
echo

echo "=============== 2. lab space"
echo "  Confirm your lab path - typically /nemo/lab/<lab_name>/home/users/<username>"
echo "  Current working directory: $(pwd)"
df -h . 2>/dev/null | tail -1 | sed 's/^/  /'
echo

echo "=============== 3. sequence-search / phylogenetics tools"
for tool in hmmer mafft iqtree raxml fasttree; do
  echo "  --- ${tool} ---"
  ml av 2>&1 | grep -i "${tool}" | sed 's/^/    /' || echo "    none found via ml av"
done
echo

echo "=============== 4. structure prediction tools"
for tool in esmfold alphafold colabfold singularity apptainer; do
  echo "  --- ${tool} ---"
  ml av 2>&1 | grep -i "${tool}" | sed 's/^/    /' || echo "    none found via ml av"
done
echo

echo "=============== 5. Python"
ml av 2>&1 | grep -oE "Python/[0-9.]+-[A-Za-z0-9.-]+" | sort -u | sed 's/^/  /'
echo

echo "=============== 6. GPU partition sanity check"
sinfo -p ga100 2>&1 | sed 's/^/  /'
echo

echo "=============== Done"
echo
echo "Fill in docs/runbooks/brca2_pipeline_runbook.md 'What's actually on Nemo'"
echo "section with whatever this printed, so the module names in"
echo "slurm_jobs/*.sbatch and workflow/rules/*.smk match reality rather than"
echo "a guess. Anything not found as a module (commonly: HMMER, MAFFT,"
echo "IQ-TREE, ESMFold on a non-bioinformatics-focused cluster) needs either"
echo "a Singularity/Apptainer container, a conda/mamba environment, or a pip"
echo "install into a venv under lab space (see setup_cluster_venv.sh) -"
echo "never under \$HOME, which is only 5GB."
