# BRCA2 homolog pipeline - Nemo runbook

Start here. This is the working guide for running the homolog search,
phylogenetics, and structural-amenability screen on the Crick Nemo cluster.
See the top-level `README.md` for the overall project rationale and
what-runs-where split; this doc is the step-by-step "how".

## 0. What's actually on Nemo (fill this in before trusting any script)

This pipeline needs HMMER (jackhmmer/hmmsearch), MAFFT, IQ-TREE, and either
ESMFold or a route to batch structure prediction. None of these are
confirmed to exist as Nemo modules - run the probe first:

```bash
bash scripts/probe_cluster_environment.sh 2>&1 | tee cluster_probe.log
```

| Tool    | Module found? | Name/version | Notes |
|---------|---------------|---------------|-------|
| HMMER   | ?             |               |       |
| MAFFT   | ?             |               |       |
| IQ-TREE | ?             |               |       |
| ESMFold | ?             |               |       |
| Singularity/Apptainer | ? |          |       |

Update this table after running the probe, and set the corresponding
`*_MODULE`/`*_SETUP` environment variables (see each `slurm_jobs/*.sbatch`
header) before submitting anything for real. If a tool has no module,
options in order of effort: a conda/mamba environment, a pip install into a
venv under lab space (`scripts/setup_cluster_venv.sh`), or a
Singularity/Apptainer container (ask #hpc-beginners if one already exists for
AlphaFold-family tools - ESMFold specifically may not, since it's not
AlphaFold).

## 1. Log in and place the project correctly

```bash
ssh <username>@login.nemo.thecrick.org
```

Your home directory (`/nemo/homefs/home/<username>`) is a 5GB landing area
only - never put the repo, a venv, or any data there. Clone into your lab's
space instead:

```bash
mkdir -p /nemo/lab/<lab_name>/home/users/<username>
cd /nemo/lab/<lab_name>/home/users/<username>
git clone <your-github-repo-url> brca2-phylo
cd brca2-phylo
```

Everything below assumes you're running from this directory on Nemo (it's
`$PROJECT_DIR` / `$SLURM_SUBMIT_DIR` throughout the job scripts).

## 2. Set up the cluster-side Python environment

```bash
bash scripts/setup_cluster_venv.sh 2>&1 | tee cluster_venv_setup.log
```

This builds a venv under `venvs/brca2phylo` (inside the project, in lab
space - never under `$HOME`) with the handful of packages the cluster-side
scripts need (pandas, biopython, requests). It also redirects `PIP_CACHE_DIR`
and `TMPDIR` into lab space, since pip's defaults land under `$HOME` and will
hit the 5GB quota otherwise.

## 3. Configure the run

Edit `config/config.yaml`:
- `database`: point at an institute-wide UniRef90/UniProtKB copy if one
  exists (ask #hpc or #hpc-beginners before downloading your own multi-GB
  copy) - or your own path under lab space.
- `slurm.account`: your lab/STP Slurm account, if required.
- Leave `slurm.partition`/`gpu_partition` as `ncpu`/`ga100` unless you've
  requested access to `gh100`/`gl40`.

## 4. Stage 1 - jackhmmer search (cluster)

```bash
sbatch slurm_jobs/01_jackhmmer.sbatch
squeue -u "$USER"
```

This is the long pole - test first on `nint` with a 2-iteration run against
a small DB slice to get a wall-clock estimate before committing to a 24h
`ncpu` allocation:

```bash
srun --partition=nint --cpus-per-task=4 --mem=16G --time=00:30:00 --pty /bin/bash
# inside: ml <hmmer module>; jackhmmer -N 2 ...
```

Once it finishes:

```bash
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS
seff <jobid>
```

## 5. Stages 2 - laptop: parse, annotate, filter

Pull the results back and run the laptop-side steps (see README's
what-runs-where table):

```bash
# from the laptop
rsync -avz <username>@login.nemo.thecrick.org:/nemo/lab/<lab_name>/home/users/<username>/brca2-phylo/results/hits/ ./results/hits/

mamba env create -f environment/environment.yml   # first time only
mamba activate brca2-phylo-laptop
snakemake --profile profiles/local --snakefile workflow/Snakefile \
    data/processed/homologs_filtered.fasta
```

Check `data/processed/homologs_filtered.tsv` before continuing - this is the
cheapest point to catch a bad coverage cutoff.

## 6. Stage 3 - MSA + tree (cluster)

```bash
# from the laptop: push the filtered fasta back up
rsync -avz data/processed/homologs_filtered.fasta <username>@login.nemo.thecrick.org:/nemo/lab/<lab_name>/home/users/<username>/brca2-phylo/data/processed/

# on Nemo
bash scripts/submit_msa_tree_pipeline.sh
squeue -u "$USER"
```

This submits MAFFT, then IQ-TREE chained with `--dependency=afterok` - IQ-TREE
only starts if MAFFT exits cleanly.

## 7. Stage 4 - laptop: sequence metrics + AlphaFold DB lookup

```bash
# pull msa + tree back
rsync -avz <username>@login.nemo.thecrick.org:.../results/msa/ ./results/msa/
rsync -avz <username>@login.nemo.thecrick.org:.../results/tree/ ./results/tree/

snakemake --profile profiles/local --snakefile workflow/Snakefile \
    results/metrics/homolog_metrics.tsv data/external/afdb_hits.tsv
```

## 8. Stage 5 - ESMFold for whatever AlphaFold DB didn't cover (cluster)

```bash
rsync -avz data/processed/needs_prediction.fasta <username>@login.nemo.thecrick.org:.../data/processed/

# on Nemo, once ESMFOLD_SETUP is confirmed in step 0:
sbatch --export=ALL,ESMFOLD_SETUP="$ESMFOLD_SETUP" slurm_jobs/04_esmfold.sbatch
```

## 9. Stage 6 - laptop: shortlist, then ColabFold for the top candidates

```bash
rsync -avz <username>@login.nemo.thecrick.org:.../results/structures/esmfold/ ./results/structures/esmfold/

snakemake --profile profiles/local --snakefile workflow/Snakefile \
    data/processed/colabfold_shortlist.fasta

rsync -avz data/processed/colabfold_shortlist.fasta <username>@login.nemo.thecrick.org:.../data/processed/

# on Nemo
sbatch --export=ALL,COLABFOLD_SETUP="$COLABFOLD_SETUP" slurm_jobs/05_colabfold.sbatch
```

## 10. Final merge and plots (laptop)

```bash
rsync -avz <username>@login.nemo.thecrick.org:.../results/structures/colabfold/ ./results/structures/colabfold/

snakemake --profile profiles/local --snakefile workflow/Snakefile \
    results/structures/structure_metrics.tsv

jupyter lab notebooks/analysis.py   # or: jupytext --to notebook && open
```

## Keeping Nemo, the laptop, and GitHub in sync

Edit scripts in exactly one place at a time, commit, then sync the other two
- don't hand-edit the same script independently on Nemo and the laptop:

```bash
git add -A
git commit -m "describe the change"
git push                     # from wherever you made the edit
git pull                     # on the other side, before you next run anything there
```

Generated data (`data/processed/`, `results/*`, `logs/`) is git-ignored by
design - it moves between laptop and Nemo via `rsync`, not git. Only
`data/raw/` (the original query fasta) is tracked.

## Troubleshooting

- **A dependency chain job stays PENDING / `DependencyNeverSatisfied`**: the
  job it depended on failed. Check that job's `.err` log under `logs/<tool>/`,
  fix, and resubmit the wrapper - it's safe to rerun.
- **`Errno 122 Disk quota exceeded`**: something is writing under `$HOME`
  (pip cache, conda cache, or the venv itself). Check `PIP_CACHE_DIR`/
  `TMPDIR` are set to a lab-space path (`setup_cluster_venv.sh` does this
  automatically - if you're pip-installing by hand, set them yourself first).
- **A `*_MODULE`/`*_SETUP` job fails immediately**: you haven't filled in
  section 0 of this runbook yet - the job scripts deliberately exit with an
  error rather than silently guessing a module name.
