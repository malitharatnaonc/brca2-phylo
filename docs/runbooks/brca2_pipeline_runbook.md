# BRCA2 homolog pipeline - Nemo runbook

Start here. This is the working guide for running the homolog search,
phylogenetics, and structural-amenability screen on the Crick Nemo cluster.
See the top-level `README.md` for the overall project rationale and
what-runs-where split; this doc is the step-by-step "how".

## 0. What's actually on Nemo

Confirmed via `scripts/probe_cluster_environment.sh` on 21 Oct 2026 (lab:
`zhangx`, login node `login006`):

| Tool    | Module found? | Name/version used by this pipeline | Notes |
|---------|---------------|--------------------------------------|-------|
| HMMER   | yes           | `HMMER/3.4-gompi-2023a`               | latest `(D)` default |
| MAFFT   | yes           | `MAFFT/7.526-GCC-14.2.0-with-extensions` | latest `(D)` default |
| IQ-TREE | **no**        | static binary instead                 | see below |
| ESMFold | **no**        | not used - see below                  | ColabFold used directly instead |
| ColabFold | yes         | `ColabFold/1.5.2-foss-2022a-CUDA-11.7.0` | GPU, needs `ga100` |
| AlphaFold | yes         | several versions, incl. a ColabFold-bundled one | not used directly - API/ColabFold cover our needs |
| Singularity/Apptainer | yes | `Singularity/3.11.3` / `Apptainer/1.4.2` | available if ever needed |
| Python  | yes           | `Python/3.11.3-GCCcore-12.3.0`        | used by `setup_cluster_venv.sh` |

These are already set as the defaults in `config/config.yaml` and in each
`slurm_jobs/*.sbatch` header - override with `*_MODULE` env vars at submit
time only if you rerun the probe and something's changed (module systems do
get upgraded).

**Two things this changes from the original plan:**

1. **No IQ-TREE module.** It ships precompiled static Linux binaries that
   need no compilation - run this once:
   ```bash
   bash scripts/install_iqtree_static.sh
   ```
   This installs it to `tools/iqtree-2.4.0/bin/iqtree2` and
   `slurm_jobs/03_iqtree.sbatch` points at it by default.

2. **No ESMFold module**, but **ColabFold is already installed with GPU
   support**. Rather than custom-installing `torch`+`fair-esm` (which also
   needs internet access to pull multi-GB model weights, uncertain on compute
   nodes), the default path skips the ESMFold triage step entirely and runs
   ColabFold directly on everything AlphaFold DB doesn't already cover:
   ```bash
   sbatch slurm_jobs/04_colabfold_direct.sbatch
   ```
   This is the right call as long as the post-filtering homolog set stays in
   the low hundreds (`config.yaml`'s `filtering.max_hits` and
   `min_domain_coverage` control this). If it turns out much larger, the
   original ESMFold-triage + top-N-ColabFold path
   (`slurm_jobs/04_esmfold.sbatch` + `05_colabfold.sbatch`) is still there as
   a fallback, but needs that custom environment built first - not done yet,
   ask if you hit this case.

**Also worth noting from the probe:**
- Your actual home directory is `/camp/home/<username>` (142K, 5GB quota) -
  not `/nemo/homefs/home/<username>` as originally assumed. Same quota
  constraint either way; just use `/camp/home/<username>` when checking
  usage (`du -sh ~`).
- `/nemo/lab` was at 89% full (1.4TB free) when last checked - worth another
  `df -h /nemo/lab` before downloading a large reference database.

## 1. Log in and place the project correctly

```bash
ssh <username>@login.nemo.thecrick.org
```

Your home directory (`/camp/home/<username>` - confirmed via the probe
script, a legacy path still in use post the camp->nemo rename) is a 5GB
landing area only - never put the repo, a venv, or any data there. Clone into
your lab's space instead:

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

The database is an existing institute copy - nothing to download. See
`docs/runbooks/uniprot_databases_on_nemo.md` for the full survey; the short
version is `config.yaml`'s `database` field already points at
`/flask/reference/Alphafold3_dataset/datasets/uniref90_2022_05.fa`, and
`slurm_jobs/01_jackhmmer.sbatch` defaults to the same path.

```bash
sbatch slurm_jobs/01_jackhmmer.sbatch
squeue -u "$USER"
```

This is the long pole - test first on `nint` with a 2-iteration run against
a small slice of the real database (it's 71.8GB; don't try to load the
whole thing interactively - the login node also caps interactive commands
at 60s, far too short for a real search):

```bash
srun --partition=nint --cpus-per-task=4 --mem=16G --time=00:30:00 --pty /bin/bash
# inside:
ml HMMER/3.4-gompi-2023a
head -c 500M /flask/reference/Alphafold3_dataset/datasets/uniref90_2022_05.fa > /tmp/small_test.fa
jackhmmer -N 2 --cpu 4 data/raw/human_BRCA2.fasta /tmp/small_test.fa
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

## 8. Stage 5 - ColabFold for whatever AlphaFold DB didn't cover (cluster)

Default path on Nemo (see section 0): ColabFold is a ready module, so it runs
directly rather than through an ESMFold triage step.

```bash
rsync -avz data/processed/needs_prediction.fasta <username>@login.nemo.thecrick.org:.../data/processed/

# on Nemo
sbatch slurm_jobs/04_colabfold_direct.sbatch
```

If `needs_prediction.fasta` turns out to be a few hundred sequences or more,
this gets expensive - check in before submitting, since the fallback
(ESMFold triage + shortlist, `slurm_jobs/04_esmfold.sbatch` +
`05_colabfold.sbatch`) needs a custom environment built first.

## 9. Final merge and plots (laptop)

```bash
rsync -avz <username>@login.nemo.thecrick.org:.../results/structures/colabfold/ ./results/structures/colabfold/
mkdir -p results/structures/esmfold   # empty is fine - the compile script handles it

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
