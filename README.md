# BRCA2 homolog phylogenetic/structural screen

Goal: find BRCA2 homologs across taxa that are more "structurally amenable"
than human BRCA2 (shorter, less disordered, higher-confidence structure
predictions) using an iterative jackhmmer search, MSA/tree building, and a
panel of sequence/structure metrics.

**Start with `docs/runbooks/brca2_pipeline_runbook.md`** for login, storage
placement, module discovery, and the full stage-by-stage submission sequence.
This README covers the project layout and rationale; the runbook covers the
commands.

## Version-control policy

Git tracks scripts, configs, the Snakefile/rules, the SLURM job scripts, and
the runbook - everything needed to reproduce the pipeline. It does not track
generated data: `data/processed/`, `data/external/`, `results/*`, and `logs/`
are git-ignored, since they're reproducible from the pipeline and move
between the laptop and Nemo via `rsync`, not git. `data/raw/` (the original
query fasta) is the one data directory that is tracked, since it defines the
experiment rather than being derived from it.

## Layout

```
docs/runbooks/     start here - login, storage, modules, full run sequence
config/            all tunable parameters live here (one YAML), nothing
                    hardcoded in scripts or rules
workflow/
  Snakefile        entry point, just includes the rule files
  rules/           one .smk file per pipeline stage
scripts/           the actual python/bash logic called by rules
profiles/slurm/    snakemake executor profile for NeMo (SLURM)
slurm_jobs/        standalone .sbatch scripts, for stages you want to run
                   by hand before wiring them into Snakemake (or if you'd
                   rather not use Snakemake's slurm executor at all)
data/
  raw/             inputs that never get modified (query fasta, raw DBs
                    or symlinks to institute copies of UniRef90 etc.)
  external/        things pulled from APIs (AlphaFold DB models, taxonomy)
  processed/       cleaned/derived data (deduped hit tables, filtered fasta)
results/           final outputs per stage (hits, msa, tree, metrics,
                    structures) - what you'd actually look at or plot from
logs/              slurm .out/.err and snakemake logs
notebooks/         laptop-side exploration/plotting, not run on cluster
```

## Where each stage runs

| Stage                              | Where      | Why                          |
|-------------------------------------|-----------|-------------------------------|
| jackhmmer search                   | cluster   | CPU/RAM heavy, big DB         |
| hit parsing / taxonomy annotation  | laptop    | pure pandas, fast             |
| MSA (mafft)                        | cluster   | slow once hit set grows       |
| tree (iqtree)                      | cluster   | bootstrap = CPU-hours         |
| disorder / composition / pI etc.   | laptop    | CPU-light, seconds per seq    |
| AlphaFold DB lookups                | laptop    | just REST calls               |
| ESMFold / ColabFold for the rest    | cluster   | needs GPU                     |
| merge, score, plot                  | laptop    | all pandas/matplotlib         |

## Cluster specifics (Crick / "NeMo")

| Partition | Hardware                          | Limits                                  | Access |
|-----------|------------------------------------|------------------------------------------|--------|
| `ncpu`    | 128 threads, 2TB RAM/node          | 7-day max, 28000mb/cpu max               | open   |
| `nint`    | same as ncpu                       | cpu<=8, 12h max - **interactive only**   | open   |
| `nhmem`   | 4TB RAM/node                       | min 2TB request, 3-day max               | open   |
| `ga100`   | 4x A100-80GB/node                  | 7-day max, 32085mb/cpu max               | open   |
| `gh100`   | 4x H100-80GB/node                  | 7-day max                                | needs HPC approval |
| `gl40`    | 4x L40/node                        | 7-day max                                | needs HPC approval |
| `vis`     | GPU, capped at 1 GPU/user          | 16h max - visualisation, not batch work  | open   |

Implications for this pipeline:
- All CPU-heavy cluster steps (jackhmmer, mafft, iqtree) go on **`ncpu`** —
  plenty of RAM/cores there, no need for `nhmem`.
- ESMFold/ColabFold default to **`ga100`** since it's open-access; `gh100`/
  `gl40` are faster but need an access request first (its-helpdesk@crick.ac.uk
  or the #hpc Slack channel) - switch `gpu_partition` in `config.yaml` once
  you have that.
- Don't point GPU batch jobs at `vis` - it's capped at 1 GPU/user and 16h,
  meant for interactive visualisation sessions, not production inference.
- Use `nint` for interactive prototyping/debugging (`srun --partition=nint
  --cpus-per-task=4 --mem=16G --time=01:00:00 --pty /bin/bash`) before
  submitting the real batch jobs — cheap way to catch bugs without burning
  an `ncpu` allocation.
- If an `ncpu` job's predicted walltime is under 6 hours, submitting with
  `--time` under `06:00:00` gets it into the "short" reservation for
  preferential scheduling. jackhmmer's actual runtime against a full
  UniRef90/UniProtKB is hard to predict in advance — worth a short test run
  on `nint` (e.g. `-N 2` against a small DB slice) before deciding whether
  the full run qualifies.
- After any job finishes, `seff <jobID>` tells you whether you over- or
  under-requested CPU/mem — worth checking after the first jackhmmer and
  iqtree runs so later resubmissions aren't wasteful (this matters here
  since TRES limits are shared across your whole lab/account, not just you).
- Resource limits apply per Slurm **account** (your lab/STP group) as well
  as per user — if jobs are queuing unexpectedly, `squeue -u <username>` and
  `scontrol show job <jobID>` will show whether you're hitting a limit, and
  `sacct` after the fact.

## Running the pipeline

See `docs/runbooks/brca2_pipeline_runbook.md` for the full sequence (module
discovery, environment setup, and the exact command for every stage). Each
rule file in `workflow/rules/` also has stage-specific notes.
