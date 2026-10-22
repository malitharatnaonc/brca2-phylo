# UniProt / UniRef reference databases on the Nemo HPC

**Surveyed:** 22 October 2026 · **Host:** `login006.nemo.thecrick.org` · **Root:** `/flask/reference`
**Method:** direct `ls`/`grep` inspection of the reference share and `module -t avail` on the login node.

All paths below are world-readable to user `mouills`. Nothing needs to be downloaded.

---

## 1. Plain FASTA — ready for jackhmmer, phmmer, or `diamond makedb`

Location: `/flask/reference/Alphafold3_dataset/datasets/`

| File | Size | Content |
|---|---|---|
| `uniref90_2022_05.fa` | 71.8 GB | UniRef90, release 2022_05 |
| `uniprot_all_2021_04.fa` | 108 GB | Full UniProtKB, release 2021_04 |
| `mgy_clusters_2022_05.fa` | 128 GB | MGnify clusters (not UniProt, listed for completeness) |
| `bfd-first_non_consensus_sequences.fasta` | 18 GB | BFD |
| `pdb_seqres_2022_09_28.fasta` | 233 MB | PDB sequence residues |

Header formats, read directly from the files:

- `uniref90_2022_05.fa` → `>UniRef90_A0A5A9P0L4 peptidylprolyl isomerase n=1 Tax=Triplophysa tibetana TaxID=1572043 RepID=A0A5A9P0L4_9TELE`
- `uniprot_sprot.fasta` (see §2) → `>sp|Q6GZX4|001R_FRG3G Putative transcription factor 001R OS=Frog virus 3 (isolate Goorha) GN=FV3-001R PE=4 SV=1`

A sharded copy of the same datasets (`<name>-00000-of-00063` … ) exists at
`/flask/reference/Alphafold3_dataset_sharded/datasets/`.

These are the databases AlphaFold 3 already searches: a 3.0.1 data-pipeline run launches four
jackhmmer subprocesses at `--cpu 8` each, against `uniref90_2022_05`, `mgy_clusters_2022_05`,
`bfd-first_non_consensus_sequences` and `uniprot_all_2021_04`.

---

## 2. Pre-formatted NCBI BLAST protein databases

Location: `/flask/reference/NCBI_BLAST/` — usable immediately with `module load BLAST+/2.17.0-gompi-2025a`.

| Database | Sequences | Built | Form |
|---|---|---|---|
| `uniref90` | 32,168,296 | Apr 2015 (build log) | multi-volume `blastp` DB (`.phr/.pin/.psq/...`) |
| `uniref50` | 10,335,613 | 2015–2016 | multi-volume DB + alias `uniref50.fasta.pal` |
| `uniref100.fasta` | — | Jul 2019 | 93 GB FASTA only, **no** BLAST index |
| `uniprot_sprot.fasta` | — | Nov 2016 | 250 MB Swiss-Prot FASTA |

Also present in the same directory: `nr`, `nr90`, `pdbaa`, `uniscop`, plus the standard scoring
matrices (BLOSUM45/62/80, PAM30/70/250) and `taxlist.txt`/`lineages.txt`.

`uniref90.fasta` carries an Easel index (`uniref90.fasta.ssi`), so single-sequence retrieval with
`esl-sfetch` works without re-indexing.

**No DIAMOND databases exist** — a `find /flask/reference -maxdepth 2 -name '*.dmnd'` returns nothing.
Build one from `uniref90_2022_05.fa` if DIAMOND is wanted.

---

## 3. Profile / MSA-format copies

| Format | Location | Contents |
|---|---|---|
| MMseqs2 | `/flask/reference/Colabfold/` | `uniref30_2302_db*` (incl. 229 GB `.idx`, `_mapping`, `_taxonomy`), `colabfold_envdb_202108_db*`, `pdb100_230517`; `UNIREF30_READY` / `COLABDB_READY` sentinel files |
| HH-suite | `/flask/reference/hhdb/uniclust30/` | `uniclust30_2018_08` a3m / hhm / cs219 ffindex set, 329 GB total |
| HH-suite | `/flask/reference/hhdb/` | also `bfd_metaclust`, `pdb70` |

---

## 4. Search tools available as modules

| Tool | Newest module |
|---|---|
| BLAST+ | `BLAST+/2.17.0-gompi-2025a` (also 2.3.0 → 2.14.0) |
| DIAMOND | `DIAMOND/2.1.14-GCC-14.2.0` |
| HMMER | `HMMER/3.4-gompi-2023a` |
| HH-suite | `HH-suite/3.3.0-gompi-2023a` |
| MMseqs2 | `MMseqs2/13-45111-gompi-2021b` (ColabFold 1.5.2 bundles its own, newer) |

---

## 5. Caveats

1. **Currency.** The BLAST-formatted UniRef databases were built in 2015–2016 and the FASTA
   `uniref100.fasta` in 2019. For anything where current sequence coverage matters, start from
   `uniref90_2022_05.fa` in the AF3 tree.
2. **No 2024+ UniProtKB.** The newest complete UniProtKB on the share is `uniprot_all_2021_04.fa`.
   A current release would have to be fetched from UniProt and indexed locally.
3. **Run searches as batch jobs.** Interactive commands on the login node are capped at 60 s, which
   is far below any real search against a 70–110 GB FASTA. Submit to `ncpu` with an explicit
   `#SBATCH --cpus-per-task=` line — `submit_job` otherwise allocates only 2 cores.
4. **Share capacity.** Check `df -h /nemo/stp` before writing large MSA or alignment output; it was
   at 99 % use in September 2026. Do not write large outputs to home (`/nemo/homefs`, 5 GB quota).

---

## 6. How this pipeline uses it

`config/config.yaml`'s `database` field points directly at
`/flask/reference/Alphafold3_dataset/datasets/uniref90_2022_05.fa` - no download needed.
`scripts/parse_domtbl.py` handles both this database's `UniRef90_<accession>` header style and
plain UniProtKB `sp|<accession>|<name>` style (in case the database is ever switched to
`uniprot_all_2021_04.fa`), stripping whichever prefix applies to recover a UniProtKB accession for
the taxonomy/sequence-fetch steps that follow.

All of this pipeline's outputs land under `/nemo/lab/<lab>/home/users/<username>/brca2-phylo/`,
not `/nemo/stp` - the 99%-full share above is a different filesystem from the one this project
actually writes to (confirmed via the probe script: `/nemo/lab` was at 89% full / 1.4TB free), but
worth re-checking occasionally since both are shared, lab-wide resources.
