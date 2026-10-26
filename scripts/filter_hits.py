"""
Collapse per-domain hits to per-protein records, filter by domain coverage
(better signal than raw E-value for a protein as divergent as BRCA2), cap
the hit count, and fetch full sequences for whatever survives.

Run locally (laptop).
"""
import sys
import time
import requests
import pandas as pd

UNIPROT_FASTA_API = "https://rest.uniprot.org/uniprotkb/{}.fasta"

# Core conserved domains of human BRCA2 (P51587, 3418 aa), from UniProt/
# InterPro feature annotations: 8 BRC repeats (RAD51-binding) + the
# C-terminal DNA-binding domain (helical + OB-fold domains + Tower).
# ~989 of 3418 residues (~29%) - the rest is the intrinsically disordered
# majority of the protein that motivated this whole project. A hit's
# `approx_coverage` (span relative to the FULL 3418-residue query) penalises
# any hit that doesn't also align those disordered regions - which, for a
# genuinely divergent ortholog, it correctly won't, since disordered linkers
# don't carry alignable sequence signal even in a real ortholog. That
# systematically discards exactly the distant, structurally-interesting
# homologs this project is trying to find. `core_domain_coverage` instead
# measures what fraction of the CONSERVED, FUNCTIONALLY MEANINGFUL residues
# a hit actually covers - the metric that should drive filtering.
CORE_DOMAINS = [
    (1002, 1036), (1212, 1246), (1421, 1454), (1517, 1550),
    (1664, 1696), (1837, 1869), (1972, 2005), (2051, 2084),  # BRC1-8
    (2474, 3190),  # DBD (helical + OB1 + Tower + OB2/OB3)
]
CORE_TOTAL = sum(end - start + 1 for start, end in CORE_DOMAINS)


def _interval_union(intervals: list[tuple[int, int]]) -> list[list[int]]:
    """Merge overlapping/adjacent (start, end) intervals into a minimal set."""
    if not intervals:
        return []
    merged = [list(intervals[0])]
    for start, end in sorted(intervals[1:]):
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def _overlap_length(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def _core_domain_coverage(hit_intervals: list[tuple[int, int]]) -> float:
    """Fraction of the CORE_TOTAL conserved-domain residues actually covered
    by this hit's aligned segments (one (hmm_from, hmm_to) pair per domain
    hit row for this protein - BRCA2's own multi-domain architecture means a
    real hit often legitimately has several separate aligned segments, not
    one contiguous span)."""
    hit_union = _interval_union(hit_intervals)
    covered = sum(
        _overlap_length(core_start, core_end, h_start, h_end)
        for core_start, core_end in CORE_DOMAINS
        for h_start, h_end in hit_union
    )
    return covered / CORE_TOTAL


def aggregate_coverage(df: pd.DataFrame) -> pd.DataFrame:
    """Per uniprot_id: both the naive whole-length span coverage (kept for
    reference/debugging) and the core-domain coverage (the one filtering
    should actually use - see CORE_DOMAINS comment above)."""
    grouped = df.groupby("uniprot_id").agg(
        organism_name=("organism_name", "first"),
        taxon_id=("taxon_id", "first"),
        lineage=("lineage", "first"),
        sequence_length=("sequence_length", "first"),
        n_domains=("dom_num", "count"),
        best_evalue=("i_evalue", "min"),
        min_hmm_from=("hmm_from", "min"),
        max_hmm_to=("hmm_to", "max"),
        qlen=("qlen", "first"),
    ).reset_index()
    grouped["approx_coverage"] = (
        (grouped["max_hmm_to"] - grouped["min_hmm_from"] + 1) / grouped["qlen"]
    ).clip(upper=1.0)

    core_cov = (
        df.groupby("uniprot_id")
        .apply(lambda g: _core_domain_coverage(list(zip(g["hmm_from"], g["hmm_to"]))),
               include_groups=False)
        .rename("core_domain_coverage")
    )
    grouped = grouped.merge(core_cov, on="uniprot_id")
    return grouped


def fetch_sequence(accession: str) -> str | None:
    try:
        resp = requests.get(UNIPROT_FASTA_API.format(accession), timeout=15)
        resp.raise_for_status()
        lines = resp.text.strip().splitlines()
        return "".join(lines[1:])  # drop header
    except Exception:
        return None


def select_top_per_species(agg: pd.DataFrame, per_species_n: int) -> pd.DataFrame:
    """Within each species (taxon_id), keep only the top `per_species_n` hits
    by coverage (best_evalue as tiebreak). Without this, a global coverage
    ranking can let a cluster of near-identical close relatives (e.g. many
    primates all with near-complete, high-coverage BRCA2) crowd out weaker
    -but-real hits from genuinely divergent, taxonomically interesting
    species - exactly the ones this project cares most about. Hits with no
    resolved taxon_id (missing taxonomy lookup) are grouped together as one
    "unknown" bucket and capped the same way, as a safety net rather than an
    ideal outcome - worth spot-checking that bucket if it ends up large."""
    agg = agg.sort_values(["core_domain_coverage", "best_evalue"], ascending=[False, True])
    return agg.groupby("taxon_id", dropna=False, group_keys=False).head(per_species_n)


def filter_hits(df: pd.DataFrame, min_cov: float, max_hits: int,
                 per_species_n: int | None = None) -> pd.DataFrame:
    agg = aggregate_coverage(df)
    filtered = agg[agg["core_domain_coverage"] >= min_cov].copy()
    if per_species_n:
        filtered = select_top_per_species(filtered, per_species_n)
    filtered = filtered.sort_values("core_domain_coverage", ascending=False)
    filtered = filtered.head(max_hits)
    return filtered


if __name__ == "__main__":
    if "snakemake" in globals():
        in_path = snakemake.input["tax_table"]
        out_fasta = snakemake.output["filtered"]
        out_table = snakemake.output["filtered_table"]
        min_cov = snakemake.params["min_cov"]
        max_hits = snakemake.params["max_hits"]
        per_species_n = snakemake.params.get("per_species_n")
    else:
        in_path, out_fasta, out_table, min_cov, max_hits = (
            sys.argv[1], sys.argv[2], sys.argv[3],
            float(sys.argv[4]), int(sys.argv[5])
        )
        per_species_n = int(sys.argv[6]) if len(sys.argv) > 6 else None

    df = pd.read_csv(in_path, sep="\t")
    filtered = filter_hits(df, min_cov, max_hits, per_species_n)

    with open(out_fasta, "w") as fh:
        for i, row in filtered.iterrows():
            seq = fetch_sequence(row["uniprot_id"])
            if seq:
                fh.write(f">{row['uniprot_id']} {row['organism_name']}\n{seq}\n")
            if i % 50 == 0:
                print(f"  fetched sequence {i}/{len(filtered)}")
            time.sleep(0.1)

    filtered.to_csv(out_table, sep="\t", index=False)
    print(f"Kept {len(filtered)} homologs (coverage >= {min_cov}) -> {out_fasta}")
