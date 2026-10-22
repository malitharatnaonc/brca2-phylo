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


def aggregate_coverage(df: pd.DataFrame) -> pd.DataFrame:
    """Per uniprot_id: union span of hmm_from/hmm_to across all domain hits,
    as a rough proxy for total query coverage (not exact union-of-intervals,
    but fine as a first filter - tighten later if you need precision)."""
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
    agg = agg.sort_values(["approx_coverage", "best_evalue"], ascending=[False, True])
    return agg.groupby("taxon_id", dropna=False, group_keys=False).head(per_species_n)


def filter_hits(df: pd.DataFrame, min_cov: float, max_hits: int,
                 per_species_n: int | None = None) -> pd.DataFrame:
    agg = aggregate_coverage(df)
    filtered = agg[agg["approx_coverage"] >= min_cov].copy()
    if per_species_n:
        filtered = select_top_per_species(filtered, per_species_n)
    filtered = filtered.sort_values("approx_coverage", ascending=False)
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
