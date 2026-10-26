"""
Fill in organism/taxonomy gaps left by parse_domtbl.py's local header
parsing, and fetch full lineage per unique SPECIES rather than per hit.
Run locally (laptop) - network calls only for what's actually missing.

Most hits from a UniRef90 search already have organism_name/taxon_id
populated locally (see parse_domtbl.py's TAX_PATTERN) - this script only
hits the UniProt API as a fallback for the remainder, and for lineage,
which isn't available in the UniRef header at all.

Usage:
    python fetch_taxonomy.py <hits_raw.tsv> <out.tsv>
"""
import sys
import time
import requests
import pandas as pd

UNIPROTKB_API = "https://rest.uniprot.org/uniprotkb/{}.json"
TAXONOMY_API = "https://rest.uniprot.org/taxonomy/{}.json"


def fetch_organism_fallback(accession: str) -> dict:
    """Per-accession UniProtKB lookup, only for hits local header parsing
    couldn't resolve (non-UniRef databases, or any malformed description)."""
    try:
        resp = requests.get(UNIPROTKB_API.format(accession), timeout=15)
        resp.raise_for_status()
        data = resp.json()
        organism = data.get("organism", {})
        return {"uniprot_id": accession,
                "organism_name": organism.get("scientificName"),
                "taxon_id": organism.get("taxonId")}
    except Exception:
        return {"uniprot_id": accession, "organism_name": None, "taxon_id": None}


def fetch_lineage(taxon_id) -> str | None:
    """One call per unique SPECIES, not per hit - the taxonomy endpoint's
    response structure is based on UniProt's documented taxonomy API (not
    independently verified live - re-check the first few real results
    against https://www.uniprot.org/taxonomy/<id> in a browser before
    trusting this at scale)."""
    try:
        resp = requests.get(TAXONOMY_API.format(int(taxon_id)), timeout=15)
        resp.raise_for_status()
        data = resp.json()
        names = [node.get("scientificName") for node in data.get("lineage", [])]
        return ";".join(n for n in names if n) or None
    except Exception:
        return None


def fetch_taxonomy(hits_df: pd.DataFrame, sleep: float = 0.1) -> pd.DataFrame:
    df = hits_df.copy()

    missing_mask = df["taxon_id"].isna()
    missing_accessions = df.loc[missing_mask, "uniprot_id"].unique()
    if len(missing_accessions):
        print(f"Local header parsing resolved {df['uniprot_id'].nunique() - len(missing_accessions)}"
              f"/{df['uniprot_id'].nunique()} accessions; falling back to the"
              f" UniProt API for the remaining {len(missing_accessions)}")
        fallback = {}
        for i, acc in enumerate(missing_accessions):
            fallback[acc] = fetch_organism_fallback(acc)
            if i % 50 == 0:
                print(f"  fallback fetch {i}/{len(missing_accessions)}")
            time.sleep(sleep)
        fallback_df = pd.DataFrame(fallback.values()).set_index("uniprot_id")
        df.loc[missing_mask, "organism_name"] = df.loc[missing_mask, "uniprot_id"].map(fallback_df["organism_name"])
        df.loc[missing_mask, "taxon_id"] = df.loc[missing_mask, "uniprot_id"].map(fallback_df["taxon_id"])
    else:
        print("Every accession resolved via local header parsing - no API fallback needed")

    unique_taxa = df["taxon_id"].dropna().unique()
    print(f"Fetching lineage for {len(unique_taxa)} unique species "
          f"(vs {df['uniprot_id'].nunique()} total hits)")
    lineage_map = {}
    for i, taxon_id in enumerate(unique_taxa):
        lineage_map[taxon_id] = fetch_lineage(taxon_id)
        if i % 50 == 0:
            print(f"  lineage fetch {i}/{len(unique_taxa)}")
        time.sleep(sleep)
    df["lineage"] = df["taxon_id"].map(lineage_map)

    return df


if __name__ == "__main__":
    if "snakemake" in globals():
        in_path = snakemake.input["raw_table"]
        out_path = snakemake.output["tax_table"]
    else:
        in_path, out_path = sys.argv[1], sys.argv[2]

    hits = pd.read_csv(in_path, sep="\t")
    annotated = fetch_taxonomy(hits)
    annotated.to_csv(out_path, sep="\t", index=False)
    print(f"Annotated {annotated['uniprot_id'].nunique()} unique accessions -> {out_path}")
