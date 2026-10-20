"""
Annotate hits with organism/taxonomy via the UniProt REST API.
Run locally (laptop) - network calls, not compute.

Usage:
    python fetch_taxonomy.py <hits_raw.tsv> <out.tsv>
"""
import sys
import time
import requests
import pandas as pd

UNIPROT_API = "https://rest.uniprot.org/uniprotkb/{}.json"


def fetch_one(accession: str) -> dict:
    try:
        resp = requests.get(UNIPROT_API.format(accession), timeout=15)
        resp.raise_for_status()
        data = resp.json()
        organism = data.get("organism", {})
        return {
            "uniprot_id": accession,
            "organism_name": organism.get("scientificName"),
            "taxon_id": organism.get("taxonId"),
            "lineage": ";".join(organism.get("lineage", [])),
            "sequence_length": data.get("sequence", {}).get("length"),
        }
    except Exception as e:
        return {"uniprot_id": accession, "organism_name": None,
                "taxon_id": None, "lineage": None, "sequence_length": None,
                "fetch_error": str(e)}


def fetch_taxonomy(hits_df: pd.DataFrame, sleep: float = 0.1) -> pd.DataFrame:
    accessions = hits_df["uniprot_id"].unique()
    records = []
    for i, acc in enumerate(accessions):
        records.append(fetch_one(acc))
        if i % 50 == 0:
            print(f"  fetched {i}/{len(accessions)}")
        time.sleep(sleep)  # be polite to the API
    tax_df = pd.DataFrame(records)
    return hits_df.merge(tax_df, on="uniprot_id", how="left")


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
