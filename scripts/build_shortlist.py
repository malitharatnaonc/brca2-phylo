"""
After ESMFold's fast first-pass screen, pick the top-N candidates (by
combined sequence + ESMFold-confidence score) to send through full
ColabFold/AF2 - the slow, expensive, high-quality step. Run locally.
"""
import sys
import glob
import json
import pandas as pd
from Bio import SeqIO


def load_esmfold_plddt(esmfold_dir: str) -> pd.DataFrame:
    rows = []
    for meta_path in glob.glob(f"{esmfold_dir}/*.json"):
        with open(meta_path) as fh:
            data = json.load(fh)
        rows.append({
            "uniprot_id": data["uniprot_id"],
            "esmfold_mean_plddt": data["mean_plddt"],
        })
    return pd.DataFrame(rows)


if __name__ == "__main__":
    if "snakemake" in globals():
        metrics_path = snakemake.input["metrics"]
        esmfold_dir = snakemake.input["esmfold"]
        out_path = snakemake.output[0]
        n = snakemake.params["n"]
        fasta_path = "data/processed/homologs_filtered.fasta"
    else:
        metrics_path, esmfold_dir, out_path, n, fasta_path = sys.argv[1:6]
        n = int(n)

    metrics = pd.read_csv(metrics_path, sep="\t")
    esm = load_esmfold_plddt(esmfold_dir)
    merged = metrics.merge(esm, on="uniprot_id", how="left")

    merged["combined_score"] = (
        merged["seq_amenability_z"].fillna(merged["seq_amenability_z"].min())
        + merged["esmfold_mean_plddt"].fillna(0) / 100.0
    )
    shortlist_ids = merged.sort_values("combined_score", ascending=False)["uniprot_id"].head(n)

    sequences = {r.id.split()[0]: r for r in SeqIO.parse(fasta_path, "fasta")}
    with open(out_path, "w") as fh:
        for acc in shortlist_ids:
            if acc in sequences:
                SeqIO.write(sequences[acc], fh, "fasta")

    print(f"Shortlisted {len(shortlist_ids)} candidates for ColabFold -> {out_path}")
