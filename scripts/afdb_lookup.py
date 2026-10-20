"""
Check AlphaFold DB for each homolog. Free structural data, zero compute -
run locally. Splits the set into "already covered" (just pull pLDDT/PAE)
vs "needs prediction" (write out a fasta for ESMFold on the cluster).
"""
import sys
import time
import requests
import pandas as pd
from Bio import SeqIO

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{}"


def check_afdb(accession: str) -> dict | None:
    try:
        resp = requests.get(AFDB_API.format(accession), timeout=15)
        if resp.status_code != 200:
            return None
        data = resp.json()
        if not data:
            return None
        entry = data[0]
        return {
            "uniprot_id": accession,
            "afdb_model_url": entry.get("pdbUrl"),
            "afdb_pae_url": entry.get("paeDocUrl"),
            "mean_plddt": entry.get("confidenceAvgLocalScore"),
        }
    except Exception:
        return None


if __name__ == "__main__":
    if "snakemake" in globals():
        metrics_path = snakemake.input[0]
        covered_out = snakemake.output["covered"]
        missing_out = snakemake.output["missing"]
        fasta_path = "data/processed/homologs_filtered.fasta"
    else:
        metrics_path, covered_out, missing_out, fasta_path = sys.argv[1:5]

    metrics = pd.read_csv(metrics_path, sep="\t")
    sequences = {r.id.split()[0]: r for r in SeqIO.parse(fasta_path, "fasta")}

    covered_rows, missing_ids = [], []
    for i, acc in enumerate(metrics["uniprot_id"]):
        result = check_afdb(acc)
        if result:
            covered_rows.append(result)
        else:
            missing_ids.append(acc)
        if i % 50 == 0:
            print(f"  checked {i}/{len(metrics)}")
        time.sleep(0.1)

    pd.DataFrame(covered_rows).to_csv(covered_out, sep="\t", index=False)

    with open(missing_out, "w") as fh:
        for acc in missing_ids:
            if acc in sequences:
                SeqIO.write(sequences[acc], fh, "fasta")

    print(f"AFDB covered: {len(covered_rows)} | needs prediction: {len(missing_ids)}")
