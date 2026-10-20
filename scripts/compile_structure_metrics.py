"""
Pull together structure-derived metrics from whichever source each homolog
came from (AlphaFold DB / ESMFold / ColabFold) into one table: mean pLDDT,
fraction confidently folded (pLDDT>70), source flag. Run locally.
"""
import sys
import glob
import json
import pandas as pd


def load_esmfold(esmfold_dir: str) -> pd.DataFrame:
    rows = []
    for meta_path in glob.glob(f"{esmfold_dir}/*.json"):
        with open(meta_path) as fh:
            data = json.load(fh)
        rows.append({
            "uniprot_id": data["uniprot_id"],
            "mean_plddt": data["mean_plddt"],
            "frac_plddt_gt70": data.get("frac_plddt_gt70"),
            "structure_source": "esmfold",
        })
    return pd.DataFrame(rows)


def load_colabfold(colabfold_dir: str) -> pd.DataFrame:
    # ColabFold writes a *_scores.json per model; adapt the glob/keys to
    # whatever version you have installed on the cluster.
    rows = []
    for score_path in glob.glob(f"{colabfold_dir}/*_scores.json"):
        with open(score_path) as fh:
            data = json.load(fh)
        uniprot_id = score_path.split("/")[-1].split("_scores")[0]
        rows.append({
            "uniprot_id": uniprot_id,
            "mean_plddt": sum(data.get("plddt", [0])) / max(len(data.get("plddt", [1])), 1),
            "mean_pae": data.get("pae_mean"),
            "structure_source": "colabfold",
        })
    return pd.DataFrame(rows)


if __name__ == "__main__":
    if "snakemake" in globals():
        afdb_path = snakemake.input["afdb"]
        esmfold_dir = snakemake.input["esmfold"]
        colabfold_dir = snakemake.input["colabfold"]
        out_path = snakemake.output[0]
    else:
        afdb_path, esmfold_dir, colabfold_dir, out_path = sys.argv[1:5]

    afdb = pd.read_csv(afdb_path, sep="\t")
    afdb["structure_source"] = "afdb"
    afdb = afdb.rename(columns={"mean_plddt": "mean_plddt"})

    esm = load_esmfold(esmfold_dir)
    colab = load_colabfold(colabfold_dir)

    # ColabFold result (highest quality) overrides ESMFold, which overrides AFDB
    # only in the sense that we prefer the best available source per protein.
    combined = pd.concat([afdb, esm, colab], ignore_index=True)
    combined = combined.sort_values("structure_source", key=lambda s: s.map(
        {"colabfold": 0, "afdb": 1, "esmfold": 2}
    ))
    combined = combined.drop_duplicates(subset="uniprot_id", keep="first")

    combined.to_csv(out_path, sep="\t", index=False)
    print(f"Compiled structure metrics for {len(combined)} homologs -> {out_path}")
