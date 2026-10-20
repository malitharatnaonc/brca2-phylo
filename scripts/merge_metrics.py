"""
Merge disorder/composition metrics, BRC-repeat counts, and taxonomy into
one per-homolog table. Run locally.
"""
import sys
import pandas as pd


def merge_all(disorder_path, brc_path, tax_path) -> pd.DataFrame:
    disorder = pd.read_csv(disorder_path, sep="\t")
    brc = pd.read_csv(brc_path, sep="\t")
    tax = pd.read_csv(tax_path, sep="\t")

    tax_small = tax[["uniprot_id", "organism_name", "taxon_id", "lineage"]].drop_duplicates()

    merged = disorder.merge(brc, on="uniprot_id", how="left")
    merged = merged.merge(tax_small, on="uniprot_id", how="left")

    merged["n_brc_repeats"] = merged["n_brc_repeats"].fillna(0).astype(int)

    # Quick composite "amenability" proxy at the sequence level only -
    # structure-derived terms get folded in later once AFDB/ESMFold data lands.
    merged["seq_amenability_z"] = (
        -zscore(merged["length"])
        - zscore(merged["pct_disorder"].fillna(merged["pct_disorder"].mean()))
        - zscore(merged["frac_disorder_long_stretch"].fillna(0))
    )
    return merged


def zscore(series: pd.Series) -> pd.Series:
    return (series - series.mean()) / series.std(ddof=0)


if __name__ == "__main__":
    if "snakemake" in globals():
        disorder_path = snakemake.input["disorder"]
        brc_path = snakemake.input["brc"]
        tax_path = snakemake.input["tax"]
        out_path = snakemake.output[0]
    else:
        disorder_path, brc_path, tax_path, out_path = sys.argv[1:5]

    merged = merge_all(disorder_path, brc_path, tax_path)
    merged = merged.sort_values("seq_amenability_z", ascending=False)
    merged.to_csv(out_path, sep="\t", index=False)
    print(f"Merged metrics for {len(merged)} homologs -> {out_path}")
