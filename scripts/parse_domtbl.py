"""
Parse a HMMER domtblout file into a tidy TSV.
Run locally (laptop) - this is pure parsing, no heavy compute.

Usage (standalone):
    python parse_domtbl.py <in.domtbl> <out.tsv>

Also callable from Snakemake via the `script:` directive, which injects a
`snakemake` object with .input/.output attributes.
"""
import sys
import pandas as pd

DOMTBL_COLS = [
    "target_name", "target_accession", "tlen",
    "query_name", "query_accession", "qlen",
    "e_value_full", "score_full", "bias_full",
    "dom_num", "dom_total", "c_evalue", "i_evalue",
    "dom_score", "dom_bias",
    "hmm_from", "hmm_to", "ali_from", "ali_to",
    "env_from", "env_to", "acc", "description"
]


def parse_domtbl(path: str) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split(None, len(DOMTBL_COLS) - 1)
            rows.append(fields)
    df = pd.DataFrame(rows, columns=DOMTBL_COLS)

    numeric_cols = [
        "tlen", "qlen", "e_value_full", "score_full", "bias_full",
        "dom_num", "dom_total", "c_evalue", "i_evalue", "dom_score",
        "dom_bias", "hmm_from", "hmm_to", "ali_from", "ali_to",
        "env_from", "env_to", "acc"
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # UniProt accession is usually embedded in target_name, e.g. sp|P51587|BRCA2_HUMAN
    df["uniprot_id"] = df["target_name"].str.split("|").str[1].fillna(df["target_name"])

    # query coverage for this single domain hit (fraction of the query HMM matched)
    df["query_coverage"] = (df["hmm_to"] - df["hmm_from"] + 1) / df["qlen"]

    return df


if __name__ == "__main__":
    if "snakemake" in globals():
        domtbl_path = snakemake.input["domtbl"]
        out_path = snakemake.output["raw_table"]
    else:
        domtbl_path, out_path = sys.argv[1], sys.argv[2]

    df = parse_domtbl(domtbl_path)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Parsed {len(df)} domain hits -> {out_path}")
