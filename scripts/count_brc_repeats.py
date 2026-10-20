"""
Count BRC-repeat (Pfam PF00634) domain hits per protein from an hmmsearch
domtblout. More BRC repeats retained usually tracks with "more RAD51-binding
capacity preserved", a useful orthology-confidence signal alongside raw
sequence similarity.

Usage: python count_brc_repeats.py <scan.domtbl> <out.tsv>
"""
import sys
import pandas as pd
from parse_domtbl import parse_domtbl, DOMTBL_COLS  # reuse the same parser


def count_brc(domtbl_path: str) -> pd.DataFrame:
    df = parse_domtbl(domtbl_path)
    counts = df.groupby("uniprot_id").agg(
        n_brc_repeats=("dom_num", "nunique"),
        best_brc_evalue=("i_evalue", "min"),
    ).reset_index()
    return counts


if __name__ == "__main__":
    domtbl_path, out_path = sys.argv[1], sys.argv[2]
    counts = count_brc(domtbl_path)
    counts.to_csv(out_path, sep="\t", index=False)
    print(f"BRC-repeat counts for {len(counts)} proteins -> {out_path}")
