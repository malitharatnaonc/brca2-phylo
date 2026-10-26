"""
Parse a HMMER domtblout file into a tidy TSV.
Run locally (laptop) - this is pure parsing, no heavy compute.

Usage (standalone):
    python parse_domtbl.py <in.domtbl> <out.tsv>

Also callable from Snakemake via the `script:` directive, which injects a
`snakemake` object with .input/.output attributes.
"""
import re
import sys
import pandas as pd

# UniRef90/50/100 headers carry organism and taxon ID directly in the
# description (e.g. "... Tax=Triplophysa tibetana TaxID=1572043 RepID=...").
# Extracting this locally avoids ~8000 UniProt REST calls that were failing
# for roughly half of all hits (rate limiting and/or UniParc-only entries,
# never fully diagnosed) - and more importantly, avoids a real correctness
# bug: select_top_per_species() in filter_hits.py groups by taxon_id, so
# unresolved taxon_id was silently collapsing many genuinely distinct
# species into one shared "unknown" bucket during per-species dedup.
TAX_PATTERN = re.compile(r"Tax=(.+?)\s+TaxID=(\d+)")

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

    # ID extraction depends on which database this searched. Two formats in
    # play on Nemo (see docs/runbooks/uniprot_databases_on_nemo.md):
    #   - UniProtKB-style:  sp|P51587|BRCA2_HUMAN           -> accession = P51587
    #   - UniRef90-style:   UniRef90_A0A5A9P0L4              -> accession = A0A5A9P0L4
    # The UniRef cluster ID is normally the accession of the cluster's
    # representative sequence, so stripping the prefix recovers a usable
    # UniProtKB accession in the common case. Occasionally the representative
    # is a UniParc entry (ID starts with "UPI") rather than a UniProtKB
    # accession - the UniProtKB REST API calls in fetch_taxonomy.py and
    # filter_hits.py will fail cleanly for those (caught, logged, dropped),
    # not crash, but it does mean a handful of hits can be silently lost if
    # their only representative is UniParc-only.
    def extract_id(name: str) -> str:
        for prefix in ("UniRef100_", "UniRef90_", "UniRef50_"):
            if name.startswith(prefix):
                return name[len(prefix):]
        if "|" in name:
            parts = name.split("|")
            if len(parts) >= 2:
                return parts[1]
        return name

    df["uniprot_id"] = df["target_name"].apply(extract_id)

    # Extract organism/taxon_id locally where the header provides it (see
    # TAX_PATTERN comment above). Rows where this doesn't match (non-UniRef
    # databases, or any malformed entry) are left NaN - fetch_taxonomy.py
    # falls back to the UniProt API only for those, a much smaller set.
    tax_matches = df["description"].str.extract(TAX_PATTERN)
    df["organism_name"] = tax_matches[0]
    df["taxon_id"] = pd.to_numeric(tax_matches[1], errors="coerce")

    # tlen is the hit's own full sequence length, already in domtblout -
    # no need for a separate API call to get this.
    df["sequence_length"] = df["tlen"]

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
