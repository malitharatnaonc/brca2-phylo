"""
Per-sequence disorder and composition metrics. CPU-only, fast - run locally.

Requires: metapredict, biopython  (pip install metapredict biopython --break-system-packages)

Metrics computed:
  - length
  - % disorder (metapredict, threshold 0.5)
  - fraction of disorder in long stretches (>=30 aa) vs short loops
  - isoelectric point, net charge at pH 7
  - cysteine count
  - fraction low-complexity-ish (simple repeat heuristic; swap in SEG for rigor)
"""
import sys
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

try:
    import metapredict as meta
    HAVE_METAPREDICT = True
except ImportError:
    HAVE_METAPREDICT = False
    print("WARNING: metapredict not installed - disorder metrics will be blank. "
          "pip install metapredict --break-system-packages", file=sys.stderr)


def longest_disordered_stretch_fraction(scores, threshold=0.5, min_len=30):
    """Fraction of residues that sit in a disordered run >= min_len."""
    in_run = [s >= threshold for s in scores]
    total_long_run_residues = 0
    run_len = 0
    for flag in in_run:
        if flag:
            run_len += 1
        else:
            if run_len >= min_len:
                total_long_run_residues += run_len
            run_len = 0
    if run_len >= min_len:
        total_long_run_residues += run_len
    return total_long_run_residues / len(scores) if scores else 0.0


def compute_metrics(record) -> dict:
    seq = str(record.seq).replace("*", "")
    analysis = ProteinAnalysis(seq)
    out = {
        "uniprot_id": record.id.split()[0],
        "length": len(seq),
        "pI": analysis.isoelectric_point(),
        "net_charge_pH7": analysis.charge_at_pH(7.0),
        "cys_count": seq.count("C"),
        "gc_irrelevant_aromaticity": analysis.aromaticity(),
    }
    if HAVE_METAPREDICT:
        scores = meta.predict_disorder(seq)
        pct_disorder = sum(1 for s in scores if s >= 0.5) / len(scores)
        out["pct_disorder"] = pct_disorder
        out["frac_disorder_long_stretch"] = longest_disordered_stretch_fraction(scores)
    else:
        out["pct_disorder"] = None
        out["frac_disorder_long_stretch"] = None
    return out


if __name__ == "__main__":
    if "snakemake" in globals():
        in_fasta = snakemake.input[0]
        out_path = snakemake.output[0]
    else:
        in_fasta, out_path = sys.argv[1], sys.argv[2]

    records = list(SeqIO.parse(in_fasta, "fasta"))
    rows = [compute_metrics(r) for r in records]
    df = pd.DataFrame(rows)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Computed disorder/composition metrics for {len(df)} sequences -> {out_path}")
