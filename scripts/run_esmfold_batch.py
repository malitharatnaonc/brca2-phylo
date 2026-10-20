"""
Batch ESMFold structure prediction. Run on the CLUSTER (GPU node).

Writes one <uniprot_id>.pdb and one <uniprot_id>.json (mean pLDDT etc.)
per sequence into the output directory, so downstream scripts can glob
over them without needing to re-parse PDB files.

Usage:
    python run_esmfold_batch.py <in.fasta> <out_dir> --batch-size 8
"""
import argparse
import json
import os
import torch
import esm
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta")
    parser.add_argument("out_dir")
    parser.add_argument("--batch-size", type=int, default=8)
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    model = esm.pretrained.esmfold_v1()
    model = model.eval().cuda()
    model.set_chunk_size(128)  # trade memory for speed as needed

    records = list(SeqIO.parse(args.fasta, "fasta"))
    print(f"Predicting {len(records)} structures...")

    for i, record in enumerate(records):
        uniprot_id = record.id.split()[0]
        seq = str(record.seq).replace("*", "")
        try:
            with torch.no_grad():
                output = model.infer_pdb(seq)
            pdb_path = os.path.join(args.out_dir, f"{uniprot_id}.pdb")
            with open(pdb_path, "w") as fh:
                fh.write(output)

            # per-residue pLDDT is embedded in the B-factor column of the PDB;
            # a quick re-parse keeps this script self-contained.
            plddts = [
                float(line[60:66]) for line in output.splitlines()
                if line.startswith("ATOM") and line[12:16].strip() == "CA"
            ]
            mean_plddt = sum(plddts) / len(plddts) if plddts else None
            frac_gt70 = (
                sum(1 for p in plddts if p > 70) / len(plddts) if plddts else None
            )

            meta_path = os.path.join(args.out_dir, f"{uniprot_id}.json")
            with open(meta_path, "w") as fh:
                json.dump({
                    "uniprot_id": uniprot_id,
                    "mean_plddt": mean_plddt,
                    "frac_plddt_gt70": frac_gt70,
                }, fh)

        except Exception as e:
            print(f"  FAILED on {uniprot_id}: {e}")

        if i % args.batch_size == 0:
            print(f"  {i}/{len(records)} done")

    print("ESMFold batch complete.")


if __name__ == "__main__":
    main()
