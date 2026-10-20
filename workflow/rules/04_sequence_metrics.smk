# Run on the LAPTOP. CPU-only, seconds per sequence even for hundreds of
# homologs (metapredict/IUPred3 are small models, not deep-learning-heavy).

rule disorder_and_composition:
    input:
        "data/processed/homologs_filtered.fasta"
    output:
        "results/metrics/disorder_composition.tsv"
    script:
        "../../scripts/disorder_composition.py"

rule brc_repeat_scan:
    input:
        fasta="data/processed/homologs_filtered.fasta",
        hmm=config["brc_repeat_hmm"]
    output:
        "results/metrics/brc_repeat_counts.tsv"
    shell:
        # Needs hmmer, so runs wherever hmmer is actually available: on the
        # laptop if installed locally (`mamba install -c bioconda hmmer`), or
        # as a quick interactive job on Nemo's `nint` partition otherwise -
        # a few hundred sequences against one small HMM is seconds of work,
        # not worth a full batch submission.
        """
        ml purge 2>/dev/null || true
        if [ -n "${{HMMER_MODULE:-}}" ]; then ml "${{HMMER_MODULE}}"; fi
        hmmsearch --domtblout results/metrics/brc_scan.domtbl {input.hmm} {input.fasta} > /dev/null
        python scripts/count_brc_repeats.py results/metrics/brc_scan.domtbl {output}
        """

rule merge_sequence_metrics:
    input:
        disorder="results/metrics/disorder_composition.tsv",
        brc="results/metrics/brc_repeat_counts.tsv",
        tax="data/processed/hits_with_taxonomy.tsv"
    output:
        "results/metrics/homolog_metrics.tsv"
    script:
        "../../scripts/merge_metrics.py"
