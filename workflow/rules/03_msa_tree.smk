# Run on the CLUSTER once the homolog set is past a couple hundred sequences;
# fine to test locally on a small subset first.

rule mafft_align:
    input:
        "data/processed/homologs_filtered.fasta"
    output:
        "results/msa/homologs_aln.fasta"
    threads: config["mafft"]["cpus"]
    resources:
        mem_mb=32000,
        runtime=240,
        partition=config["slurm"]["partition"]
    log:
        "logs/mafft.log"
    shell:
        # Secondary path - see slurm_jobs/02_mafft.sbatch for the primary one.
        """
        ml purge
        ml "${{MAFFT_MODULE:?Set MAFFT_MODULE - see docs/runbooks/brca2_pipeline_runbook.md}}"
        mafft --{config[mafft][mode]} --thread {threads} {input} > {output} 2> {log}
        """

rule iqtree_build:
    input:
        "results/msa/homologs_aln.fasta"
    output:
        "results/tree/homologs.treefile"
    threads: config["iqtree"]["cpus"]
    resources:
        mem_mb=32000,
        runtime=480,
        partition=config["slurm"]["partition"]
    params:
        prefix="results/tree/homologs"
    log:
        "logs/iqtree.log"
    shell:
        # No IQ-TREE module on Nemo - uses the static binary from
        # scripts/install_iqtree_static.sh. Secondary path overall; see
        # slurm_jobs/03_iqtree.sbatch for the primary one.
        """
        "${{IQTREE_BIN:?Run scripts/install_iqtree_static.sh and set IQTREE_BIN}}" \
            -s {input} -m {config[iqtree][model_selection]} \
            -bb {config[iqtree][bootstrap]} -nt {threads} \
            -pre {params.prefix} > {log} 2>&1
        """
