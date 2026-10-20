# Run on the CLUSTER. Heaviest single step - big DB, multi-iteration search.
# Typically launched via slurm_jobs/01_jackhmmer.sbatch directly, but wired
# in here too so the whole thing can run under `snakemake --profile slurm`.

rule jackhmmer_search:
    input:
        query=config["query_fasta"],
        db=config["database"]
    output:
        domtbl="results/hits/brca2.domtbl",
        tbl="results/hits/brca2_tbl.txt",
        sto="results/hits/brca2_align.sto"
    threads: config["jackhmmer"]["cpus"]
    resources:
        mem_mb=64000,
        runtime=720,
        partition=config["slurm"]["partition"]
    log:
        "logs/jackhmmer.log"
    shell:
        # HMMER_MODULE must be confirmed via scripts/probe_cluster_environment.sh
        # first. This rule is a secondary path; the primary one is submitting
        # slurm_jobs/01_jackhmmer.sbatch directly, per the runbook.
        """
        ml purge
        ml "${{HMMER_MODULE:?Set HMMER_MODULE - see docs/runbooks/brca2_pipeline_runbook.md}}"
        jackhmmer -N {config[jackhmmer][iterations]} \
            --incE {config[jackhmmer][incE]} \
            --cpu {threads} \
            --domtblout {output.domtbl} \
            --tblout {output.tbl} \
            -A {output.sto} \
            {input.query} {input.db} > {log} 2>&1
        """
