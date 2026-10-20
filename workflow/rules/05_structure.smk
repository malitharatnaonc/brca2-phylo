# Mixed stage:
#   - alphafold_db_lookup   -> LAPTOP  (just REST calls)
#   - esmfold_predict       -> CLUSTER (GPU, batch over everything not in AFDB)
#   - colabfold_shortlist   -> CLUSTER (GPU, slow, only top-N candidates)

rule alphafold_db_lookup:
    input:
        "results/metrics/homolog_metrics.tsv"
    output:
        covered="data/external/afdb_hits.tsv",
        missing="data/processed/needs_prediction.fasta"
    script:
        "../../scripts/afdb_lookup.py"

rule esmfold_predict:
    input:
        "data/processed/needs_prediction.fasta"
    output:
        directory("results/structures/esmfold/")
    threads: 8
    resources:
        mem_mb=64000,
        runtime=360,
        partition=config["slurm"]["gpu_partition"],
        slurm_extra="--gres=gpu:1"
    log:
        "logs/esmfold.log"
    shell:
        # ESMFOLD_SETUP must be confirmed via scripts/probe_cluster_environment.sh
        # first - it's unlikely to be a plain module (see slurm_jobs/04_esmfold.sbatch).
        # This rule is a secondary path; the primary one is submitting
        # slurm_jobs/04_esmfold.sbatch directly, per the runbook.
        """
        ml purge
        eval "${{ESMFOLD_SETUP:?Set ESMFOLD_SETUP - see docs/runbooks/brca2_pipeline_runbook.md}}"
        python scripts/run_esmfold_batch.py {input} {output} \
            --batch-size {config[structure][esmfold_batch_size]} > {log} 2>&1
        """

rule shortlist_for_colabfold:
    input:
        metrics="results/metrics/homolog_metrics.tsv",
        esmfold="results/structures/esmfold/"
    output:
        "data/processed/colabfold_shortlist.fasta"
    params:
        n=config["structure"]["colabfold_shortlist_size"]
    script:
        "../../scripts/build_shortlist.py"

rule colabfold_predict:
    input:
        "data/processed/colabfold_shortlist.fasta"
    output:
        directory("results/structures/colabfold/")
    threads: 8
    resources:
        mem_mb=64000,
        runtime=720,
        partition=config["slurm"]["gpu_partition"],
        slurm_extra="--gres=gpu:1"
    log:
        "logs/colabfold.log"
    shell:
        # Secondary path - see slurm_jobs/05_colabfold.sbatch for the primary one.
        """
        ml purge
        eval "${{COLABFOLD_SETUP:?Set COLABFOLD_SETUP - see docs/runbooks/brca2_pipeline_runbook.md}}"
        colabfold_batch {input} {output} > {log} 2>&1
        """

rule compile_structure_metrics:
    input:
        afdb="data/external/afdb_hits.tsv",
        esmfold="results/structures/esmfold/",
        colabfold="results/structures/colabfold/"
    output:
        "results/structures/structure_metrics.tsv"
    script:
        "../../scripts/compile_structure_metrics.py"
