# Run on the LAPTOP. Pure pandas/requests, fast even for thousands of hits.

rule parse_domtbl:
    input:
        domtbl="results/hits/brca2.domtbl"
    output:
        raw_table="data/processed/hits_raw.tsv"
    script:
        "../../scripts/parse_domtbl.py"

rule fetch_taxonomy:
    input:
        raw_table="data/processed/hits_raw.tsv"
    output:
        tax_table="data/processed/hits_with_taxonomy.tsv"
    script:
        "../../scripts/fetch_taxonomy.py"

rule filter_and_dedup:
    input:
        tax_table="data/processed/hits_with_taxonomy.tsv"
    output:
        filtered="data/processed/homologs_filtered.fasta",
        filtered_table="data/processed/homologs_filtered.tsv"
    params:
        min_cov=config["filtering"]["min_domain_coverage"],
        max_hits=config["filtering"]["max_hits"],
        per_species_n=config["filtering"]["per_species_top_n"]
    script:
        "../../scripts/filter_hits.py"
