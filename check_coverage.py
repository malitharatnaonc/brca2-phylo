import pandas as pd

df = pd.read_csv('data/processed/hits_with_taxonomy.tsv', sep='\t')

agg = df.groupby('uniprot_id').agg(
    organism_name=('organism_name', 'first'),
    min_from=('hmm_from', 'min'),
    max_to=('hmm_to', 'max'),
    qlen=('qlen', 'first'),
).reset_index()
agg['coverage'] = ((agg['max_to'] - agg['min_from'] + 1) / agg['qlen']).clip(upper=1.0)

print(agg['coverage'].describe())
print()

near_miss = agg[(agg['coverage'] >= 0.1) & (agg['coverage'] < 0.4)]
print("Count just below the 0.4 cutoff (0.1-0.4):", len(near_miss))
print()
print("Example organisms in that near-miss bucket:")
print(near_miss[['organism_name', 'coverage']].head(15).to_string(index=False))
