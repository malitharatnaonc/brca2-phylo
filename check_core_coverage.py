"""
Recompute coverage using the core-domain metric against your EXISTING
hits_with_taxonomy.tsv - no need to re-run jackhmmer or re-fetch taxonomy,
this only needs columns already in that file.
"""
import sys
sys.path.insert(0, 'scripts')
import pandas as pd
from filter_hits import aggregate_coverage

df = pd.read_csv('data/processed/hits_with_taxonomy.tsv', sep='\t')
agg = aggregate_coverage(df)

print("Core-domain coverage distribution:")
print(agg['core_domain_coverage'].describe())
print()

print("Side-by-side: old metric vs new metric (sorted by old metric, showing")
print("the cases where they disagree most - old low, new high):")
rescued = agg[(agg['approx_coverage'] < 0.4) & (agg['core_domain_coverage'] >= 0.4)]
print(f"  {len(rescued)} hits would be RESCUED by switching (old <0.4, new >=0.4)")
print(rescued[['organism_name', 'approx_coverage', 'core_domain_coverage']].head(20).to_string(index=False))
print()

for threshold in [0.1, 0.15, 0.2, 0.25, 0.3, 0.4]:
    n = (agg['core_domain_coverage'] >= threshold).sum()
    print(f"  threshold {threshold}: {n} hits kept, {agg['taxon_id'][agg['core_domain_coverage'] >= threshold].nunique()} distinct species")
