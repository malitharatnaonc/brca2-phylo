# %% [markdown]
# # BRCA2 homolog amenability analysis
# Laptop-side: merge sequence + structure metrics, rank, visualize.
# (Saved as a plain .py with cell markers - open in Jupyter via
# jupytext, or just run top to bottom.)

# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

seq_metrics = pd.read_csv("results/metrics/homolog_metrics.tsv", sep="\t")
struct_metrics = pd.read_csv("results/structures/structure_metrics.tsv", sep="\t")

df = seq_metrics.merge(struct_metrics, on="uniprot_id", how="left")

# %% [markdown]
# ## Composite amenability score
# Weight these however matches your downstream goal - crystallography vs.
# cryo-EM vs. "just want a trustworthy computational model" will weight
# length/disorder vs. pLDDT differently.

# %%
def zscore(s):
    return (s - s.mean()) / s.std(ddof=0)

df["amenability_score"] = (
    -zscore(df["length"])
    - zscore(df["pct_disorder"])
    - zscore(df["frac_disorder_long_stretch"])
    + zscore(df["mean_plddt"].fillna(df["mean_plddt"].mean()))
)

ranked = df.sort_values("amenability_score", ascending=False)
ranked.head(20)[
    ["uniprot_id", "organism_name", "length", "pct_disorder",
     "n_brc_repeats", "mean_plddt", "amenability_score"]
]

# %% [markdown]
# ## Length vs. disorder, colored by domain/structure confidence

# %%
fig, ax = plt.subplots(figsize=(8, 6))
sc = ax.scatter(
    df["length"], df["pct_disorder"] * 100,
    c=df["mean_plddt"], cmap="viridis", s=40, alpha=0.8
)
ax.set_xlabel("Length (aa)")
ax.set_ylabel("% predicted disorder")
plt.colorbar(sc, label="Mean pLDDT")
ax.set_title("BRCA2 homologs: length vs. disorder vs. structure confidence")
plt.tight_layout()
plt.savefig("results/metrics/length_vs_disorder.png", dpi=150)

# %% [markdown]
# ## Tree with metric overlay (ete3)
# Load the IQ-TREE output and annotate leaves with amenability_score /
# pct_disorder / n_brc_repeats as a heatmap strip. Fill in once the tree
# and merged table both exist.

# %%
# from ete3 import Tree, TreeStyle, faces
# t = Tree("results/tree/homologs.treefile")
# ... annotate t's leaves from `ranked`, render with a TreeStyle ...
