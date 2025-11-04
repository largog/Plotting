import pandas as pd
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import glob
import os

#Sample_names

#Read metadata

metadata = pd.read_excel('../metadata/metadata_species.xlsx')
df_metadata = metadata
indir = '../data/summaries_sample'

files = glob.glob(os.path.join(indir,"*_cluster_summary.tsv"))


# Initialize a list for all dataframes
dfs = []
for file in files:
    sample_name = os.path.basename(file).replace("_cluster_summary.tsv", "")
    df = pd.read_csv(file, sep="\t", usecols=["dominant_virus", "total_coverm_mean"])
    df = df.groupby("dominant_virus", as_index=False)["total_coverm_mean"].sum()
    df.rename(columns={"dominant_virus": "name", "total_coverm_mean": sample_name}, inplace=True)
    dfs.append(df)

# merge on 'name'
from functools import reduce
df_merged = reduce(lambda left, right: pd.merge(left, right, on="name", how="outer"), dfs)

# Replace NaN with zeros
df_merged.fillna(0, inplace=True)
#Calculate relative abundances
df_rel = df_merged.copy()
df_rel.iloc[:, 1:] = df_rel.iloc[:, 1:].div(df_rel.iloc[:, 1:].sum(axis=0), axis=1)

sp_palette = sns.color_palette("colorblind", n_colors=df_metadata['Species'].nunique())
sp_colors = dict(zip(df_metadata['Species'].unique(), sp_palette))

#Rename columns
sp_dict = dict(zip(metadata["Sample"], metadata["Code"]))
df_rel = df_rel.rename(columns=sp_dict)

# Color Species
metadata_ordered = df_metadata.set_index("Code").loc[df_rel.columns[1:]]
metadata_ordered = metadata_ordered.rename(columns={'Species': 'Host species'})
col_colors = metadata_ordered["Host species"].map(sp_colors)


#correct height
n_viruses = df_rel.shape[0]
height = max(8, n_viruses * 0.3)

cluster = sns.clustermap(
    df_rel.set_index("name"),
    cmap="rocket_r",
    metric="correlation",
    method="average",
    standard_scale=1,
    col_cluster=True,
    row_cluster=True,
    col_colors=col_colors,
    figsize=(10, height),
    dendrogram_ratio=(0.01, 0.05)
)
cluster.ax_heatmap.collections[0].colorbar.set_label("Relative abundance", fontsize=8)

# Legend species
species_present = metadata_ordered["Host species"].unique()
for sp in species_present:
    cluster.ax_col_dendrogram.bar(0, 0, color=sp_colors[sp], label=sp, linewidth=0)

# Move legend outside
cluster.ax_col_dendrogram.legend(
    title='Host species',
    loc='center left',
    bbox_to_anchor=(1.02, 0.5),  # moves legend outside to the right
    borderaxespad=0,
    fontsize=8,
    title_fontsize=9,
    frameon=False
)

plt.subplots_adjust(right=0.75)  # leave space on the right for legend
plt.tight_layout()
plt.savefig("virus_blastx_heatmap.png", dpi=600, bbox_inches="tight")
plt.show()