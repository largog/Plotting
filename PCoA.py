import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa

#Read input files

indir = '../data/summaries_sample'
metadata = pd.read_excel('../metadata/metadata_species.xlsx')
files = glob.glob(os.path.join(indir,"*_cluster_summary.tsv"))  # selects all files with this format

#Data preparation

collapsed_series = []
for file in files:
    #Extract sample name
    sample = os.path.basename(file).replace("_cluster_summary.tsv", "")

    df = pd.read_csv(file, sep="\t")

    df_subset = df[["dominant_virus", "family", "class", "order", "total_coverm_mean"]].copy()

    df_subset["total_coverm_mean"] = pd.to_numeric(df["total_coverm_mean"], errors="coerce").fillna(0)

    #Sum all family counts

    collapsed = df_subset.groupby("family")["total_coverm_mean"].sum().sort_index()
    collapsed.name = sample
    collapsed_series.append(collapsed)

#merge final
merged = pd.concat(collapsed_series, axis=1).fillna(0).sort_index()
# Normalize to proportions
df_prop = merged.div(merged.sum(axis=0), axis=1)

#calculate bray-curtis distance matrix

bray_curtis = beta_diversity("braycurtis", df_prop.T.values,ids=df_prop.columns)

ordination = pcoa(bray_curtis)

# Convert scikit_bio object to dataframe

dataframe = ordination.samples.copy()
dataframe["sample"] = dataframe.index
dataframe = dataframe.reset_index()
dataframe = pd.merge(dataframe, metadata, left_on="sample", right_on="Sample", how="inner")

#axis explanation
explained = ordination.proportion_explained * 100
print(explained)

plt.figure(figsize=(7,6))

#Set custom color palette

genus_colors = {
    "Aedes": "#F0E442",
    "Culex": "#0072B2",
    "Psorophora": "#D55E00",
    "Coquillettidia": "#CC79A7"
}

sns.scatterplot(
    data=dataframe,
    x='PC1',
    y='PC2',
    hue='Genus',
    s=200,
    palette=genus_colors,
    style="DEPARTAMENT",
    alpha=0.8,
)
#Legend parameters

plt.legend(
    title='Genus',
    bbox_to_anchor=(1.05, 1),
    loc='upper left'
)


plt.xlabel(f"PC1 ({explained['PC1']:.1f}%)")
plt.ylabel(f"PC2 ({explained['PC2']:.1f}%)")
plt.title("PCoA - Bray–Curtis distances (viral communities)")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

sns.set(font="Arial")
#Extra
sns.set_theme(style="white", context="talk")
sns.despine()
plt.tight_layout()
plt.savefig("PCoA_plot.png", dpi=600)
plt.show()