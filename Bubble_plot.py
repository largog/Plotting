import pandas as pd
import os
import glob
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

#Read input files

indir = '../data/summaries_sample'
files = glob.glob(os.path.join(indir,"*_cluster_summary.tsv"))  # selects all files with this format

#Data preparation

master_series = []
for file in files:
    #Extract sample name
    sample = os.path.basename(file).replace("_cluster_summary.tsv", "")

    df = pd.read_csv(file, sep="\t")

    df_subset = df[["dominant_virus", "representative_length","family", "class", "order", "total_coverm_mean"]].copy()
    df_subset["representative_identity"] = pd.to_numeric(df["representative_identity"], errors="coerce").fillna(0)
    df_subset["representative_length"] = pd.to_numeric(df["representative_length"], errors="coerce").fillna(0)
    # Add sample name
    df_subset["sample"] = sample

    master_series.append(df_subset)

# All to a master df
master_df = pd.concat(master_series, ignore_index=True)

#legend colors

ICTV = {   "ssRNA(-)": ["Chuviridae", "Lispiviridae", "Phasmaviridae",
                        "Rhabdoviridae", "Xinmoviridae"],
            "ssRNA(+)": ["Closteroviridae", "Flaviviridae"],
            "ssRNA(RT)": ["Metaviridae", "Retroviridae"],
            "dsRNA": ["Inseviridae"],
            "dsDNA": ["Nudiviridae","Phycodnaviridae"],
            "dsDNA(RT)": ["Caulimoviridae"],
            "Unasigned": ["Unasigned family"]
}
color_map = {}
palette_names = ["Reds_r", "Wistia", "ocean", "magma", "Greens_r", "summer_r", "Greys"]

for (subgroup, members), pal in zip(ICTV.items(), palette_names):
    n = len(members)
    shades = sns.color_palette(pal, n)
    for m, c in zip(members, shades):
        color_map[m] = c


plt.figure(figsize=(12, 8))
ax = sns.scatterplot(
    data=master_df,
    x="representative_identity",    # identity
    y="dominant_virus",             # dominant virus
    size="representative_length",   # length
    hue="family",
    palette=color_map,
    sizes=(20, 500),
    alpha=0.7,
    edgecolor="none"
)
#Legend position
ax.legend(bbox_to_anchor=(1.05, 1),
          loc="upper left",
          fontsize=8,
          title="Family")
#Labels
plt.xlabel("Blastx identity")
plt.ylabel("Virus")
plt.title("Viruses across samples")
plt.tight_layout()
sns.set(font="Arial")
plt.savefig("blast_identity_length.png", dpi=600)
plt.show()