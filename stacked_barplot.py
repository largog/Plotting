import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

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

#Rename columns
metadata['unique_name'] = metadata['Species'] + ' (' + metadata['Code'] + ')'
sp_dict = dict(zip(metadata["Sample"],metadata["unique_name"]))
df_prop = df_prop.rename(columns=sp_dict)
#reorder dataframe
df_prop = df_prop[sorted(df_prop.columns)]

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
#Assign color to viral groups
for (subgroup, members), pal in zip(ICTV.items(), palette_names):
    n = len(members)
    shades = sns.color_palette(pal, n)
    for m, c in zip(members, shades):
        color_map[m] = c

ax = df_prop.T.plot(
    kind="barh",
    stacked=True,
    figsize=(10,6),
    color=[color_map[c] for c in df_prop.index]
)


#Legend formatting
# Grouped legend

legend_handles = []
for group, members in ICTV.items():
    legend_handles.append(Patch(facecolor='none', edgecolor='none', label=group))
    for fam in members:
        if fam in df_prop.index:  # only if present in the data
            legend_handles.append(Patch(facecolor=color_map[fam], label=fam))

ax.legend(handles=legend_handles, bbox_to_anchor=(1.05, 1),
          loc="upper left", fontsize=8, title="Family")


ax.invert_yaxis()
plt.xlabel("Relative abundance")
plt.ylabel("Host species")
plt.title("Viral abundance per sample")
sns.despine()
sns.set(font="Arial")
plt.tight_layout()
plt.savefig("family_stacked_barplot1.png", dpi=600)
plt.show()