import pandas as pd
import plotly.graph_objects as go
import matplotlib.colors as mcolors

#Read input files

clusters = pd.read_csv("../data/global_culicidae_clusters.tsv", sep ="\t")
metadata = pd.read_excel("../metadata/metadata_species.xlsx")

#Transform tables to required data format

clusters["members_list"] = clusters["members"].str.split(",")
clusters_exp = clusters.explode("members_list", ignore_index=True)
clusters_exp["sample"] = clusters_exp["members_list"].str.split("|").str[0]
#Use genus instead of samples
clusters_exp = clusters_exp.merge(metadata[["Sample", "Genus"]],
                                  left_on="sample", right_on="Sample",
                                  how="left")
clusters_unique = clusters_exp.drop_duplicates(subset=["unique_cluster_name", "sample"])
counts_fam = clusters_unique.groupby(["sample", "family",]).size().reset_index(name="count")
counts_ord = clusters_unique.groupby(["sample", "order"]).size().reset_index(name="count")
counts_cls = clusters_unique.groupby(["sample", "class"]).size().reset_index(name="count")
complete_hierarchy = clusters_unique.groupby(["sample","class","order","family"]).size().reset_index(name="count")




counts_fam_genus = clusters_unique.groupby(["Genus", "family"]).size().reset_index(name="count")
counts_ord_genus = clusters_unique.groupby(["Genus", "order"]).size().reset_index(name="count")



#Create plotly links
nodes = pd.Series(pd.concat([counts_fam_genus["Genus"], counts_fam_genus["family"]])).unique()
#add labels to each element

#dictionary comprehension
name_to_id = {name: i for i,name in enumerate(nodes)}

sources = counts_fam_genus["Genus"].map(name_to_id)
targets = counts_fam_genus["family"].map(name_to_id)
values = counts_fam_genus["count"]

#Generate the plot

genus_colors = {
    "Aedes": "#F0E442",
    "Culex": "#0072B2",
    "Psorophora": "#D55E00",
    "Coquillettidia": "#CC79A7"
}

family_colors = {
    "Metaviridae": "#000055",
    "Nudiviridae": "#98d594",
    "Flaviviridae": "#ffaa00",
    "Unasigned family": "#959595",
    "Rhabdoviridae": "#fca082",
    "Inseviridae":"#b73779",
    "Caulimoviridae":"#7fbf66",
    "Lispiviridae":"#e32f27",
    "Xinmoviridae":"#fdd4c2",
    "Chuviridae":"#b11218",
    "Phycodnaviridae":"#d3eecd",
    "Phasmaviridae":"#fb6b4b"
}

#assigns color based on dictionary
def assign_node_color(label):
    #node is a family
    if label in family_colors:
        return family_colors[label]
    #node corresponds to a sample — infer its genus
    for genus, color in genus_colors.items():
        if genus in label:
            return color
    # Fallback color
    return "#CCCCCC"

#converts hex to a faded version of the color in rgb

def fade_hex(hex_color, alpha=0.3):
    rgb = mcolors.to_rgb(hex_color)
    return f"rgba({int(rgb[0]*255)}, {int(rgb[1]*255)}, {int(rgb[2]*255)}, {alpha})"

#asign colors

node_colors = [assign_node_color(label) for label in nodes]
source_colors = [node_colors[s] for s in sources]
link_colors = [fade_hex(c, 0.35) for c in source_colors]

# === BUILD SANKEY PLOT ===
plot = go.Figure(data=[go.Sankey(
    arrangement="snap",
    node=dict(
        pad=25,
        thickness=25,
        line=dict(color="black", width=0.5),
        color=node_colors,
        label=list(name_to_id.keys())
    ),
    link=dict(
        source=sources,
        target=targets,
        value=values,
        color=link_colors
    )
)])


plot.update_layout(
    title_text="<b>Host Genus– Viral Family Connections</b>",
    font=dict(size=20, color="black", family="Arial"),
    paper_bgcolor="white",
    plot_bgcolor="white",
    margin=dict(l=60, r=60, t=90, b=60),
    height=700,
    width=1200,
)

#Save the plot
plot.write_image("Sankey_family_genus.png", width=1000, height=900, scale=3)
plot.show()





