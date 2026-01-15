#Load all the dependencies
import celloracle as co
co.check_python_requirements()
from celloracle.data_conversion import seurat_object_to_anndata
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import numpy as np
import scipy.sparse as sp
import copy
import glob
import time
import os
import shutil
import sys
from tqdm.auto import tqdm
import plotly.express as px
from celloracle.applications import Pseudotime_calculator
co.__version__
print(os.getcwd())

NK_subset_adata = sc.read_h5ad("NK_subset.h5ad")

sc.pp.neighbors(NK_subset_adata, n_neighbors=4, n_pcs=20)
sc.tl.diffmap(NK_subset_adata)
# Calculate neihbors again based on diffusionmap which is not used in the CellOracle
sc.pp.neighbors(NK_subset_adata, n_neighbors=10, use_rep='X_diffmap')
sc.tl.louvain(NK_subset_adata, resolution=0.8)
# PAGA graph construction
sc.tl.paga(NK_subset_adata, groups='louvain')
plt.rcParams["figure.figsize"] = [6, 4.5]
sc.pl.paga(NK_subset_adata)
sc.tl.draw_graph(NK_subset_adata, init_pos='paga', random_state=123)
sc.pl.draw_graph(NK_subset_adata, color='louvain', legend_loc='on data')
NK_subset_adata.raw.var.rename(columns={"_index": "index_tmp"}, inplace=True)
NK_subset_adata.write_h5ad(
    "/disk2/user/yizhsu/Wanmeng/cervical_scRNA_RT/Cell_oracle/NK_subset_for_co.h5ad"
)
NK_subset_adata = sc.read_h5ad("NK_subset_for_co.h5ad")

####Pseudotime analysis
plt.rcParams["figure.figsize"] = [5,5]
plt.rcParams["savefig.dpi"] = 300
NK_subset_adata.obsm.keys()
pt = Pseudotime_calculator(adata=NK_subset_adata,
                           obsm_key="X_umap", # Dimensional reduction data name
                           cluster_column_name="NK_cell_type" # Clustering data name
                           )
sc.pl.umap(NK_subset_adata, color=['NK_cell_type'],legend_loc='on data')

print("Clustering name: ", pt.cluster_column_name)
print("Cluster list", pt.cluster_list)
pt.plot_cluster(fontsize=8)
embedding = pt.adata.obsm[pt.obsm_key]
df = pd.DataFrame(embedding, columns=["x", "y"])
df["cluster"] = pt.adata.obs[pt.cluster_column_name].values
df["label"] = pt.adata.obs.index.values
fig = px.scatter(df, x="x", y="y", hover_name=df["label"], color="cluster")
fig.write_html("plotly_scatter.html")
print("Saved plotly_scatter.html")



clusters_in_aNK_lineage = ['Adaptive_NK_cells', 'NK_cluster_6', 'NK_cluster_4', 'NK_cluster_7','NK_cluster_0','NK_cluster_1']
clusters_in_cNK_lineage = ['NK_cluster_1','NK_cluster_8','NK_cluster_9','NK_cluster_3']

# Make a dictionary
lineage_dictionary = {"Lineage_aNK": clusters_in_aNK_lineage,
           "Lineage_cNK": clusters_in_cNK_lineage}

# Input lineage information into pseudotime object
pt.set_lineage(lineage_dictionary=lineage_dictionary)

# Visualize lineage information
pt.plot_lineages()

root_cells = {"Lineage_aNK": "TGAGGGACAAATAAGC-1_12","Lineage_cNK": "TGAGGGACAAATAAGC-1_12"}
pt.set_root_cells(root_cells=root_cells)

# Check root cell and lineage
pt.plot_root_cells()

# Check diffusion map data.
"X_diffmap" in pt.adata.obsm

# Calculate pseudotime
pt.get_pseudotime_per_each_lineage()
# Check results
pt.plot_pseudotime(cmap="rainbow")
# Check result
pt.adata.obs[["Pseudotime"]].head()

# Add calculated pseudotime data to the oracle object
NK_subset_adata.obs = pt.adata.obs

# Save updated anndata object
NK_subset_adata.write_h5ad(
    "/disk2/user/yizhsu/Wanmeng/cervical_scRNA_RT/Cell_oracle/NK_subset_with_pseudo_for_co.h5ad"
)
NK_subset_adata = sc.read_h5ad("NK_subset_with_pseudo_for_co.h5ad")


#################################################
######################Start Celloracle###########

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import celloracle as co
co.__version__

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300

NK_subset_adata = sc.read_h5ad("NK_subset_with_pseudo_for_co.h5ad")

print(f"Cell number is :{NK_subset_adata.shape[0]}")
print(f"Gene number is :{NK_subset_adata.shape[1]}")

#Load the basic GRN 
base_GRN = co.data.load_human_promoter_base_GRN()
base_GRN.head()
oracle = co.Oracle()

NK_subset_adata.obs.columns
NK_subset_adata.obsm.keys()

print("Metadata columns :", list(NK_subset_adata.obs.columns))
print("Dimensional reduction: ", list(NK_subset_adata.obsm.keys()))

# In this notebook, we use the unscaled mRNA count for the nput of Oracle object.
#Load the converted h5ad file from R file with the scale data replaced by the count data
NK_subset_adata_2 = sc.read_h5ad("NK_subset_2.h5ad")
NK_subset_adata.X=NK_subset_adata_2.X.copy()
NK_subset_adata.layers["raw_count"]=NK_subset_adata_2.X.copy()

# Instantiate Oracle object.
oracle.import_anndata_as_raw_count(adata=NK_subset_adata,
                                   cluster_column_name="NK_cell_type",
                                   embedding_name="X_umap")

# You can load TF info dataframe with the following code.
oracle.import_TF_data(TF_info_matrix=base_GRN)

# Perform PCA
oracle.perform_PCA()

# Select important PCs
plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
plt.show()
print(n_comps)
n_comps = min(n_comps, 50)

n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")

k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")

oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=4)
# Save oracle object.
oracle.to_hdf5("NK_subset.celloracle.oracle")

# Load file.
oracle = co.load_hdf5("NK_subset.celloracle.oracle")

# Check clustering data
sc.pl.umap(oracle.adata, color="NK_cell_type")



########################construct link (NK_cell_type) ############################

links = oracle.get_links(cluster_name_for_GRN_unit="NK_cell_type", alpha=10,
                         verbose_level=10)

links.links_dict.keys()
links.links_dict["Adaptive_NK_cells"]

cluster = "Adaptive_NK_cells"
links.links_dict[cluster].to_csv(f"raw_GRN_for_{cluster}.csv")

links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)

plt.rcParams["figure.figsize"] = [9, 4.5]
links.plot_degree_distributions(plot_model=True,
                                               #save=f"{save_folder}/degree_distribution/",
                                               )
###Calculate the score
plt.rcParams["figure.figsize"] = [6, 4.5]
links.get_network_score()
links.merged_score.head()

# Save Links object.
links.to_hdf5(file_path="NK_subset_with_score.celloracle.links")
links = co.load_hdf5(file_path="NK_subset_with_score.celloracle.links")


#7.Network analysis; Network score for each gene
links.cluster

# Visualize top n-th genes with high scores.
plt.rcParams.update({
    "figure.figsize": (6, 8),
    "font.size": 13,
    "axes.titlesize": 16,
    "axes.labelsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 14,
    "font.family": "sans-serif",   
    "axes.linewidth": 1.2
})

###Top 30 TFs in aNK GRN
links.plot_scores_as_rank(cluster="Adaptive_NK_cells", n_gene=30, save=f"/disk2/user/yizhsu/Wanmeng/cervical_scRNA_RT/Cell_oracle/ranked_score")

##construct link (CRT_status) ############################


links_CRT = oracle.get_links(cluster_name_for_GRN_unit="category", alpha=10,
                         verbose_level=10)


links_CRT.links_dict.keys()

links_CRT.links_dict["onRT1"]
links_CRT.links_dict["onRT2"]
links_CRT.links_dict["preRT"]
cluster = "onRT1"
links_CRT.links_dict["onRT1"].to_csv(f"raw_GRN_for_onRT1.csv")
links_CRT.links_dict["onRT2"].to_csv(f"raw_GRN_for_onRT2.csv")
links_CRT.links_dict["preRT"].to_csv(f"raw_GRN_for_preRT.csv")
links_CRT.palette
links_CRT.to_hdf5(file_path="NK_subset_CRT_states.celloracle.links")
# Load links.
links_CRT = co.load_hdf5("NK_subset_CRT_states.celloracle.links")
links_CRT.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)

plt.rcParams["figure.figsize"] = [9, 4.5]

links_CRT.plot_degree_distributions(plot_model=True,
                                               #save=f"{save_folder}/degree_distribution/",
                                               )
###计算score
plt.rcParams["figure.figsize"] = [6, 4.5]
links_CRT.get_network_score()
links_CRT.merged_score.head()
# Save Links object.
links_CRT.to_hdf5(file_path="NK_subset_CRT_states.celloracle.links")

links_CRT = co.load_hdf5("NK_subset_CRT_states.celloracle.links")

#7.Network analysis; Network score for each gene
links_CRT.cluster

plt.rcParams.update({
    "figure.figsize": (6, 8),
    "font.size": 13,
    "axes.titlesize": 16,
    "axes.labelsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 14,
    "font.family": "sans-serif",  
    "axes.linewidth": 1.2
})


links_CRT.plot_scores_as_rank(cluster="onRT1", n_gene=30, save=f"/disk2/user/yizhsu/Wanmeng/cervical_scRNA_RT/Cell_oracle/ranked_score")
links_CRT.plot_scores_as_rank(cluster="onRT2", n_gene=30, save=f"/disk2/user/yizhsu/Wanmeng/cervical_scRNA_RT/Cell_oracle/ranked_score")
links_CRT.plot_scores_as_rank(cluster="preRT", n_gene=30, save=f"/disk2/user/yizhsu/Wanmeng/cervical_scRNA_RT/Cell_oracle/ranked_score")



# Compare GRN score between two states
plt.clf() 
links_CRT.plot_score_comparison_2D(value="betweenness_centrality",
                               cluster1="onRT2", cluster2="preRT",
                               percentile=98,
                              save=f"/disk2/user/yizhsu/Wanmeng/cervical_scRNA_RT/Cell_oracle/score_comparison")

plt.clf()
links_CRT.plot_score_comparison_2D(
    value="betweenness_centrality",
    cluster1="onRT2",
    cluster2="preRT",
    percentile=98,
    save=None
)
ax = plt.gca()   

for txt in ax.texts:
    x, y = txt.get_position()
    if y > x:
        txt.set_color("royalblue")   # preRT dominant
    else:
        txt.set_color("firebrick")   # onRT2 dominant
for txt in ax.texts:
    txt.set_fontsize(17)
    txt.set_fontfamily("Liberation Sans")
# 散点
sc = ax.collections[0]
offsets = sc.get_offsets()

colors = []
for x, y in offsets:
    if y > x:
        colors.append("royalblue")
    else:
        colors.append("firebrick")

sc.set_color(colors)

lims = [
    min(ax.get_xlim()[0], ax.get_ylim()[0]),
    max(ax.get_xlim()[1], ax.get_ylim()[1])
]

ax.plot(lims, lims, '--', color='black', linewidth=1)
ax.set_xlim(lims)
ax.set_ylim(lims)

ax.set_title("Betweenness centrality comparison",
             fontsize=20)

ax.set_xlabel("onRT2", fontsize=25)
ax.set_ylabel("preRT", fontsize=25)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(labelsize=25)

# 图例文字
ax.text(lims[1]*0.05, lims[1]*0.9,
        "preRT dominant", color="royalblue", fontsize=25)

ax.text(lims[1]*0.05, lims[1]*0.85,
        "onRT2 dominant", color="firebrick", fontsize=25)

plt.tight_layout()

#plt.savefig(
#    "/disk2/user/yizhsu/Wanmeng/cervical_scRNA_RT/Cell_oracle/score_comparison/score_comparison_beautified.pdf",
#    dpi=300
#)
plt.show()

##Second comparison
plt.clf() 
links_CRT.plot_score_comparison_2D(value="betweenness_centrality",
                               cluster1="onRT1", cluster2="preRT",
                               percentile=98,
                               save=f"/disk2/user/yizhsu/Wanmeng/cervical_scRNA_RT/Cell_oracle/score_comparison")
                               
plt.clf()
links_CRT.plot_score_comparison_2D(
    value="betweenness_centrality",
    cluster1="onRT1",
    cluster2="preRT",
    percentile=98,
    save=None
)
ax = plt.gca()   
for txt in ax.texts:
    x, y = txt.get_position()
    if y > x:
        txt.set_color("royalblue")   # preRT dominant
    else:
        txt.set_color("firebrick")   # onRT2 dominant
for txt in ax.texts:
    txt.set_fontsize(17)
    txt.set_fontfamily("Liberation Sans")
sc = ax.collections[0]
offsets = sc.get_offsets()

colors = []
for x, y in offsets:
    if y > x:
        colors.append("royalblue")
    else:
        colors.append("firebrick")

sc.set_color(colors)
lims = [
    min(ax.get_xlim()[0], ax.get_ylim()[0]),
    max(ax.get_xlim()[1], ax.get_ylim()[1])
]

ax.plot(lims, lims, '--', color='black', linewidth=1)
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_title("Betweenness centrality comparison",
             fontsize=20)

ax.set_xlabel("onRT1", fontsize=25)
ax.set_ylabel("preRT", fontsize=25)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(labelsize=25)

# 图例文字
ax.text(lims[1]*0.05, lims[1]*0.9,
        "preRT dominant", color="royalblue", fontsize=25)

ax.text(lims[1]*0.05, lims[1]*0.85,
        "onRT1 dominant", color="firebrick", fontsize=25)

plt.tight_layout()

#plt.savefig(
#    "/disk2/user/yizhsu/Wanmeng/cervical_scRNA_RT/Cell_oracle/score_comparison/score_comparison_beautified.pdf",
#    dpi=300
#)
plt.show()


##########################################################################################
##################In silico perturbation ##################################################
############################################################################################

import os
import sys

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import celloracle as co
co.__version__

plt.rcParams["figure.figsize"] = [6,6]
plt.rcParams["savefig.dpi"] = 600

save_folder = "/disk2/user/yizhsu/Wanmeng/cervical_scRNA_RT/Cell_oracle/figures"

os.makedirs(save_folder, exist_ok=True)

#load oracle file
oracle = co.load_hdf5("NK_subset.celloracle.oracle")
links = co.load_hdf5(file_path="NK_subset_with_score.celloracle.links")

links.filter_links()
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10,
                              use_cluster_specific_TFdict=True)

# Check gene expression
plt.clf()
goi = "PRDM1"
sc.pl.umap(oracle.adata, color=[goi, oracle.cluster_column_name],
                 layer="imputed_count", use_raw=False, cmap="viridis")
axs = sc.pl.umap(
    oracle.adata,
    color=[goi, oracle.cluster_column_name],
    layer="imputed_count",
    use_raw=False,
    cmap="viridis",
    size=40,
    alpha=0.8,
    vmax="p99",
    legend_loc="on data",
    show=False
)
ax_left, ax_right = axs[0], axs[1]

for ax in [ax_left, ax_right]:
    ax.title.set_size(20)
    ax.set_xlabel("UMAP_1", fontsize=20)
    ax.set_ylabel("UMAP_2", fontsize=20)
    ax.tick_params(axis="both", labelsize=20)
    for txt in ax.texts:
        txt.set_fontsize(14)
        txt.set_fontweight("bold")
fig = plt.gcf()
cbar_ax = min(fig.axes, key=lambda a: a.get_position().width)

cbar_ax.tick_params(labelsize=18)

plt.tight_layout()
plt.show()

#########################################
############PRDM1 KO#####################
oracle.simulate_shift(perturb_condition={"PRDM1": 0.0},
                      n_propagation=3)


# Get transition probability
oracle.estimate_transition_prob(n_neighbors=200,
                                knn_random=True,
                                sampled_fraction=1)

# Calculate embedding
oracle.calculate_embedding_shift(sigma_corr=0.05)

################################################
###############Visualization ############
###################################################
####Optimize the scale and n_grid and min_mass
fig, ax = plt.subplots(1, 2,  figsize=[13, 6])
scale = 30
# Show quiver plot
oracle.plot_quiver(scale=scale, ax=ax[0])
ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")

# Show quiver plot that was calculated with randomized graph.
oracle.plot_quiver_random(scale=scale, ax=ax[1])
ax[1].set_title(f"Randomized simulation vector")
plt.show()
# n_grid = 40 is a good starting value.
plt.clf() 
n_grid = 35
oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
# Search for best min_mass.
oracle.suggest_mass_thresholds(n_suggestion=12)
plt.show()
min_mass = 0.02
oracle.calculate_mass_filter(min_mass=min_mass, plot=True)
plt.show()


###############Plot vector fields #####
############### Plot vector fields#####
goi="PRDM1"
fig, ax = plt.subplots(1, 1, figsize=[6, 6])

scale_simulation = 33
oracle.plot_simulation_flow_on_grid(
    scale=scale_simulation,
    ax=ax
)

ax.set_title(f"Simulated cell identity shift vector: {goi} KO",fontsize=20,)
for spine in ax.spines.values():
    spine.set_linewidth(0.8)

for col in ax.collections:
    col.set_alpha(0.8)
for line in ax.lines:
    line.set_alpha(0.3)
for col in ax.collections:
    sizes = col.get_sizes()
    if sizes is not None and len(sizes) > 0:
        if np.mean(sizes) < 10:  
            col.set_alpha(0.6)    
            col.set_color("gray")
plt.show()

# Plot vector field with cell cluster
fig, ax = plt.subplots(figsize=[8, 8])

oracle.plot_cluster_whole(ax=ax, s=18)
oracle.plot_simulation_flow_on_grid(scale=33, ax=ax, show_background=False)
plt.show()

#################################################
##############Visualize pseudotime###############
fig, ax = plt.subplots(figsize=[6,6])

sc.pl.embedding(adata=oracle.adata, basis=oracle.embedding_name, ax=ax, cmap="rainbow",
                color=["Pseudotime"])

from celloracle.applications import Gradient_calculator

# Instantiate Gradient calculator object
gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
gradient.calculate_p_mass(smooth=0.8, n_grid=35, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=0.01, plot=True)

plt.show()
gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":4}, plot=True)
#gradient.transfer_data_into_grid(args={"method": "knn", "n_knn":50},plot=True)
plt.show()

# Calculate graddient
gradient.calculate_gradient()
# Show results
scale_dev = 42
fig, ax = plt.subplots(figsize=[6, 6])

gradient.plot_dev_flow_on_grid(
    scale=scale_dev,
    ax=ax
)

for spine in ax.spines.values():
    spine.set_linewidth(0.8)
for col in ax.collections:
    col.set_alpha(0.8)
for col in ax.collections:
    sizes = col.get_sizes()
    if sizes is not None and len(sizes) > 0:
        if np.mean(sizes) < 10:   
            col.set_alpha(0.7) 
            col.set_color("gray")

plt.show()


#################################################################
#############Compare the development and perturbation############

from celloracle.applications import Oracle_development_module
# Make Oracle_development_module to compare two vector field
dev = Oracle_development_module()

# Load development flow
dev.load_differentiation_reference_data(gradient_object=gradient)

# Load simulation result
dev.load_perturb_simulation_data(oracle_object=oracle)

# Calculate inner produc scores
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)


# Show perturbation scores
# Show perturbation score only
vm = 0.1

fig, ax = plt.subplots(1, 1, figsize=[6, 6])

dev.plot_inner_product_on_grid(
    vm=0.02,
    s=50,
    ax=ax
)
ax.set_title("Perturbation score", fontsize=16)
plt.show()

# Show perturbation scores with perturbation simulation vector field
fig, ax = plt.subplots(figsize=[6, 6])
dev.plot_inner_product_on_grid(vm=0.1, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=30, show_background=False, ax=ax)

plt.show()


# Let's visualize the results
dev.visualize_development_module_layout_0(s=5,
                                          scale_for_simulation=25,
                                          s_grid=50,
                                          scale_for_pseudotime=25,
                                          vm=0.1)
plt.show()

#last two plots
vm = 0.1
s_grid = 50

fig, ax = plt.subplots(figsize=(6, 6))
dev.plot_inner_product_on_pseudotime(vm=vm, s=s_grid, ax=ax)

ax.set_xlabel("Pseudotime", fontsize=20)
ax.set_ylabel("Perturbation score", fontsize=20)

ax.tick_params(axis="both", labelsize=20)

# colorbar
fig = plt.gcf()
cbar_ax = min(fig.axes, key=lambda a: a.get_position().width)
cbar_ax.tick_params(labelsize=20)

plt.tight_layout()

#out = "/disk2/user/yizhsu/Wanmeng/cervical_scRNA_RT/Cell_oracle/figures/inner_product_on_pseudotime.pdf"
#plt.savefig(out, dpi=300, bbox_inches="tight")
plt.show()

####last one plot 
import numpy as np
import matplotlib.pyplot as plt

vm = 0.1
fig, ax = plt.subplots(figsize=(7, 6))

dev.plot_inner_product_as_box(vm=vm, ax=ax)
ax.set_ylabel("Perturbation score", fontsize=20)
ax.set_xlabel("Digitized pseudotime", fontsize=20)
ax.tick_params(axis="both", labelsize=20)
ax.axhline(0, color="black", linewidth=1.2, linestyle="--", alpha=0.8)

boxes = ax.patches

h_lines = []
for l in ax.lines:
    y = l.get_ydata()
    if len(y) >= 2 and np.isclose(y[0], y[1]):
        h_lines.append(l)

for box in boxes:
    verts = box.get_path().vertices
    x_min, x_max = np.min(verts[:, 0]), np.max(verts[:, 0])
    y_min, y_max = np.min(verts[:, 1]), np.max(verts[:, 1])
    candidates = []
    for l in h_lines:
        xd = l.get_xdata()
        yd = l.get_ydata()
        y0 = yd[0]
        if (y_min <= y0 <= y_max) and (np.min(xd) <= x_max) and (np.max(xd) >= x_min):
            candidates.append(l)

    if len(candidates) == 0:
        continue
    med = max(candidates, key=lambda l: np.max(l.get_xdata()) - np.min(l.get_xdata()))
    median_y = med.get_ydata()[0]

    color = "#2F7D32" if median_y > 0 else "#B1125B" 
    box.set_facecolor(color)
    box.set_alpha(0.55)
    box.set_edgecolor("black")
    box.set_linewidth(1.2)
for line in ax.lines:
    line.set_linewidth(1.2)
    line.set_color("black")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_ylim(-0.3, 0.6)
plt.tight_layout()
plt.show()




#####################################################################################
#################Focus on a single part (aNK) to interpret the results in detail#####
from celloracle.visualizations.config import CONFIG
CONFIG["cmap_ps"] = "PiYG"
print(CONFIG)
# Get cell index list for the cells of interest: adaptive_NK_cells
clusters = ['Adaptive_NK_cells']
cell_idx = np.where(oracle.adata.obs["NK_cell_type"].isin(clusters))[0]

# Check index
print(cell_idx)

dev = Oracle_development_module()

# Load development flow
dev.load_differentiation_reference_data(gradient_object=gradient)

# Load simulation result
dev.load_perturb_simulation_data(oracle_object=oracle,
                                 cell_idx_use=cell_idx, # Enter cell id list
                                 name="lineage_aNK" # Name of this cell group. You can enter any name.
                                 )

# Calculation
dev.calculate_inner_product()
dev.calculate_digitized_ip(n_bins=10)

# Let's visualize the results
dev.visualize_development_module_layout_0(s=5,
                                          scale_for_simulation=25,
                                          s_grid=50,
                                          scale_for_pseudotime=25,
                                          vm=0.1)

plt.show()



#################Six plots
################First plot
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(6,6))
dev.plot_cluster_cells_use(ax=ax,s=20)   # 

ax.set_title("Adaptive NK cells", fontsize=20)

plt.tight_layout()
#plt.savefig(f"{outdir}/01_cluster.pdf", dpi=300, bbox_inches="tight")
plt.show()



###Second plot
fig, ax = plt.subplots(figsize=(6,6))
dev.plot_reference_flow_on_grid(scale=24, ax=ax)  
for col in ax.collections:
    sizes = col.get_sizes()
    if sizes is not None and len(sizes) > 0:
        if np.mean(sizes) < 10:   
            col.set_alpha(0.6) 
            col.set_color("gray")
plt.tight_layout()
#plt.savefig(f"{outdir}/02_development_flow.pdf", dpi=300, bbox_inches="tight")
plt.show()


#Third plot

fig, ax = plt.subplots(figsize=(6,6))
dev.plot_simulation_flow_on_grid(scale=25, ax=ax)  # 

for col in ax.collections:
    sizes = col.get_sizes()
    if sizes is not None and len(sizes) > 0:
        if np.mean(sizes) < 10:  
            col.set_alpha(0.6) 
            col.set_color("gray")
plt.tight_layout()
#plt.savefig(f"{outdir}/03_perturb_simulation.pdf", dpi=300, bbox_inches="tight")
plt.show()


#Fourth plot

fig, ax = plt.subplots(figsize=(6,6))
dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax)
dev.plot_simulation_flow_on_grid(scale=30, show_background=False, ax=ax)

plt.tight_layout()
#plt.savefig(f"{outdir}/04_inner_product_grid.pdf", dpi=300, bbox_inches="tight")
plt.show()



#Fifth plot
vm = 0.1
s_grid = 50

fig, ax = plt.subplots(figsize=(6, 6))

dev.plot_inner_product_on_pseudotime(vm=vm, s=s_grid, ax=ax)

ax.set_xlabel("Pseudotime", fontsize=20)
ax.set_ylabel("Perturbation score", fontsize=20)

ax.tick_params(axis="both", labelsize=20)

# colorbar
fig = plt.gcf()
cbar_ax = min(fig.axes, key=lambda a: a.get_position().width)
cbar_ax.tick_params(labelsize=20)

plt.tight_layout()
#out = "/disk2/user/yizhsu/Wanmeng/cervical_scRNA_RT/Cell_oracle/figures/inner_product_on_pseudotime.pdf"
#plt.savefig(out, dpi=300, bbox_inches="tight")
plt.show()



####Sixth plot
import numpy as np
import matplotlib.pyplot as plt

vm = 0.1
fig, ax = plt.subplots(figsize=(7, 6))

dev.plot_inner_product_as_box(vm=vm, ax=ax)

ax.set_ylabel("Perturbation score", fontsize=20)
ax.set_xlabel("Digitized pseudotime", fontsize=20)
ax.tick_params(axis="both", labelsize=20)

ax.axhline(0, color="black", linewidth=1.2, linestyle="--", alpha=0.8)
boxes = ax.patches

h_lines = []
for l in ax.lines:
    y = l.get_ydata()
    if len(y) >= 2 and np.isclose(y[0], y[1]):
        h_lines.append(l)

for box in boxes:
    verts = box.get_path().vertices
    x_min, x_max = np.min(verts[:, 0]), np.max(verts[:, 0])
    y_min, y_max = np.min(verts[:, 1]), np.max(verts[:, 1])
    candidates = []
    for l in h_lines:
        xd = l.get_xdata()
        yd = l.get_ydata()
        y0 = yd[0]
        if (y_min <= y0 <= y_max) and (np.min(xd) <= x_max) and (np.max(xd) >= x_min):
            candidates.append(l)

    if len(candidates) == 0:
        continue
    med = max(candidates, key=lambda l: np.max(l.get_xdata()) - np.min(l.get_xdata()))
    median_y = med.get_ydata()[0]

    color = "#2F7D32" if median_y > 0 else "#B1125B" 
    box.set_facecolor(color)
    box.set_alpha(0.55)
    box.set_edgecolor("black")
    box.set_linewidth(1.2)
for line in ax.lines:
    line.set_linewidth(1.2)
    line.set_color("black")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.set_ylim(-0.6, 0.3)

plt.tight_layout()
plt.show()



