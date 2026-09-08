#clustering of integrated data

import os
import sys
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import harmonypy as hm
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
import umap
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap


np.random.seed(1234)




file_path = "xxx"
PATH_O = os.path.join(file_path, "ALS_review_comments_output11052025/6_Integration_snRNA_with_public_GSE_01052026_syn51105515_results_01302026")
PATH_O_data = os.path.join(PATH_O, "Data")
PATH_O_data_seq = os.path.join(PATH_O_data, "data_seq")
PATH_O_fig = os.path.join(PATH_O, "Figures")
PATH_O_fig_polish = os.path.join(PATH_O_fig, "Figure_polished")
PATH_O_fig_polish_data = os.path.join(PATH_O_fig_polish, "Data")
PATH_syn51105515 = "xxx/0_Data_collection_integration_results/0_data_download_preparation_results/ALS_FTD/Syn51105515/Data/data_seq"
PATH_syn51105515_meta="xxx/0_Data_collection_integration_results/0_data_download_preparation_results/ALS_FTD/Syn51105515/Data/data_seq/Syn51105515_metadata_08152025.tsv"


my36colors = [
    '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
    '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
    '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
    '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
    '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
    '#968175'
]

pal1_subcelltype2=[ "#90CAF9","#6C9AC2" ,"#69729E", "#667FE1"]
color_list_python = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
    "#5254a3", "#9c9ede", "#6b6ecf", "#de9ed6", "#e6550d",
    "#31a354", "#756bb1", "#636363", "#e7ba52", "#ad494a",
    "#9e9ac8", "#cedb9c", "#e7969c", "#6baed6", "#fd8d3c",
    "#74c476"
]


adata_merged=sc.read_h5ad(os.path.join(PATH_O_data_seq,"1_2_ALS_ourRNA_merge_GSE219821_merge_syn51105515_cluster_01302026.h5ad"))  #i want to dataframe results


#####################################################################################################################
adata_merged.obsm["X_umap"]=adata_merged.obsm["X_harmonly_umap"]
fig = sc.pl.umap(
    adata_merged,
    color=["leiden_r0.1","Region","data_source","Celltype","SubCelltype","Cellclass"],
    cmap='Blues',
    ncols=4,
    frameon=False,
    vmax='p99',
    vmin="p50",
    size=2,
    show=False,
    legend_loc="on data",
    return_fig=True
)
plt.savefig(
    os.path.join(
        PATH_O_fig,
        "2_ALS_integration_UMAP_overall_BY_SCANPY_12152025.png"
    ),
    dpi=300,
    bbox_inches="tight",
)
plt.close()






#####################################################################################################################
################################ subcluster for astrcoytes################################################
#####################################################################################################################


adata_astro = adata_merged[
    adata_merged.obs["Celltype"] == "Astro",
    :
].copy()
adata_astro.X = adata_astro.layers["counts"].copy()
sc.pp.normalize_total(adata_astro, target_sum=1e4)
sc.pp.log1p(adata_astro)


valid_batches = adata_astro.obs['IndiviID_Region_dataset'].value_counts()
valid_batches = valid_batches[valid_batches > 50].index
adata_astro = adata_astro[adata_astro.obs['IndiviID_Region_dataset'].isin(valid_batches),:].copy()
sc.pp.highly_variable_genes(adata_astro, n_top_genes=2500,batch_key="IndiviID_Region_dataset")

adata_astro.raw = adata_astro
sc.pp.scale(adata_astro, max_value=10)
sc.tl.pca(
    adata_astro,
    n_comps=50,
    svd_solver="arpack"
)
sc.pl.pca_variance_ratio(adata_astro, log=True)


ho = hm.run_harmony(
    adata_astro.obsm["X_pca"],
    adata_astro.obs,
    vars_use=["IndiviID_Region_dataset"]
)


if ho.Z_corr.shape[0] == adata_astro.shape[0]:
    adata_astro.obsm["X_pca_harmony_astro"] =ho.Z_corr
else:
    adata_astro.obsm["X_pca_harmony_astro"] = ho.Z_corr.T




sc.pp.neighbors(
    adata_astro,
    use_rep="X_pca_harmony_astro",
    n_neighbors=30,
    n_pcs=30,
)


sc.tl.umap(adata_astro)
adata_astro.obsm["X_harmonly_umap_astro"]=adata_astro.obsm["X_umap"]


for res in [0.1,0.2,0.3]:
    sc.tl.leiden(adata_astro, resolution=res, key_added=f'leiden_r{res}',n_iterations=2)




sc.pp.neighbors(
    adata_astro,
    use_rep="X_pca_harmony_astro",
    n_neighbors=30,
    n_pcs=50,
# metric='cosine'
)

sc.tl.umap(adata_astro)
adata_astro.obsm["X_harmonly_umap_astro_knnpcs50"]=adata_astro.obsm["X_umap"]

for res in [0.1]:
    sc.tl.leiden(adata_astro, resolution=res, key_added=f'leiden_r{res}_knnpcs50',n_iterations=2)

adata_astro.write_h5ad(os.path.join(PATH_O_data_seq,"2_ALS_ourRNA_merge_GSE219821_merge_syn51105515_Astrocytes_subcluster_01302026.h5ad"))



astrocyte_markers_display_clear = {
    # Homeostatic/resting astrocytes
    'Homeostatic': [
        'SLC1A2',  # GLT-1, glutamate transporter
        'SLC1A3',  # GLAST
        "GRM3",
        "ARHGAP24",'CTNNA2',"WIF1"
    ],
    'Pan-reactive': [
        'GFAP',  # Glial fibrillary acidic protein (most classic
        'CD44',  # Cell adhesion molecule
        'CRYAB',  # Crystallin alpha B
        "C3",
        "DPP10", "DCLK1", "AQP1",  # "UBASH3B","ADGRV1" #,"CISH"
    ],
    'DAA': [
        'VIM',  # Vimentin
        "CHI3L1", "OSMR", 'SERPINA3',  # Serpin family A member 3
    ],
}


sc.pl.dotplot(
    adata_astro,
    var_names=astrocyte_markers_display_clear,
    groupby='leiden_r0.2',
    dendrogram=True,         # Reorders clusters by similarity
    standard_scale='var',    # Normalizes 0-1 per gene
    cmap='RdYlBu_r',         # Red-Yellow-Blue palette
swap_axes=False,
    # figsize=(3.5, 6)         # Wide figure because you have many genes
figsize = (8,2)
)

plt.savefig(
    os.path.join(
        PATH_O_fig,
        "2_ALS_integration_Astrocytes_subcluster_markers_dotplot_panreactive01302026.png"
    ),
    dpi=300,
    bbox_inches="tight",
    # bbox_extra_artists=(legend1,)  # Ensure the legend object is included in the saved area
)
plt.close()


meta=adata_astro.obs
meta.to_csv(os.path.join(PATH_O_data_seq,"2_ALS_ourRNA_merge_GSE219821_merge_syn51105515_Astrocytes_subcluster_01302026_metadata.tsv"),sep="\t")

