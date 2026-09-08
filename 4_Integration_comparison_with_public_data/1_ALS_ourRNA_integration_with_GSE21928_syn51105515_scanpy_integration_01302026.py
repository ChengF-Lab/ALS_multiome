#INTEGRATION OF OUR SNRNA-seq and GSE21928 and syn51105515

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





adata_ourRNA_with_GSE21918 = sc.read_h5ad(os.path.join(PATH_O_data_seq,"0_2_ALS_ourRNA_merge_GSE219821_01302026.h5ad"))
adata_Syn51105515 = sc.read_h5ad(os.path.join(PATH_syn51105515,"Syn51105515_clean_07252025.h5ad"))
meta_Syn51105515 = pd.read_csv(PATH_syn51105515_meta,sep="\t")
meta_Syn51105515 = meta_Syn51105515.drop(columns=["Unnamed: 0"], errors="ignore")
meta_Syn51105515 = meta_Syn51105515.drop_duplicates(subset=["Indivi_ID"]).set_index("Indivi_ID")
for col in meta_Syn51105515.columns:
    if col not in adata_Syn51105515.obs.columns:
        adata_Syn51105515.obs[col] = adata_Syn51105515.obs["Indivi_ID"].map(meta_Syn51105515[col])

adata_Syn51105515=adata_Syn51105515[adata_Syn51105515.obs.Disease_group.isin(["Normal","ALS"]),:]
adata_Syn51105515.obs["ALS_group_sf"]=adata_Syn51105515.obs["Disease_subgroup"]
adata_Syn51105515.obs["ALS_group_sf"] = adata_Syn51105515.obs["ALS_group_sf"].replace("ALS", "sALS")
adata_Syn51105515.obs["ALS_group_sf"] = adata_Syn51105515.obs["ALS_group_sf"].replace("Normal", "Control")
adata_Syn51105515.obs["ALS_group_sf"] = adata_Syn51105515.obs["ALS_group_sf"].replace("C9ALS", "fALS")

adata_Syn51105515.obs["ALS_diagnosis"]=adata_Syn51105515.obs["Disease_group"]
adata_Syn51105515.obs["ALS_diagnosis"] = adata_Syn51105515.obs["ALS_diagnosis"].replace("Normal", "Control")
adata_Syn51105515.obs["ALS_type"]=adata_Syn51105515.obs["Disease_subgroup"]
adata_Syn51105515.obs["ALS_type"] = adata_Syn51105515.obs["ALS_type"].replace("Normal", "Control")
adata_Syn51105515.obs["ALS_type"] = adata_Syn51105515.obs["ALS_type"].replace("C9ALS", "C9-ALS")
adata_Syn51105515.obs["ALS_type_Fibrinogen"]="Unknown"
adata_Syn51105515.obs["Fibrinogen_level"]="Unknown"
adata_Syn51105515.obs["disease_furation_yrs"]="Unknown"
adata_Syn51105515.obs["data_source"]="Syn51105515"


##cell type information
adata_Syn51105515.obs["Celltype"]=adata_Syn51105515.obs["L2_celltype_by_dou"]
adata_Syn51105515.obs["Celltype"] = adata_Syn51105515.obs["Celltype"].replace("Fib", "Fibroblast")
adata_Syn51105515.obs["SubCelltype"]=adata_Syn51105515.obs["L2_celltype_by_dou"]
adata_Syn51105515.obs["Cellclass"]=adata_Syn51105515.obs["Celltype"]
adata_Syn51105515.obs["Cellclass"] = adata_Syn51105515.obs["Cellclass"].replace("Endo", "Vascular")
adata_Syn51105515.obs["Cellclass"] = adata_Syn51105515.obs["Cellclass"].replace("Mural", "Vascular")
adata_Syn51105515.obs["Cellclass"] = adata_Syn51105515.obs["Cellclass"].replace("Fibroblast", "Vascular")


adata_ourRNA_with_GSE21918.obs["Cellclass"]=adata_ourRNA_with_GSE21918.obs["Celltype"]
adata_ourRNA_with_GSE21918.obs["Cellclass"] = adata_ourRNA_with_GSE21918.obs["Cellclass"].replace("Endo", "Vascular")
adata_ourRNA_with_GSE21918.obs["Cellclass"] = adata_ourRNA_with_GSE21918.obs["Cellclass"].replace("Pericytes", "Vascular")
adata_ourRNA_with_GSE21918.obs["Cellclass"] = adata_ourRNA_with_GSE21918.obs["Cellclass"].replace("VLMC", "Vascular")
adata_ourRNA_with_GSE21918.obs["Cellclass"] = adata_ourRNA_with_GSE21918.obs["Cellclass"].replace("Fibroblast", "Vascular")


cols = ["Indivi_ID",
    "Fibrinogen_level", "Age", "disease_furation_yrs", "PMI", "Sex",
    "Region", "ALS_group_sf", "ALS_diagnosis", "ALS_type_Fibrinogen",
    "data_source", "SubCelltype", "Celltype", "Cellclass", "ALS_type"
]

adata_Syn51105515.obs = adata_Syn51105515.obs.loc[:, cols].copy()
adata_ourRNA_with_GSE21918.obs=adata_ourRNA_with_GSE21918.obs.loc[:, cols].copy()

adata_merged=sc.contact(adata_Syn51105515,adata_ourRNA_with_GSE21918)
adata_merged.obs["Region"] = adata_merged.obs["Region"].replace("MC", "PMC")
adata_merged.obs["Region"] = adata_merged.obs["Region"].replace("ALS motor cortex", "PMC")
adata_merged.obs["Region"] = adata_merged.obs["Region"].replace("medial frontal cortex", "PFC")
adata_merged.obs["Region"] = adata_merged.obs["Region"].replace("motor cortex", "PMC")


adata_merged = sc.concat(
    [adata_Syn51105515, adata_ourRNA_with_GSE21918],
    join="outer",
    index_unique=None
)


age = adata_merged.obs["Age"]
age_num = pd.to_numeric(
    age.astype(str).str.strip().str.replace(r"\+$", "", regex=True),
    errors="coerce"
)
adata_merged.obs["Age"] = age_num.round(0).astype("Int32")

PMI = adata_merged.obs["PMI"]
PMI_num = pd.to_numeric(
    PMI.astype(str).str.strip().str.replace(r"\+$", "", regex=True),
    errors="coerce"
)
adata_merged.obs["PMI"] = PMI_num.round(0).astype("Int32")

adata_merged.write_h5ad(os.path.join(PATH_O_data_seq,"1_1_ALS_ourRNA_merge_GSE219821_merge_syn51105515_01302026.h5ad"))

adata_merged.layers["counts"] = adata_merged.X.copy()
adata_merged.var["mt"] = adata_merged.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata_merged, qc_vars=["mt"], inplace=True)
sc.pl.violin(
    adata_merged,
    keys=["total_counts", "n_genes_by_counts", "pct_counts_mt"],
    groupby="data_source",
    jitter=0.4,
    multi_panel=True,
    size=0
)
plt.savefig(
    os.path.join(
        PATH_O_fig,
        "1_ALS_QC_metrics_01302026.png"
    ),
    dpi=100,
    bbox_inches="tight",
)
plt.close()

min_genes = 200
min_counts = 200
sc.pp.filter_cells(adata_merged, min_genes=min_genes)
adata_merged = adata_merged[adata_merged.obs["n_genes_by_counts"] > min_counts, :].copy()
sc.pp.normalize_total(adata_merged, target_sum=1e4)
sc.pp.log1p(adata_merged)

adata_merged.obs["IndiviID_Region_dataset"] = (
    adata_merged.obs["Indivi_ID"].astype(str) + "_" +
    adata_merged.obs["Region"].astype(str) + "_" +
    adata_merged.obs["data_source"].astype(str)
)

valid_batches = adata_merged.obs['IndiviID_Region_dataset'].value_counts()
valid_batches = valid_batches[valid_batches > 50].index
adata_merged = adata_merged[adata_merged.obs['IndiviID_Region_dataset'].isin(valid_batches)].copy()
sc.pp.highly_variable_genes(adata_merged, n_top_genes=2500,batch_key="IndiviID_Region_dataset")

adata_merged.raw = adata_merged
sc.pp.scale(adata_merged, max_value=10)
sc.tl.pca(
    adata_merged,
    n_comps=50,
    svd_solver="arpack"
)
sc.pl.pca_variance_ratio(adata_merged, log=True)

ho = hm.run_harmony(
    adata_merged.obsm["X_pca"],
    adata_merged.obs,
    vars_use=["IndiviID_Region_dataset"]
)

if ho.Z_corr.shape[0] == adata_merged.shape[0]:
    adata_merged.obsm["X_pca_harmony"] =ho.Z_corr
else:
    adata_merged.obsm["X_pca_harmony"] = ho.Z_corr.T
sc.pp.neighbors(
    adata_merged,
    use_rep="X_pca_harmony",
    n_neighbors=30,
    n_pcs=30,
)

sc.tl.umap(adata_merged)
adata_merged.obsm["X_harmonly_umap"]=adata_merged.obsm["X_umap"]

for res in [0.1]:
    sc.tl.leiden(adata_merged, resolution=res, key_added=f'leiden_r{res}',n_iterations=2)

adata_merged.write_h5ad(os.path.join(PATH_O_data_seq,"1_2_ALS_ourRNA_merge_GSE219821_merge_syn51105515_cluster_01302026.h5ad"))






