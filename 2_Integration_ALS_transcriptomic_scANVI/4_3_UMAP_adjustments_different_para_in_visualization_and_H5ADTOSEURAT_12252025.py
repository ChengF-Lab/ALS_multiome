#convert scANVI-h5ad to seurat for downstream analysis |

import scanpy as sc
import scvi
import torch
import os
import matplotlib.pyplot as plt
import pandas as pd
PATH_O = "xxx "
PATH_O_data = os.path.join(PATH_O, "Data")
PATH_O_data_seq = os.path.join(PATH_O_data, "data_seq")
PATH_O_data_scvi_results = os.path.join(PATH_O_data, "scvi_results12252025")
PATH_O_data_scvi_model = os.path.join(PATH_O_data_scvi_results, "model_save")
PATH_O_figures = os.path.join(PATH_O, "Figures")

adata = sc.read_h5ad(os.path.join(PATH_O_data_scvi_results, "4_3_ALS_multiomeRNA_merge_extraRNA_all_nuclei_scvi_scanvi_label_transfer_all_nuclei_mapping_ODCAstrocytesMGNeurons_12252025_final.h5ad"))


adata.obsm["X_umap"]=adata.obsm["X_scANVI_umap"]
umap = pd.DataFrame(
    adata.obsm["X_scANVI_umap"],
    index=adata.obs_names,
    columns=["UMAP_1", "UMAP_2"]
)
umap.to_csv(os.path.join(PATH_O_data_scvi_results,"4_4_ALS_multiomeRNA_merge_extraRNA_all_nuclei_scanvi_umap_12302025.csv"),sep="\t")
celltype = adata.obs[['C_scANVI', 'SubCelltype_merge_scANVI']]
celltype.to_csv("celltype_from_h5ad.csv")
celltype.to_csv(os.path.join(PATH_O_data_scvi_results,"4_4_ALS_multiomeRNA_merge_extraRNA_all_nuclei_scanvi_celltype_12302025.csv"),sep="\t")










