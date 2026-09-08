# keep all nuclei but only perform label transfer for ODCs, Astrocytes and MG |

import scanpy as sc
import scvi
import torch
import os
import matplotlib.pyplot as plt


PATH_O = "xxx"
PATH_O_data = os.path.join(PATH_O, "Data")
PATH_O_data_seq = os.path.join(PATH_O_data, "data_seq")
PATH_O_data_scvi_results = os.path.join(PATH_O_data, "scvi_results12252025")
PATH_O_data_scvi_model = os.path.join(PATH_O_data_scvi_results, "model_save")
PATH_O_figures = os.path.join(PATH_O, "Figures")

adata = sc.read_h5ad(os.path.join(PATH_O_data_seq, "4_1_ALS_multiomeRNA_merge_extraRNA_scanvi_input_12252025_all_nuclei_mapping_ODCAstrocytesMGNeurons.h5ad"))
adata.layers["counts"] = adata.X.copy()
adata.obs["seed_label"]=adata.obs["SubCelltype_scanvi"]
sc.pp.filter_genes(adata,min_counts=None, min_cells=100, max_counts=None, max_cells=None)
adata.var["mt"] = adata.var_names.str.contains("^MT-", case=False, regex=True)
adata.var["ribo"] = adata.var_names.str.contains("^(RPS|RPL)", case=False, regex=True)
adata.var["hb"] = adata.var_names.str.contains("^HB(?!.*\(P\))", case=False, regex=True)
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt","ribo"], inplace=True)



scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="Indivi_ID",
labels_key="seed_label",
                            continuous_covariate_keys=["n_genes_by_counts", "total_counts", "pct_counts_mt"])
scvi_model = scvi.model.SCVI(adata, n_layers=2, n_latent=30)
scvi_model.train(max_epochs = 200, plan_kwargs={'lr':1e-4},devices=1)
scvi_model.save(os.path.join(PATH_O_data_scvi_model),overwrite=True)
os.rename(os.path.join(PATH_O_data_scvi_model,"model.pt"), os.path.join(PATH_O_data_scvi_model,"4_3_ALS_multiomeRNA_merge_extraRNA_scvi_label_transfer_unsupervised_all_nuclei_mapping_ODCAstrocytesMGNeurons_12252025_final.pt"))
adata.obsm["X_scVI"] = scvi_model.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=30, metric='cosine')
sc.tl.umap(adata)
adata.obsm["X_scVI_umap"]=adata.obsm["X_umap"]




#########################################################
##################scANVI model training#####################
###🔹 scANVI--> based on the scvi model to train
scanvi_model = scvi.model.SCANVI.from_scvi_model(scvi_model, unlabeled_category="Unknown")
scanvi_model.train(max_epochs=100, plan_kwargs={'lr':1e-4},devices=1)
##MODEL SAVE
scanvi_model.save(PATH_O_data_scvi_model, overwrite=True)
os.rename(os.path.join(PATH_O_data_scvi_model,"model.pt"), os.path.join(PATH_O_data_scvi_model,"4_3_ALS_multiomeRNA_merge_extraRNA_scanvi_label_transfer_semi-supervised_all_nuclei_mapping_ODCAstrocytesMGNeurons_12252025_final.pt"))


adata.obs["C_scANVI"] = scanvi_model.predict(adata)
adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata)
adata.raw = None
adata.write_h5ad(os.path.join(PATH_O_data_scvi_results,"4_3_ALS_multiomeRNA_merge_extraRNA_all_nuclei_scvi_scanvi_label_transfer_all_nuclei_mapping_ODCAstrocytesMGNeurons_12252025_final.h5ad"))
adata= sc.read_h5ad(os.path.join(PATH_O_data_scvi_results,"4_3_ALS_multiomeRNA_merge_extraRNA_all_nuclei_scvi_scanvi_label_transfer_all_nuclei_mapping_ODCAstrocytesMGNeurons_12252025_final.h5ad"))

#UMAP
sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=30, metric='cosine')
sc.tl.umap(adata)
adata.obsm["X_scANVI_umap"]=adata.obsm["X_umap"]

sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=30) #, metric='cosine')
sc.tl.umap(adata)
adata.obsm["X_scANVI_umap_metriceuclei"]=adata.obsm["X_umap"]



adata.obs["SubCelltype_merge_scANVI"] = (
    adata.obs["seed_label"].astype(str)
)
update_ct = ["Unknown"]
mask = adata.obs["seed_label"].isin(update_ct)
adata.obs.loc[
    mask,
    "SubCelltype_merge_scANVI"
] = adata.obs.loc[mask, "C_scANVI"]
adata.raw = None
adata.write_h5ad(os.path.join(PATH_O_data_scvi_results,"4_3_ALS_multiomeRNA_merge_extraRNA_all_nuclei_scvi_scanvi_label_transfer_all_nuclei_mapping_ODCAstrocytesMGNeurons_12252025_final.h5ad"))





# UMAP visulization of scVI
adata=sc.read_h5ad(os.path.join(PATH_O_data_scvi_results,"4_3_ALS_multiomeRNA_merge_extraRNA_all_nuclei_scvi_scanvi_label_transfer_all_nuclei_mapping_ODCAstrocytesMGNeurons_12252025_final.h5ad"))
adata.obsm["X_umap"]=adata.obsm["X_scVI_umap"]
plt.figure(figsize=(10, 10))
sc.pl.umap(
    adata,
    color=["Indivi_ID","data_source","Celltype"],
    frameon=False,
    size=4,
    legend_fontsize='medium')
plt.savefig(os.path.join(PATH_O_figures,"4_3_ALS_multiomeRNA_merge_extraRNA_scvi_label_transfer_unsupervised_all_nuclei_mapping_ODCAstrocytesMGNeurons_UMAP_12252025_final.png"))



adata.obsm["X_umap"]=adata.obsm["X_scANVI_umap"]
plt.figure(figsize=(15, 10))
sc.pl.umap(
    adata,
    color=["Indivi_ID","data_source","Celltype","SubCelltype","C_scANVI","SubCelltype_merge_scANVI"],
    frameon=False,
    size=4,
    legend_fontsize='medium')
plt.savefig(os.path.join(PATH_O_figures,"4_3_ALS_multiomeRNA_merge_extraRNA_scanvi_label_transfer_semi-supervised_all_nuclei_mapping_ODCAstrocytesMGNeurons_UMAP_12252025_final.png"))


adata.obsm["X_umap"]=adata.obsm["X_scANVI_umap_metriceuclei"]
plt.figure(figsize=(15, 10))
sc.pl.umap(
    adata,
    color=["Indivi_ID","data_source","Celltype","SubCelltype","C_scANVI","SubCelltype_merge_scANVI"],
    frameon=False,
    size=4,
    legend_fontsize='medium')
plt.savefig(os.path.join(PATH_O_figures,"4_3_ALS_multiomeRNA_merge_extraRNA_scanvi_label_transfer_semi-supervised_all_nuclei_mapping_ODCAstrocytesMGNeurons_UMAP_12252025_metriccosine.png")) 
