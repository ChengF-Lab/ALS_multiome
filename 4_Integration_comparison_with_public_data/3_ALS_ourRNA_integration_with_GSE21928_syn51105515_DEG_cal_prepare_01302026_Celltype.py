#prepare psedo-bulk data for DEGs calucaltion

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
Celltypes_updated_va=['ExN', 'ODC', 'Astro', 'OPC', 'MG', 'InN', 'T cells', 'Other vascular', 'Endo']
Celltypes_updated_va_need=['ExN', 'ODC', 'Astro', 'OPC', 'MG', 'InN', 'Endo']



################################################################################################
#################################################prepare psedobulk data################################################
################################################################################################
adata =sc.read_h5ad(os.path.join(PATH_O_data_seq,"1_2_ALS_ourRNA_merge_GSE219821_merge_syn51105515_cluster_01302026.h5ad"))  #i want to dataframe results

adata.obs["Celltype_updated_vas"]=adata.obs["Celltype"]
adata.obs["Celltype_updated_vas"] = adata.obs["Celltype_updated_vas"].replace({"Mural": "Other vascular", "Fibroblast": "Other vascular","Pericytes": "Other vascular","VLMC": "Other vascular"})

adata.obs["data_source_upgrade"]=adata.obs["data_source"]
adata.obs["data_source_upgrade"] = adata.obs["data_source_upgrade"].replace({"multiome_new": "This work", "extra_snRNA-seq": "This work"})

def load_protein_coding_genes(dir_file, header = True):
    protein_coding_genes = []
    with open(dir_file, mode='r') as f:
        if header:
            next(f)
        for line in f:
            gene_sym, ncbi_ID, ensid, is_pc = line.strip("\n").split("\t")
            if gene_sym != '' and is_pc == 'protein-coding gene':
                protein_coding_genes.append(gene_sym)
    return protein_coding_genes

protein_coding_genes = set(load_protein_coding_genes(dir_file = os.path.join("/home/doul2/isilon/Cheng-Dou/LijunDou/Work/Workstation/lri-uapps-1/Work/database/snRNA_gene_symbol_ID_QC_Mapping/protein_coding_genes_ensembleID_and_genesym.tsv")))
print('There are totally {:,} protein coding genes.'.format(len(protein_coding_genes)))


pc_genes = set(adata.var_names) & protein_coding_genes
adata = adata[:, list(pc_genes)]


cell_meta = adata.obs
genes = adata.var_names.to_list()
celltypes = adata.obs["Celltype_updated_vas"].unique()

celltypes = adata.obs.Celltype_updated_vas.unique()
genes = adata.var_names
meta_keep_cols = [
    "Indivi_ID", "Fibrinogen_level", "Age", "disease_furation_yrs", "PMI",
    "Sex", "Region", "ALS_group_sf", "ALS_diagnosis", "ALS_type_Fibrinogen",
    "data_source", "SubCelltype", "Celltype", "Cellclass", "ALS_type",
    "IndiviID_Region_dataset", "Celltype_updated_vas", "data_source_upgrade"
]
for celltype in Celltypes_updated_va_need:
    print(f"Processing {celltype}...")
    cell_meta_sub = cell_meta[cell_meta["Celltype_updated_vas"] == celltype]
    cell_meta_sub["barcode"] = cell_meta_sub.index
    unique_samples = cell_meta_sub['IndiviID_Region_dataset'].unique()
    pseudobulk_dict = {}
    for sample_id in unique_samples:
        barcodes_sample = cell_meta_sub[cell_meta_sub['IndiviID_Region_dataset'] == sample_id]['barcode'].tolist()
        adata_sample = adata[barcodes_sample]
    print(f"Processing {celltype}...")
    cell_meta_sub = cell_meta[cell_meta["Celltype_updated_vas"] == celltype]
    cell_meta_sub["barcode"] = cell_meta_sub.index
    unique_samples = cell_meta_sub['IndiviID_Region_dataset'].unique()
    pseudobulk_dict = {}
    for sample_id in unique_samples:
        barcodes_sample = cell_meta_sub[cell_meta_sub['IndiviID_Region_dataset'] == sample_id]['barcode'].tolist()
        adata_sample = adata[barcodes_sample]
        count_X = adata_sample.layers["counts"] if "counts" in adata_sample.layers else adata_sample.X
        if not sp.issparse(count_X):
            count_X = sp.csr_matrix(count_X)
        sample_counts = np.asarray(count_X.sum(axis=0)).ravel()
        print(celltype,sample_id,adata_sample.shape[0],sample_counts.shape)
        pseudobulk_dict[sample_id] = sample_counts
    pseudobulk = pd.DataFrame.from_dict(
        pseudobulk_dict,
        orient='index',
        columns=genes
    )
    metadata = (
        cell_meta_sub[meta_keep_cols]
        .drop_duplicates("IndiviID_Region_dataset")
        .set_index('IndiviID_Region_dataset')
    )
    metadata_aligned = metadata.loc[pseudobulk.index]
    assert all(metadata_aligned.index == pseudobulk.index), "Index mismatch between metadata and pseudobulk!"
    pseudobulk_with_meta = pd.concat([metadata_aligned, pseudobulk], axis=1)
    # Save file
    celltype_safe = celltype.replace("/", "_").replace(" ", "_")
    out_path = os.path.join(PATH_O_data_DEGs,"Psedobulk",
                            f"3_integrated_ALS_DEGs_Step1_pseudobulk_counts_Celltype_updated_vas_{celltype_safe}_01302026.tsv")
    pseudobulk_with_meta.to_csv(out_path, sep="\t")
    print(f"{celltype} is done. Samples: {len(pseudobulk.index)}")

