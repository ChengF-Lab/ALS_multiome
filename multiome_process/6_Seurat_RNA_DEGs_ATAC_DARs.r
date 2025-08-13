library(ArchR)
library('Cairo')
addArchRGenome("hg38")
addArchRThreads(50)
addArchRLocking(locking = TRUE)
set.seed(123)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(svglite)
 library(Signac)
library(Seurat)
library(SummarizedExperiment)



input_path="xxx"
output_path="xxx"

allsamples = c("ALS1","ALS9","ALS11","ALS17","ALS3","ALS4","ALS8","ALS16")
celltypes=c("ODC","Astro","OPC","MG","Neuron","Endo")
subcelltypes=c("ODC_1",   "ODC_2",   "ODC_3","ODC_4","ODC_5","Astro_1", "Astro_2","Astro_3","OPC",
"MG_1","MG_2","Neuron","Endo")
num_celltypes=length(celltypes)
num_subcelltypes=length(subcelltypes)
cluster_inte="Clusters_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined_reso0.2"
umap_inte="UMAP_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined"
cluster_RNA="Clusters_RNA_0.2_Har_reso0.2"
umap_RNA="UMAP_LSI_RNA_0.2_Har"
cluster_atac="Clusters_ATAC_0.2_Har_reso0.2"
umap_ATAC="UMAP_LSI_ATAC_0.2_Har"
umap_RNA_beforeHar="UMAP_LSI_RNA_0.2"
umap_ATAC_beforeHar="UMAP_LSI_ATAC_0.2"



#####part II RNA_DEGs calculation
so=readRDS(file.path("xxx","xxx.rds"))
so$Celltype_Group=paste(so$SubCelltype,so$Fibrinigen,sep="_")
Idents(so)=so$Celltype_Group
DEGs_data_all=data.frame()
for (i in 1:num_subcelltypes){
    for (j in 1:num_com_groups){
celltype=subcelltypes[i]
group1=compare_group1[j]
group2=compare_group2[j]
Celltype_Group1=paste(celltype,group1,sep="_")
Celltype_Group2=paste(celltype,group2,sep="_")
DEGs <- FindMarkers(so,ident.1 =Celltype_Group2, ident.2 = Celltype_Group1,logfc.threshold=0,test.use="MAST",latent.vars=c("nCount_RNA","nFeature_RNA"))  #usually ident.)  #usually ident. vs. ident.2 , so we need to give right ident.1 and ident.2
DEGs$celltype=celltype
DEGs$group1=group2
DEGs$group2=group1
DEGs$compare_group =compare_groups[j]
DEGs$gene=rownames(DEGs)
DEGs$pct.all = DEGs$pct.1 + DEGs$pct.2
DEGS_file=file.path(output_path,paste0("6_ALS_RNA_DEGs_subcelltypes_",celltype,"_",group2,"_vs_",group1,"240920_seurat_MAST_var_nCount_nFeature.tsv"))
write.table(DEGs, file = DEGS_file, quote = FALSE, sep = "\t", col.names = NA)
if (dim(DEGs_data_all)[1]==0){DEGs_data_all=DEGs} else{DEGs_data_all=rbind(DEGs_data_all,DEGs)}
cat(celltype,"is done")
}
}


#####part II ATAC-DARs calculation
so_atac_peak=readRDS("xxx.rds")
for (i in 4:5){
    for (j in 1:num_com_groups){
celltype=subcelltypes[i]
group1=compare_group1[j]
group2=compare_group2[j]
Celltype_Group1=paste(celltype,group1,sep="_")
Celltype_Group2=paste(celltype,group2,sep="_")
DARs <- FindMarkers(
  object = so_atac_peak,
  ident.1 = Celltype_Group2,
  ident.2 = Celltype_Group1,
  test.use = 'LR',
  latent.vars = c('nCount_peaks','nFeature_peaks'),
  min.pct = 0.01,
  logfc.threshold = 0.0
)
DARs$celltype=celltype
DARs$group1=group2
DARs$group2=group1
DARs$compare_group =compare_groups[j]
DARs$gene=rownames(DARs)
DARs$pct.all = DARs$pct.1 + DARs$pct.2
DARs_file=file.path(output_path,paste0("6_ALS_ATAC_peaks_DARs_subcelltypes_big_max_Peaks_",celltype,"_",group2,"_vs_",group1,"240920_seurat_binarize_RunTFIDF_LR_var_nCount_nFeatures_peaks.tsv"))
write.table(DARs, file = DARs_file, quote = FALSE, sep = "\t", col.names = NA)
cat(celltype,"is done")
}
}





