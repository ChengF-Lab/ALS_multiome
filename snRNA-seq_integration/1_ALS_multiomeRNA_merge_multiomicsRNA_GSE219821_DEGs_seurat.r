##merge RNA data with public data within motor cortex
library(Seurat)
library(harmony)
library(cowplot)
library(ggplot2)
library(svglite)



muliome_ALS_RNA_file="xxx.rds"
muliomics_ALS_RNA_file="xxx.rds"
out_path="xxx"


celltypes=c("ODC","Astro","OPC","MG","Neuron","Endo")
num_celltypes=length(celltypes)

groups= c("Control_Unknown","ALS_LOW_Fibrinigen", "ALS_HIGH_Fibrinigen")
compare_group1=c("Control_Unknown","Control_Unknown","ALS_LOW_Fibrinigen")
compare_group2=c("ALS_LOW_Fibrinigen","ALS_HIGH_Fibrinigen","ALS_HIGH_Fibrinigen")
compare_groups=paste0(compare_group2,"_vs_",compare_group1)
num_groups=length(groups)
num_com_groups=length(compare_group2)




#####MOTOR CORTEX
#######5. FIB ALS  vs, C9-ALS
#######6. FIB ALS vs, Control
so=readRDS(file.path(out_path_data_seq,"0_2_ALS_multiomeRNA_merge_multiomicsRNA_GSE219821_afterQC_241120.rds"))
so_subset=subset(so,subset=Region != "medial frontal cortex")
so_subset$ALS_type_Celltype=paste(so_subset$Celltype,so_subset$ALS_type,sep="_")   ###---> totally have ALS_HIGH_Fibrinigen ALS_LOW_Fibrinigen ALS_Unknown Control_Unknown

groups= c("Control","C9-ALS", "Fibrinigen_ALS")
compare_group1=c("Control","C9-ALS")
compare_group2=c("Fibrinigen_ALS","Fibrinigen_ALS")
compare_groups=paste0(compare_group2,"_vs_",compare_group1)
num_groups=length(groups)
num_com_groups=length(compare_group2)
vars_need=c(c("nCount_RNA","nFeature_RNA","data_source"),c("nCount_RNA","nFeature_RNA","data_source"))

Idents(so_subset)=so_subset$ALS_type_Celltype
for (i in 1:num_celltypes){
    for (j in 1:num_com_groups){
celltype=celltypes[i]
group1=compare_group1[j]
group2=compare_group2[j]
Celltype_Group1=paste(celltype, group1,sep="_")
Celltype_Group2=paste(celltype, group2, sep="_")
var_need=vars_need[[j]]
DEGs <- FindMarkers(so_subset,ident.1 =Celltype_Group2, ident.2 = Celltype_Group1,logfc.threshold=0,test.use="MAST",latent.vars=var_need) #usually ident. vs. ident.2 , so we need to give right ident.1 and ident.2)  #usually ident. vs. ident.2 , so we need to give right ident.1 and ident.2
DEGs$celltype=celltype
DEGs$group1=group2
DEGs$group2=group1
DEGs$compare_group =compare_groups[j]
DEGs$gene=rownames(DEGs)
DEGs$pct.all = DEGs$pct.1 + DEGs$pct.2
DEGs$region = "MotorCortex"
DEGS_file=file.path(out_path_data_DEGs,paste0("1_ALS_RNA_DEGs_celltypes_",celltype,"_",group2,"_vs_",group1,"_seurat_MAST_var_nCount_nFeature_datasource_250415_MotorCortex.tsv"))
write.table(DEGs, file = DEGS_file, quote = FALSE, sep = "\t", col.names = NA)
cat(celltype,"is done")
}
}

#####frontal CORTEX
#######7. C9-ALS  vs, Control
so=readRDS(file.path(out_path_data_seq,"0_2_ALS_multiomeRNA_merge_multiomicsRNA_GSE219821_afterQC_241120.rds"))
#here only check motor cortex data
so_subset=subset(so,subset=Region == "medial frontal cortex")
so_subset$ALS_type_Celltype=paste(so_subset$Celltype,so_subset$ALS_type,sep="_")   ###---> totally have ALS_HIGH_Fibrinigen ALS_LOW_Fibrinigen ALS_Unknown Control_Unknown

groups= c("Control")
compare_group1=c("Control")
compare_group2=c("C9-ALS")
compare_groups=paste0(compare_group2,"_vs_",compare_group1)
num_groups=length(groups)
num_com_groups=length(compare_group2)
vars_need=c(c("nCount_RNA","nFeature_RNA"))

Idents(so_subset)=so_subset$ALS_type_Celltype
for (i in 1:num_celltypes){
    for (j in 1:num_com_groups){
celltype=celltypes[i]
group1=compare_group1[j]
group2=compare_group2[j]
Celltype_Group1=paste(celltype, group1,sep="_")
Celltype_Group2=paste(celltype, group2, sep="_")
var_need=vars_need[[j]]
DEGs <- FindMarkers(so_subset,ident.1 =Celltype_Group2, ident.2 = Celltype_Group1,logfc.threshold=0,test.use="MAST",latent.vars=var_need) #usually ident. vs. ident.2 , so we need to give right ident.1 and ident.2)  #usually ident. vs. ident.2 , so we need to give right ident.1 and ident.2
DEGs$celltype=celltype
DEGs$group1=group2
DEGs$group2=group1
DEGs$compare_group =compare_groups[j]
DEGs$gene=rownames(DEGs)
DEGs$pct.all = DEGs$pct.1 + DEGs$pct.2
DEGs$region = "FrontalCortex"
DEGS_file=file.path(out_path_data_DEGs,paste0("1_ALS_RNA_DEGs_celltypes_",celltype,"_",group2,"_vs_",group1,"_seurat_MAST_var_nCount_nFeature_datasource_250415_FrontalCortex.tsv"))
write.table(DEGs, file = DEGS_file, quote = FALSE, sep = "\t", col.names = NA)
cat(celltype,"is done")
}
}


