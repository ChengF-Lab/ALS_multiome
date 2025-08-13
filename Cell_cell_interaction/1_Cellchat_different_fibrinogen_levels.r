# cell-cell communication analysis
library(Seurat)
library(SeuratDisk)
library(tidyverse)
set.seed(123)
library(CellChat)
library(svglite)





input_path="xxx"
output_path="xxx"



groups <- c("LOW_Fibrinigen", "HIGH_Fibrinigen")
compare_group1 <- c("LOW_Fibrinigen")
compare_group2 <- c("HIGH_Fibrinigen")
compare_groups <- paste0(compare_group2, "_vs_", compare_group1)
num_groups <- length(groups)
num_com_groups <- length(compare_group2)

celltypes=c("ODC","Astro","OPC","MG","Neuron","Endo")
subcelltypes=c("ODC_1",   "ODC_2",   "ODC_3","ODC_4","ODC_5","Astro_1", "Astro_2","Astro_3","OPC",
"MG_1","MG_2","Neuron","Endo")
num_subcelltypes=length(subcelltypes)




so<-readRDS(file.path(input_path,"xxx.rds"))
Idents(so)=so$SubCelltype
so$SubCelltype=as.character(so$SubCelltype)  #remove level



#######################################################filter low-expressed genes
gene_list=list()
for (celltype  in subcelltypes){
so_sub=subset(so,subset=SubCelltype==celltype)
data_matrix <- GetAssayData(so_sub, slot = "counts")
num_nuclei=dim(so_sub)[2]
binary_matrix <- data_matrix > 0
nuclei_count_per_gene <- Matrix::rowSums(binary_matrix)
nuclei_count_per_gene=nuclei_count_per_gene[nuclei_count_per_gene>num_nuclei*0.05]
gene_list[[celltype]]=names(nuclei_count_per_gene)
}
saveRDS(gene_list, file =file.path(output_path,"gene_list_keep.rds"))

all_genes_keep=unique(unlist(gene_list))
so_filter=so[all_genes_keep,]
saveRDS(so_filter, file =file.path(output_path,"xxx_for_cellchat.rds"))


#####################################################START CELLCHAT ANALYSIS
so<-readRDS(file.path(input_path,"xxx_for_cellchat.rds"))
Idents(so)=so$SubCelltype
so$SubCelltype=as.character(so$SubCelltype)  #remove level

CellChatDB <- CellChatDB.human
library(ggplot2)
p=showDatabaseCategory(CellChatDB)
svg(file.path(PATH_O_fig,"0_ALS_cellchat_human_database_category.svg"))
p
dev.off()


for (i in 1:num_groups){
cat("Cellchat analysis for different groups")
group=groups[i]
so_subset=subset(so,subset=Fibrinigen==group)
data <- GetAssayData(object = so_subset, layer = 'data')
meta <- so_subset@meta.data

######################################################creat cellchat object##################################################
seurat_annotations="SubCelltype"  #give the colname which give celltypes
cellchat <- createCellChat(object = data,
                           meta = meta,
                           group.by = seurat_annotations)

#################################################data preparation#################
cellchat@DB = CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

##########################################infer cell-cell interaction################
cellchat <- computeCommunProb(cellchat,raw.use=FALSE,population.size = TRUE)   # raw.use: whether use the raw data (i.e., `object@data.signaling`) or
cellchat <- filterCommunication(cellchat, min.cells = 10) #Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))#
cellchat_data_file=file.path(output_path,paste0("0_ALS_cellchat_",group,"_allhumandatabase250220_diff_group_subcelltype_based.rds"))
saveRDS(cellchat,cellchat_data_file)
cat("Group",group, "cell-cell interaction based on subcelltype is done","\n''" )
}


#COMPARE INTERACTIONS BETWEEN DIFFERENT GROUPS
group_colors_for_sort_fill = c("#A7CADF","#4880B8")
gg1 <- compareInteractions(cellchat, show.legend = T,title.name = paste0("ALS_group_compare_count"),group = c(1,2),color.use=group_colors_for_sort_fill)+theme(legend.position="none",axis.text.x=element_text(angle=90))+labs(title="")
gg2 <- compareInteractions(cellchat, show.legend = T,title.name = paste0("ALS_group_compare_weight"), group = c(1,2),measure = "weight",color.use=group_colors_for_sort_fill)+theme(legend.position="none",axis.text.x=element_text(angle=90))+labs(title="")
ggsave(gg1+gg2,filename=file.path(output_path,paste0("2_2_ALS_group_comparison_barplot_allhumandatabase250220.svg")),bg="white",height=3.8,width=5)

