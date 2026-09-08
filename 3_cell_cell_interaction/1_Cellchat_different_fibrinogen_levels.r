# cell-cell communication analysis

library(svglite)
library(Seurat)
library(harmony)
library(cowplot)
library(ggplot2)
library(svglite)
library(Seurat)
library(harmony)
library(cowplot)
library(ggplot2)
library(svglite)
library(Seurat)
library(SeuratDisk)
set.seed(123)
library(DoubletFinder)
library(tidyverse)
library(Matrix)
library(matrixStats)
library(tidyverse)
library(harmony)
library(cowplot)
library(ggplot2)
library(clustree)
library(Seurat)
library(paletteer)
#Library(harmony)
library(ggrepel)
library(tidydr)  #add mini-axis
library(cowplot)
library(ggpubr)
library(scCustomize)
library(viridis)
library(patchwork)  # or cowplot depending on what FeaturePlot_scCustom returns
library(Seurat)
library(edgeR)
library(dplyr)
library(tidyverse)
library(corrplot)
library(ggplot2)
library(ggpubr)


file_path <- "xxx"

PATH_O <- file.path(file_path, "ALS_review_comments_output11052025/3_Integration_CCI_explore_01052026_results")
PATH_O_data <- file.path(PATH_O,"Data")
PATH_O_data_seq <- file.path(PATH_O_data,"data_seq")
PATH_O_fig <- file.path(PATH_O,"Figures")
PATH_O_fig_polish <- file.path(PATH_O_fig,"Figure_polished")
PATH_O_fig_polish_data <- file.path(PATH_O_fig_polish,"Data")



subcelltypes_merged=c("ODC_1",   "ODC_2",   "ODC_3","ODC_4","ODC_5","Astro_1", "Astro_2","Astro_3","OPC",
"MG_1","MG_2","ExN","InN","Endo","Pericytes","Fibroblast")
pal_subcelltypes_merged=c("#667FE1", "#ACD4EC", "#6F99AD","#64AAD2", "#8AB6E9FF","#FE9586" ,"#D87B8B","#D898B9","#EBCC78", "#5DA59E" ,"#56BC9B","#9e6f8a", "#A3939D" ,  "#8C7F5F", "#C17E73", "#476D87")   ##E95C59",

subcelltypes_merged_need=c("ODC_1",   "ODC_2",   "ODC_3","ODC_4","ODC_5","Astro_1", "Astro_2","Astro_3","OPC",
"MG_1","MG_2","ExN","InN","Endo")  #peri adn fib, not shared by multiome and snRNA-seq, so ignore
pal_subcelltypes_merged_need=c("#667FE1", "#ACD4EC", "#6F99AD","#64AAD2", "#8AB6E9FF","#FE9586" ,"#D87B8B","#D898B9","#EBCC78", "#5DA59E" ,"#56BC9B","#9e6f8a", "#A3939D" ,  "#8C7F5F")   ##E95C59",

SubCelltype_merge_scANVI_Celltype=c("ODC","Astro","OPC","MG","Neuron","Vascular")

pal_SubCelltype_merge_scANVI_Celltype=c("#8AB6E9FF","#D87B8B","#EBCC78", "#5DA59E" ,"#A48BCA", "#8C7F5F")
names(pal_SubCelltype_merge_scANVI_Celltype)=SubCelltype_merge_scANVI_Celltype

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)

groups= c("LOW", "HIGH")
compare_group1=c("LOW")
compare_group2=c("HIGH")
num_groups=length(groups)
num_com_groups=length(compare_group2)

pal_groups=c("#A7C9DF","#4880B8")
names(pal_groups)=groups

Integrated_Celltypes=c("ODC","Astro","OPC","MG","ExN","InN","Endo","Fib_Peri_VLMC")
Integrated_Celltypes_need=c("ODC","Astro","OPC","MG","ExN","InN","Endo")
num_celltypes=length(Integrated_Celltypes_need)





so=readRDS(file =file.path(PATH_O_data_seq,"0_ALS_RNA_filter_gene_subcelltype_pct5_01052026.rds"))
so= subset(so,SubCelltype %in% subcelltypes_merged_need)


Idents(so)=so$SubCelltype
so$SubCelltype=as.character(so$SubCelltype)

CellChatDB <- CellChatDB.human
library(ggplot2)

Idents(so)=so$SubCelltype_merge_scANVI
so$Celltype=as.character(so$SubCelltype_merge_scANVI_Celltype)

for (group in groups){
seurat_object=subset(so,subset=Fibrinogen==group)

data.input <- seurat_object[["RNA"]]@data # normalized data matrix
Idents(seurat_object)=seurat_object$Celltype
labels <- Idents(seurat_object)
meta <- data.frame(labels = labels, row.names = names(labels))

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use


# library(future)
# plan(strategy = 'multiprocess', workers = 5)  #no use on HPC
#subset for time-saving
cellchat <- subsetData(cellchat)  #SAVE TIME   #this step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use=TRUE,  trim = 0.1, population.size = TRUE)  #CONSIDER CELL POPULATION effect

cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
data=df.net
write.table(
    data,
    file = file.path(PATH_O_data,"Cellchatv2.2.0_diff_group",
                     paste0("1_ALS_integrated_RNA_cellchatV2_celltype_01052026_", group,
                           ".tsv")),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
cellchat <- computeCommunProbPathway(cellchat)


cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
 cellchat_data_file = file.path(PATH_O_data,"Cellchatv2.2.0_diff_group",
                     paste0("1_ALS_integrated_RNA_cellchatV2_celltype_01052026_", group,
                           ".rds"))
saveRDS(cellchat,cellchat_data_file)
cat("Group",group, "cell-cell interaction based on celltype is done","\n''" )
}



