##merged data for downstream analysis

library(Seurat)
library(SeuratDisk)
set.seed(123)
library(ggplot2)
library(cowplot)
library(scCustomize)
library(viridis)
library(patchwork)
library(Seurat)
library(edgeR)
library(dplyr)
library(tidyverse)
library(corrplot)
library(ggplot2)
library(ggpubr)

file_path <- "xxx"
PATH_O <- file.path(file_path, "ALS_review_comments_output11052025/0_Output_Integration_snRNAseq_together12122025")
PATH_O_data <- file.path(PATH_O,"Data")
PATH_O_data_seq <- file.path(PATH_O_data,"data_seq")
PATH_O_scvi <- file.path(PATH_O_data,"scvi_results12252025")
PATH_O_fig <- file.path(PATH_O,"Figures")
PATH_O_fig_polish <- file.path(PATH_O_fig,"Figure_polished")

so_merged = readRDS(file.path(PATH_O_data_seq,"2_2_ALS_multiomeRNA_merge_extraRNA_with_clustering_1225205.rds"))
umap <- read.csv(
  file.path(PATH_O_scvi,
    "4_4_ALS_multiomeRNA_merge_extraRNA_all_nuclei_scanvi_umap_12302025.csv"),
  sep = "\t",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

celltype <- read.csv(
  file.path(PATH_O_scvi,
    "4_4_ALS_multiomeRNA_merge_extraRNA_all_nuclei_scanvi_celltype_12302025.csv"),
  sep = "\t",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)


identical(colnames(so_merged)  ,rownames(celltype))
so_merged$C_scANVI <- celltype$C_scANVI
so_merged$SubCelltype_merge_scANVI <- celltype$SubCelltype_merge_scANVI

identical(colnames(so_merged)  ,rownames(umap))
umap_use <- as.matrix(umap)
umap_scanvi <- CreateDimReducObject(
  embeddings = umap_use,
  key = "SCANVIUMAP_",
  assay = DefaultAssay(so_merged)
)
so_merged[["UMAP_scANVI_batch_IndiviID"]] <- umap_scanvi


###add metedata information from Tony
so_merged$Fibrinogen=ifelse(so_merged$Indivi_ID %in% c("ALS3","ALS4","ALS8","ALS16","ALS10","ALS14"),"HIGH","LOW")
so_merged$PMI=ifelse(so_merged$Indivi_ID %in% c("ALS07"), 6.8, so_merged$PMI)
so_merged$PMI=ifelse(so_merged$Indivi_ID %in% c("ALS10"), 6.5, so_merged$PMI)
so_merged$PMI=ifelse(so_merged$Indivi_ID %in% c("ALS12"), 12.2, so_merged$PMI)
so_merged$PMI=ifelse(so_merged$Indivi_ID %in% c("ALS14"), 6.6, so_merged$PMI)
so_merged$Indivi_ID=ifelse(so_merged$Indivi_ID %in% c("ALS07"), "ALS7", so_merged$Indivi_ID)
so_merged$SubCelltype_merge_scANVI_Celltype =ifelse (so_merged$SubCelltype_merge_scANVI %in% c("ODC_1",   "ODC_2",   "ODC_3","ODC_4","ODC_5"),"ODC",so_merged$SubCelltype_merge_scANVI)
so_merged$SubCelltype_merge_scANVI_Celltype =ifelse (so_merged$SubCelltype_merge_scANVI %in% c("Astro_1", "Astro_2","Astro_3"),"Astro",so_merged$SubCelltype_merge_scANVI_Celltype)
so_merged$SubCelltype_merge_scANVI_Celltype =ifelse (so_merged$SubCelltype_merge_scANVI %in% c("MG_1","MG_2"),"MG",so_merged$SubCelltype_merge_scANVI_Celltype)
so_merged$SubCelltype_merge_scANVI_Celltype =ifelse (so_merged$SubCelltype_merge_scANVI %in% c("Endo","Pericytes","Fibroblast"),"Vascular",so_merged$SubCelltype_merge_scANVI_Celltype)
so_merged$SubCelltype_merge_scANVI_Celltype =ifelse (so_merged$SubCelltype_merge_scANVI %in% c("ExN","InN"),"Neuron",so_merged$SubCelltype_merge_scANVI_Celltype)

saveRDS(so_merged,file.path(PATH_O_scvi,"4_5_ALS_multiomeRNA_merge_extraRNA_scANVI_seurat_1230205.rds"))



subcelltypes_merged=c("ODC_1",   "ODC_2",   "ODC_3","ODC_4","ODC_5","Astro_1", "Astro_2","Astro_3","OPC",
,"ExN","InN","Endo","Pericytes","Fibroblast")

p2_4=DimPlot(so_merged,reduction="UMAP_harmony_IndiviID_datasource",group.by="SubCelltype_merge_scANVI", size=3,label=F,raster=TRUE)+
  labs(x = "UMAP_1", y = "UMAP_2", color = "donors")
p2_5=DimPlot(so_merged,reduction="UMAP_scANVI_batch_IndiviID",group.by="SubCelltype_merge_scANVI", size=3,label=F,raster=TRUE)+
  labs(x = "UMAP_1", y = "UMAP_2", color = "donors")
p=plot_grid(p2_4,p2_5, ncol = 2)
ggsave(p,filename = file.path(PATH_O_fig,"ALS_Multiome_extraRNA_integration_label_transfer_scVI","4_5_ALS_margedRNA_umap_groupby_SubCelltype_merge_scANVI_12252025.svg"),width =10,height = 4.5,dpi=100)

so_merged = readRDS(file.path(PATH_O_scvi,"4_5_ALS_multiomeRNA_merge_extraRNA_scANVI_seurat_1230205.rds"))
pal <- viridis(n = 10, option = "C")
markers <- c("MOBP",
"VCAN",
"AQP4",
"SLC17A7","GAD2",
"CSF1R",
"CLDN5",
"MYH11",
"FBLN1")

pal <- viridis(n = 10, option = "C")
plots <- lapply(markers, function(gene) {
  FeaturePlot_scCustom(
    seurat_object = so_merged,
    features = gene,  # Use each individual gene
    colors_use = pal,
    reduction = "UMAP_scANVI_batch_IndiviID",
    raster = TRUE,
    pt.size = 3)
}+labs(x="",y=""))

p_combined <- wrap_plots(plots,nrow=2)
ggsave(p_combined,filename = file.path(PATH_O_fig,"ALS_Multiome_extraRNA_integration_label_transfer_scVI","4_5_ALS_margedRNA_celltypemarker_featureplot_12252025.svg"),width =12,height = 3.8,dpi=300)


MARKERS_heatmap_show <- c(
  # Oligodendrocytes
  "PLP1", "MBP", "ST18", "MOBP",
  # Astrocytes
  "AQP4", "GFAP",
  # OPCs
  "CSPG4", "OLIG1",
  # Microglia
  "CX3CR1", "CSF1R", "P2RY12",
  # Neurons
  "SYT1", "SLC17A7", "GAD1",
  # Endothelial cells
  "CLDN5", "FLT1",
  # Vascular smooth muscle cells
  "MYH11",
  # Fibroblasts
  "DCN"
)

Idents(so_merged)=so_merged$SubCelltype_merge_scANVI

marker_dot=DotPlot(so_merged,features=MARKERS_heatmap_show,group.by="SubCelltype_merge_scANVI",)$data
marker_dot$features.plot=factor(marker_dot$features.plot,levels=(MARKERS_heatmap_show))
marker_dot$id=factor(marker_dot$id,levels=rev(subcelltypes_merged))
p1 <- ggplot(marker_dot, aes(x = features.plot, y = id, size = pct.exp)) +
  geom_point(aes(fill = avg.exp.scaled), shape = 21, color = "grey50", lwd = 2) +  # Use shape 21 for outlined circles
  scale_fill_gradient2(
    low = "#B7D4E9",
    mid = "#F5F9FC",
    high = "#D21E1C",
    midpoint = 0,
    limits = c(-2.5, 2.5)  #avg.exp.scaled
  ) +  # Set mid color
  theme_bw() +
  labs(x = NULL, y = NULL, size = "pct.exp", fill = "avg.exp") +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.text.y = element_text(size = 8),
    legend.position = "right",
    legend.title = element_text(size = 8),
    axis.text.y = element_text(size = 8),  # Set y-axis text size
    axis.title.y = element_text(size = 9),
    axis.title.x = element_text(size = 8),
    axis.ticks = element_blank(),
    axis.text.x = element_text(color = "black", angle = 90, vjust = 0.65, size = 8),  # Set x-axis text size
    panel.grid = element_blank(),  # Remove grid lines
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(0.3, "cm"),
  )+
  scale_size(range = c(0.5, 4))
#   + # Adjust the size range of the points
#   coord_flip()


ggsave(p1,filename = file.path(PATH_O_fig,"ALS_Multiome_extraRNA_integration_label_transfer_scVI","4_5_ALS_margedRNA_celltypemarker_dotplot_12252025.svg"),width =5,height =3,dpi=100,bg="white")

