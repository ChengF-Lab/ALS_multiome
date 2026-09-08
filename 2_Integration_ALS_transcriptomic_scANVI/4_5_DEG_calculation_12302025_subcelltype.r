#DEGs based on donor-aware edgeR method |

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
library(ggrepel)
library(tidydr)
library(cowplot)
library(ggpubr)
library(scCustomize)
library(viridis)
library(patchwork)  #
library(Seurat)
library(edgeR)
library(dplyr)

file_path <- "xxx"
PATH_O <- file.path(file_path, "ALS_review_comments_output11052025/1_Cellfraction_DEGs_enrichments_out12302025")
PATH_O_data <- file.path(PATH_O,"Data")
PATH_O_data_seq <- file.path(PATH_O_data,"data_seq")
PATH_O_fig <- file.path(PATH_O,"Figures")
PATH_O_fig_polish <- file.path(PATH_O_fig,"Figure_polished")
PATH_O_fig_polish_data <- file.path(PATH_O_fig_polish,"Data")


allsamples = c("ALS1","ALS9","ALS11","ALS17", "ALS7","ALS12",      "ALS3","ALS4","ALS8","ALS16","ALS10","ALS14") #low; high
celltypes=c("ODC","Astro","OPC","MG","Neuron","Endo")
subcelltypes=c("ODC_1",   "ODC_2",   "ODC_3","ODC_4","ODC_5","Astro_1", "Astro_2","Astro_3","OPC",
"MG_1","MG_2","Neuron","Endo")
num_celltypes=length(celltypes)
num_subcelltypes=length(subcelltypes)

groups= c("LOW_Fibrinigen", "HIGH_Fibrinigen")
compare_group1=c("LOW_Fibrinigen")
compare_group2=c("HIGH_Fibrinigen")
num_groups=length(groups)
num_com_groups=length(compare_group2)
cluster_inte="Clusters_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined_reso0.2"
umap_inte="UMAP_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined"
reduction_inte="LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined_reso0.2"
cluster_RNA="Clusters_RNA_0.2_Har_reso0.2"
umap_RNA="UMAP_LSI_RNA_0.2_Har"
cluster_atac="Clusters_ATAC_0.2_Har_reso0.2"
umap_ATAC="UMAP_LSI_ATAC_0.2_Har"








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

###only do label transfer for ODC, astrocytes, microglia
path_merged_data <- file.path("xxx/0_Output_Integration_snRNAseq_together12122025/Data/scvi_results12252025/4_5_ALS_multiomeRNA_merge_extraRNA_scANVI_seurat_1230205.rds") # output file dictonary
so = readRDS(file.path(path_merged_data))
so$Age = as.numeric(so$Age)
so$Age_round <- round(so$Age)
so$PMI=as.numeric(so$PMI)
so$data_source <- factor(so$data_source)
so$ALS_group_sf <- factor(so$ALS_group_sf)
so$Sex <- factor(so$Sex)





####################################################DEGS(psedobulk: sub-celltype)##########################################
########################################################################################################################

DEGs_data_all=data.frame()
for (i in 1:length(subcelltypes_merged)){
celltype_use <- subcelltypes_merged[i]
cells_use <- WhichCells(
  so,
  expression = SubCelltype_merge_scANVI == celltype_use
)
so_ct <- subset(so, cells = cells_use)
# raw counts
counts <- as.matrix(GetAssayData(so_ct, slot = "counts"))
donor  <- so_ct$Indivi_ID
pb_counts <- rowsum(
  t(counts),
  group = donor
)
pb_counts <- t(pb_counts)

#psedobulk data
meta_pb <- so_ct@meta.data %>%
  dplyr::select(
    Indivi_ID,
    Fibrinogen,
    Age,
    Sex,
    PMI,
    ALS_group_sf,
    data_source
  ) %>%
  dplyr::distinct()

meta_pb$Age=as.numeric(meta_pb$Age)
meta_pb$PMI=as.numeric(meta_pb$PMI)
meta_pb$data_source <- factor(meta_pb$data_source)
meta_pb$ALS_group_sf <- factor(meta_pb$ALS_group_sf)
meta_pb$Sex <- factor(meta_pb$Sex)
##keep donor order
meta_pb <- meta_pb[match(colnames(pb_counts), meta_pb$Indivi_ID), ]
rownames(meta_pb) <- meta_pb$Indivi_ID
meta_pb$Age_round <- round(meta_pb$Age)

library(edgeR)
dge <- DGEList(
  counts = pb_counts,
  samples = meta_pb
)
keep <- filterByExpr(
  dge,
  group = meta_pb$Fibrinogen
)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
meta_pb$Sex <- factor(meta_pb$Sex)
  meta_pb$Fibrinogen <- factor(meta_pb$Fibrinogen, levels = c("LOW", "HIGH"))
  meta_pb$ALS_group_sf <- factor(meta_pb$ALS_group_sf, levels = c("fALS", "sALS"))
design <- model.matrix(
  ~ Fibrinogen + Age_round + Sex + PMI + ALS_group_sf+data_source,
  data = meta_pb
)


dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = "FibrinogenHIGH")

res_edgeR <- topTags(qlf, n = Inf)$table
res_edgeR$gene <- rownames(res_edgeR)
res_edgeR$Celltype <- celltype_use

DEGS_file=file.path(PATH_O_data,"DEGs_psedobulk_edgeR",paste0("1_integrated_ALS_RNA_DEGs_scANVISubcelltypes_",celltype_use,"_HIGH_vs_LOW_psedobulk_edgeR_AgeSexPMI_ALSgroup_Datasource_12302025.tsv"))
write.table(res_edgeR, file = DEGS_file, quote = FALSE, sep = "\t", col.names = NA)
if (dim(DEGs_data_all)[1]==0){DEGs_data_all=res_edgeR} else{DEGs_data_all=rbind(DEGs_data_all,res_edgeR)}
cat(celltype_use,"is done")
}


DEGS_all_file=file.path(PATH_O_data,"DEGs_psedobulk_edgeR",paste0("1_integrated_ALS_RNA_DEGs_scANVISubcelltypes_ALL_HIGH_vs_LOW_psedobulk_edgeR_AgeSexPMI_ALSgroup_Datasource_12302025.tsv"))
write.table(DEGs_data_all, file = DEGS_all_file, quote = FALSE, sep = "\t", col.names = NA)





####Test interaction model (Does fibrinogen effect differ by ALS type?)
DEGs_INTERACTION_ALL <- data.frame() #(HIGH-LOW in sALS)- (HIGH-LOW in fALS)
DEGs_SIMPLEEFFECT_ALL <- data.frame()  #sALS: HIGH-LOW ; fALS: HIGH - LOW
for (i in 1:length(subcelltypes_merged)){
  celltype_use <- subcelltypes_merged[i]
  cat("Processing:", celltype_use, "\n")
  cells_use <- WhichCells(
    so,
    expression = SubCelltype_merge_scANVI == celltype_use
  )
  so_ct <- subset(so, cells = cells_use)

  # raw counts
  counts <- as.matrix(GetAssayData(so_ct, slot = "counts"))
  donor  <- so_ct$Indivi_ID
  pb_counts <- rowsum(t(counts), group = donor)
  pb_counts <- t(pb_counts)

  # pseudobulk metadata
  meta_pb <- so_ct@meta.data %>%
    dplyr::select(
      Indivi_ID,
      Fibrinogen,
      Age,
      Sex,
      PMI,
      ALS_group_sf,
      data_source
    ) %>%
    dplyr::distinct()

  # 数据类型转换
  meta_pb$Age <- as.numeric(meta_pb$Age)
  meta_pb$PMI <- as.numeric(meta_pb$PMI)

  meta_pb$Fibrinogen <- factor(meta_pb$Fibrinogen, levels = c("LOW", "HIGH"))
  meta_pb$ALS_group_sf <- factor(meta_pb$ALS_group_sf, levels = c("fALS", "sALS"))
  meta_pb$Sex <- factor(meta_pb$Sex)
  meta_pb$data_source <- factor(meta_pb$data_source)

  meta_pb <- meta_pb[match(colnames(pb_counts), meta_pb$Indivi_ID), ]
  rownames(meta_pb) <- meta_pb$Indivi_ID
  meta_pb$Age_round <- round(meta_pb$Age)

  library(edgeR)
  dge <- DGEList(counts = pb_counts, samples = meta_pb)

  keep <- filterByExpr(dge, group = meta_pb$Fibrinogen)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)

  design <- model.matrix(
    ~ Fibrinogen * ALS_group_sf + Age_round + Sex + PMI + data_source,
    data = meta_pb
  )
colnames(design) <- make.names(colnames(design))
print(colnames(design))


dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

  qlf_interaction <- glmQLFTest(fit, coef = "FibrinogenHIGH.ALS_group_sfsALS")
  res_interaction <- topTags(qlf_interaction, n = Inf)$table
  res_interaction$gene <- rownames(res_interaction)
  res_interaction$Celltype <- celltype_use
  res_interaction$Subgroup <- "fALS and sALS"
  res_interaction$Effect   <- "Interaction_sALS_vs_fALS"
  DEGs_INTERACTION_ALL <- rbind(DEGs_INTERACTION_ALL, res_interaction)

  contrast_fALS <- makeContrasts(
    Fibrinogen_in_fALS = FibrinogenHIGH,
    levels = design
  )

  qlf_fALS <- glmQLFTest(fit, contrast = contrast_fALS)
  res_fALS <- topTags(qlf_fALS, n = Inf)$table
  res_fALS$gene     <- rownames(res_fALS)
  res_fALS$Celltype <- celltype_use
  res_fALS$Subgroup <- "fALS"
  res_fALS$Effect   <- "Fibrinogen_HIGH_vs_LOW in fALS"

  ## ============================================================
  ## extract Simple effect in sALS (HIGH fibrinogen in fALS  + interaction)
  ## ============================================================
  contrast_sALS <- makeContrasts(
    Fibrinogen_in_sALS =
      FibrinogenHIGH + FibrinogenHIGH.ALS_group_sfsALS,
    levels = design
  )

  qlf_sALS <- glmQLFTest(fit, contrast = contrast_sALS)
  res_sALS <- topTags(qlf_sALS, n = Inf)$table
  res_sALS$gene     <- rownames(res_sALS)
  res_sALS$Celltype <- celltype_use
  res_sALS$Subgroup <- "sALS"
  res_sALS$Effect   <- "Fibrinogen_HIGH_vs_LOW in sALS"
  DEGs_SIMPLEEFFECT_ALL <- rbind(
    DEGs_SIMPLEEFFECT_ALL,
    res_fALS,
    res_sALS
  )
}





write.table(
  DEGs_INTERACTION_ALL,
  file = file.path(PATH_O_data, "DEGs_psedobulk_edgeR",
                  "1_integrated_ALS_RNA_DEGs_scANVISubcelltypes_ALL_HIGH_vs_LOW_psedobulk_edgeR_INTERACTION_FibrinogenALSgroup_vars_AgeSexPMIdatasource12302025_INTERACTIONEFFECT.tsv"),
  quote = FALSE, sep = "\t", col.names = NA
)

write.table(
  DEGs_SIMPLEEFFECT_ALL,
  file = file.path(PATH_O_data, "DEGs_psedobulk_edgeR",
                  "1_integrated_ALS_RNA_DEGs_scANVISubcelltypes_ALL_HIGH_vs_LOW_psedobulk_edgeR_INTERACTION_FibrinogenALSgroup_vars_AgeSexPMIdatasource12302025_SIMPLEEFFECT.tsv"),
  quote = FALSE, sep = "\t", col.names = NA
)

