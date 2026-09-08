# batch correction and cell clustering and annotation |

set.seed(123)
library(Seurat)
library(DoubletFinder)
library(harmony)
library(tidyverse)
library(Matrix)
library(matrixStats)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(clustree)
library(paletteer)
library(tidydr)




file_path <- "xxxx"  #give general dir of input file, which further includes GSE file of different sample (related details of dir, age, sex ,etc. are recorded in "sample_info.tsv"  )


PATH_O <- file.path(file_path, "ALS_review_comments_output11052025/0_Output_Integration_snRNAseq_together12122025") # output file dictonary
PATH_O_data <- file.path(PATH_O,"Data")
PATH_O_data_seq <- file.path(PATH_O_data,"data_seq")
PATH_O_fig <- file.path(PATH_O,"Figures")
PATH_O_fig_polish <- file.path(PATH_O_fig,"Figure_polished")


Indivi_IDs=c("ALS07","ALS10","ALS12","ALS14")
Sample_IDs=c("ALS-7","ALS-10","ALS-12","ALS-14")
Fibrinigen_list=c("HIGH","HIGH","LOW","LOW")
num_samp=length(Indivi_IDs)


so<-readRDS(file = file.path(PATH_O_data_seq,"xxx.rds"))


so <- subset(so, subset =
doubFindadj_res =="Singlet")

Idents(so)=so$Indivi_ID
feas = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb")
vlnplot=VlnPlot(so, pt.size=0,features = feas, ncol = 3)+NoLegend()
ggsave(vlnplot,filename = file.path(PATH_O_fig,"1_0_ALS_huamn_AfterQC_doublets_Vln12252025.png"),width =14,height = 8,dpi=100)
#scatter plot
plot1 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(plot1+plot2,filename = file.path(PATH_O_fig,"1_0_ALS_extra_AfterQC_doublets_FeaScatter12252025.png"),width =10,height =5,dpi=100)

p4=DimPlot(so,label=FALSE,group.by="Indivi_ID",reduction="UMAP_doublets",raster= FALSE)
ggsave(p4,filename = file.path(PATH_O_fig,"1_0_ALS_extra_Umap_after_doublets012252025.png"),width =7, height = 5,dpi=100)








###################################################################################################################################
###########################batch correction by harmony  based on different Indivi_ID and perform  clustering ########################
#####################################################################################################################################
dim.usage=30
so@meta.data$Indivi_ID<-as.factor(so@meta.data$Indivi_ID)  # if not, error like: Error in harmonyObj$init_cluster_cpp(0) :
so <- RunPCA(so,reduction.name = 'PCA', reduction.key = 'PCA_')
so <- RunHarmony(so, reduction="PCA",c("Indivi_ID"),reduction.save = "harmony_Indivi_ID")
so <- RunUMAP(so,  dims = 1:dim.usage,reduction = "harmony_Indivi_ID",reduction.name = 'UMAP_harmony_Indivi_ID')

saveRDS(so, file = file.path(PATH_O_data_seq,"xxx.rds"),compress=FALSE)



p2=DimPlot(so,reduction="UMAP_doublets",group.by="Indivi_ID",raster=FALSE)
p2_3=DimPlot(so,reduction="UMAP_harmony_Indivi_ID",group.by="Indivi_ID",label=F,raster=FALSE)
p=plot_grid(p2+p2_3, ncol = 1)
ggsave(p,filename = file.path(PATH_O_fig,"1_1_ALS_extra_human_data_Umap_cluster_Group_diff_batch_harmony12252025.png"),width =10,height = 5,dpi=100)


Idents(so)=so$Indivi_ID
so <- FindNeighbors(so,reduction="harmony_Indivi_ID",graph.name="harmony_Indivi_ID_knn",k.param = 30,  dims = 1:dim.usage)
so <- FindClusters(so, graph.name="harmony_Indivi_ID_knn",resolution = c(0.05,0.1,0.2,0.4))
saveRDS(so, file = file.path(PATH_O_data_seq,"1_1_ALS_extra_human_data_harmony_12252025.rds"),compress=FALSE)

p3_1=DimPlot(so,reduction="UMAP_harmony_Indivi_ID",group.by="harmony_Indivi_ID_knn_res.0.05",label=TRUE,raster=FALSE)
p3_2=DimPlot(so,reduction="UMAP_harmony_Indivi_ID",group.by="harmony_Indivi_ID_knn_res.0.1",label=TRUE,raster=FALSE)
p3_3=DimPlot(so,reduction="UMAP_harmony_Indivi_ID",group.by="harmony_Indivi_ID_knn_res.0.2",label=TRUE,raster=FALSE)
p3_4=DimPlot(so,reduction="UMAP_harmony_Indivi_ID",group.by="harmony_Indivi_ID_knn_res.0.4",label=TRUE,raster=FALSE)

p=plot_grid(p3_1,p3_2,p3_3,p3_4, ncol = 2)
ggsave(p,filename = file.path(PATH_O_fig,"1_1_ALS_extra_data_Umap_cluster_Group_diff_batch_harmony_diff_reso12252025.png"),width =12,height = 14,dpi=100)





###################################################################################################################################
###########################################celltype annotation based on the markers#############################################
####################################################################################################################################################

so<-readRDS(file = file.path(PATH_O_data_seq,"xxxx.rds"))

#choose reso: 0.1
so$KNN_cluster=so$harmony_Indivi_ID_knn_res.0.4
Idents(so)=so$KNN_cluster

## ============ Cell type identification ============
#Must give the makergene name, if not, error will be occurred; if the given maker gene is not involved in default file "scale.data", a warning will be reported.
MARKERS <- unique(c("RBFOX3","SYT1",  #Neuron
"SLC17A7","PDE1A",  #eXn
"ATP2B1","PPFIA", "ADARB2","LHX6","SST", #inn
"GAD1","GAD2","DLX6OS1","NPY", #iNN
"AQP4","GJA1", #Astro
"P2RY12","CX3CR1","C3", #MG_micropahge
"VCAN","OLIG1", #ODC
"MOBP", "MOG", #OPC
"CLDN5","FLT1","COLIA2","MGP","SPOCK2",#endo
"PDGFRB","DCN","VTN","LAMA2","ATP1A2", #PERI
"FOXJ1","CALML4","PIFO","CFAP299", #Epen
"ABCA9","FBLN1" #fibroblast
))

p0<-DotPlot(so, features =MARKERS,group.by="KNN_cluster",cluster.idents = TRUE)+RotatedAxis()
svg(file.path(PATH_O_fig, "1_2_ALS_extra_human_data_DotPlot0_markergenes_12252025.svg"),height=6,width=15)
p0
dev.off()

vascular_ependy_markers_chatgpt  <- unique(c(
  ## Endothelial cells
  "CLDN5","PECAM1","VWF","KDR","FLT1","ESAM",
  "ENG","RAMP2","PLVAP","ROBO4","CD34",
  ## Pericytes / mural cells
  "PDGFRB","RGS5","CSPG4","MCAM","ABCC9",
  "KCNJ8","NOTCH3","DES",
  ## Vascular smooth muscle cells
  "ACTA2","TAGLN","MYH11","CNN1","CALD1",
  ## Fibroblast / vascular fibroblast-like
  "COL1A1","COL1A2","DCN","LUM","COL3A1",
  "FBLN1","PDGFRA",
  ## Ependymal cells (gold standard)
  "FOXJ1","TPPP3","PIFO","CALML4",
  "DNAH5","DNAH9","IFT88","CFAP43",
  ## Ependymal epithelial / CSF interface
  "EPCAM","KRT8","KRT18","SLC12A2","AQP4"
))
p0<-DotPlot(so, features =vascular_ependy_markers_chatgpt,group.by="KNN_cluster",cluster.idents = TRUE)+RotatedAxis()
svg(file.path(PATH_O_fig, "1_2_ALS_extra_human_data_DotPlot0_markergenes_vascular_12252025.svg"),height=6,width=15)
p0
dev.off()


so$KNN_cluster=so$harmony_Indivi_ID_knn_res.0.4
Idents(so)=so$KNN_cluster
so <- RenameIdents(
  so,
`0` = "ODC",
`1` = "Astro", #
`2` = "ExN", #
 `3` = "OPC",#
`4` = "MG",
`5` = "ODC",  #
`6` = "ExN",  #HARD SAY MG OR NEURON
`7` = "ODC",  #
`8` = "InN_VIP_RELN",
`9` = "Pericytes",  #
`10` = "Endo", #
`11` = "InN_LHX6", #
`12` = "InN_LHX6",  #
`13` = "InN",
`14` = "Fibroblast",  #
`15` = "InN_LHX6", #
`16` = "ExN", #
`17` = "ExN"
)

so$SubCelltype_har_0.4=Idents(so)




p0<-DimPlot(so, reduction="UMAP_harmony_Indivi_ID",group.by="SubCelltype_har_0.4",label=TRUE,raster=TRUE)+RotatedAxis()
svg(file.path(PATH_O_fig, "1_2_ALS_extra_human_data_UMAP_12252025_subcelltype.svg"),height=6,width=6)
p0
dev.off()

p0<-DimPlot(so, reduction="UMAP_harmony_Indivi_ID",group.by="Celltype_har_0.4",label=TRUE,raster=TRUE)+RotatedAxis()
svg(file.path(PATH_O_fig, "1_2_ALS_extra_human_data_UMAP_12252025.svg"),height=6,width=6)
p0
dev.off()


MARKERS_show <- unique(c("RBFOX3","SYT1",  #Neuron
"SLC17A7","PDE1A",  #eXn
"GAD1","GAD2","LHX6","SST","VIP","RELN","SST", #inn
"AQP4","GJA1", "GFAP",#Astro
"P2RY12","CX3CR1", #MG_micropahge
"VCAN","OLIG1", #ODC
"MOBP", "MOG", #OPC
"CLDN5","FLT1",#endo
"PDGFRB","MYH11","MCAM", #PERI
"DCN","COL1A1","FBLN1",#fibroblast
"ABCA9","FBLN1" #fibroblast
))
p0<-DotPlot(so, features =MARKERS_show,group.by="Celltype_har_0.4")+RotatedAxis()
svg(file.path(PATH_O_fig, "1_2_ALS_extra_human_data_DotPlot_markergenes_12252025.svg"),height=4,width=10)
p0
dev.off()

p0<-DotPlot(so, features =MARKERS_show,group.by="SubCelltype_har_0.4")+RotatedAxis()
svg(file.path(PATH_O_fig, "1_2_ALS_extra_human_data_DotPlot_markergenes_12252025_Subcelltype.svg"),height=4,width=10)
p0
dev.off()



saveRDS(so, file = file.path(PATH_O_data_seq,"xxx.rds"))# compress, less memory, bur more time to save

































































































































































































































































































































