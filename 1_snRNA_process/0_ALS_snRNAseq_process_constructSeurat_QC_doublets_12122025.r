#For additional snRNA-seq data, performing bioinformatics analysis, including QC, doublet removal |

#Load package
library("Seurat")
library(DoubletFinder)
library(future)
set.seed(123)
library(ggplot2)
library(Seurat)
library(DoubletFinder)
library(tidyverse)
library(Matrix)
library(matrixStats)
set.seed(123)
library(ggplot2)



file_path <- "xxx"
PATH_I_sequencing=file.path(file_path,"ALS_review_comments_analysis11052025/0_Integration_snRNAseq_together12122025/Sequencing_data_GENEWIZ_12122025/Cellranger_v9.0.1")
PATH_I_sample_info=file.path(file_path,"ALS_review_comments_analysis11052025/0_Integration_snRNAseq_together12122025/Sequencing_data_GENEWIZ_12122025/Sample_info_extra_four_snRNAseq.tsv")
PATH_O <- file.path(file_path, "ALS_review_comments_output11052025/0_Output_Integration_snRNAseq_together12122025")
PATH_O_data <- file.path(PATH_O,"Data")
PATH_O_data_seq <- file.path(PATH_O_data,"data_seq")
PATH_O_fig <- file.path(PATH_O,"Figures")
PATH_O_fig_polish <- file.path(PATH_O_fig,"Figure_polished")

####################################################################################################
##############################creating the seurat object########################################

##sample info file
sample_info=read.csv(PATH_I_sample_info,sep="\t")
num_samp=dim(sample_info)[1]


Indivi_IDs=c("ALS07","ALS10","ALS12","ALS14")
Sample_IDs=c("ALS-7","ALS-10","ALS-12","ALS-14")
Fibrinigen_list=c("HIGH","HIGH","LOW","LOW")
num_samp=length(Indivi_IDs)


so_ls <- list()
for (i in 1:num_samp){
	Sample_ID <- Sample_IDs[i]
    Indivi_ID <- Indivi_IDs[i]
    Tissue <- "ALS motor cortex"
    posi=match(Indivi_ID,sample_info$subject.ID)
    Sex <- as.character(sample_info$SEX[posi])
    Group_Fibrinigen <- Fibrinigen_list[i]
    Age_onset <- as.character(sample_info$age.at.onset[posi])
    Age_death <- as.character(sample_info$age.at.death[posi])
    mutation <- as.character(sample_info$mutation[posi])
    disease_furation_yrs = as.numeric(Age_death)-as.numeric(Age_onset)
    pmi<- "unknown"
    ALS_group=as.character(sample_info$fALS.sALS[posi])
	curr_so.data<-Read10X_h5(file.path(PATH_I_sequencing,Sample_ID, "outs/filtered_feature_bc_matrix.h5"))
    curr_so <- CreateSeuratObject(counts = curr_so.data, min.cells=3, min.features=200)  #with default condition to create seurat object
    curr_so$Sample_ID <-Sample_ID
    curr_so$Indivi_ID <- Indivi_ID
    curr_so$Tissue <- Tissue
    curr_so$Sex <- Sex
    curr_so$Group_Fibrinigen <- Group_Fibrinigen
    curr_so$Age_onset <- Age_onset
    curr_so$Age_death <-Age_death
    curr_so$disease_furation_yrs <- disease_furation_yrs
    curr_so$pmi <- pmi
    curr_so$date_post_fixation_MRI <- "unknown"
    curr_so$ALS_group <- ALS_group
    curr_so$onset_location <- "unknown"
 	curr_so[["percent.mt"]] <- PercentageFeatureSet(curr_so, pattern = "^MT-")
    curr_so[["percent.ribo"]] <- PercentageFeatureSet(curr_so, pattern = "^RP[SL]")
    curr_so[["percent.hb"]] <- PercentageFeatureSet(curr_so, pattern = "^HB[^(P)]")
    so_ls[[i]] <- curr_so
}

print(so_ls)
subsetList <- function(myList, elementNames) {
  sapply(elementNames, FUN=function(x) myList[[x]])
}
generateID <- function(ele) {
  return(paste0("ALS_snRNAseq_sample_", ele))
}

so <- merge(x = so_ls[[1]], y = c(subsetList(so_ls, seq(2, length(so_ls)))), add.cell.ids = sapply(seq(1, length(so_ls)),generateID))

so$Indivi_ID=factor(so$Indivi_ID,levels=Indivi_IDs)
feas = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb")
vlnplot=VlnPlot(so, pt.size=0,features = feas, ncol = 5,group.by="Indivi_ID")+NoLegend()

ggsave(vlnplot,filename = file.path(PATH_O_fig,"0_0_ALS_extra_snRNA_BeforeQC_Vln12252025.png"),width =10,height = 4,dpi=300)

table(so$Indivi_ID) #consistent with library (maybenames by sequencing person)
saveRDS(so, file = file.path(PATH_O_data_seq,"xxx.rds"),compress=F)  #compress can accelerate the process




####################################################################################################
##############################QC process########################################
so<-readRDS(file.path(PATH_O_data_seq,"xxx.rds"))
Idents(so)=so$Indivi_ID
so$Indivi_ID=factor(so$Indivi_ID,levels=Indivi_IDs)
plot1 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(plot1+plot2,filename = file.path(PATH_O_fig,"xxx.png"),width =10,height = 4,dpi=100)



so <- subset(so, subset =
nFeature_RNA > 200 &
nFeature_RNA < 6000 &
nCount_RNA > 200 &
nCount_RNA<20000 &
percent.mt < 10 &
percent.hb<10 &
percent.ribo<5)




so$Indivi_ID=factor(so$Indivi_ID,levels=Indivi_IDs)
feas = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb")
vlnplot=VlnPlot(so, pt.size=0,features = feas, ncol = 5,group.by="Indivi_ID")+NoLegend()
ggsave(vlnplot,filename = file.path(PATH_O_fig,"0_1_ALS_extra_snRNA_AfterQC_Vln12252025.png"),width =12,height = 3,dpi=300)
#scatter plot
plot1 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(plot1+plot2,filename = file.path(PATH_O_fig,"0_1_ALS_human_data_AfterQC_FeaScatter12252025.png"),width =10,height =5,dpi=100)


saveRDS(so, file = file.path(PATH_O_data_seq,"xxx.rds"),compress=F)


####################################################################################################
###############################do doublets by DoubletFinder########################################


so=readRDS(file = file.path(PATH_O_data_seq,"xxx.rds"))
#1 define DR rate
DR_chose<-function(data,DRs){
cell_num=dim(data@meta.data)[1]
DR=cell_num*8*1e-6
print(paste("Cell number:",cell_num,"; DoubletRate:",DR,sep=""))
return(DR)
}

#2  find existed doublets
dim.usage=30
Find_doublet <- function(data){
#optimize pk value
sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
#homotypic doublet proportion estimate
annotations<-data@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  #0.3149338
DoubletRate=DR_chose(data,DRs) ##give rate according to number of cells
nExp_poi <- round(DoubletRate*nrow(data@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#nExp_poi
data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
colnames(data@meta.data)[ncol(data@meta.data)-1] = "doubFind_score"   #change col of "pANNxxxx" to "doubFind_score"
colnames(data@meta.data)[ncol(data@meta.data)] = "doubFind_res"   #change col of "DF.classificaionsxxxx" to "doubFind_res"

data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
colnames(data@meta.data)[ncol(data@meta.data)-1] = "doubFindadj_score"
colnames(data@meta.data)[ncol(data@meta.data)] = "doubFindadj_res"
return(data)
}


subsetList <- function(myList, elementNames) {   #for list merge
  sapply(elementNames, FUN=function(x) myList[[x]])
}


so.list<-SplitObject(so,split.by="Indivi_ID")
so.list<- lapply(X = so.list, FUN = function(x) {
x <- NormalizeData(x)
x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
x <- ScaleData(x)
x <- RunPCA(x)
x <- RunUMAP(x, dims = 1:dim.usage)
x <- FindNeighbors(x,reduction="pca",k.param = 30,  dims = 1:dim.usage)
x <- FindClusters(x, resolution = 0.2) # need tuning, generally larger value for larger dataset to obtain more clusters to find more details, which can be further tuned based on the presented results
x <-Find_doublet(x)
})


dim.usage=30
so <- merge(x = so.list[[1]], y = c(subsetList(so.list, seq(2, length(so.list)))))
so <- NormalizeData(so)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
so <- ScaleData(so)
so <- RunPCA(so,reduction.name = 'PCA_doublets', reduction.key = 'PCA_doublets_')
so <- RunUMAP(so, reduction ='PCA_doublets',  reduction.name = 'UMAP_doublets',dims = 1:dim.usage)
so <- FindNeighbors(so,reduction="PCA_doublets",k.param = 30,  dims = 1:dim.usage)
so <- FindClusters(so, resolution = 0.2)


p1=DimPlot(so,group.by="doubFind_res",reduction="UMAP_doublets",raster= FALSE)
p2=DimPlot(so,group.by="doubFindadj_res",reduction="UMAP_doublets",raster= FALSE)
p3=DimPlot(so,label=T,reduction="UMAP_doublets",raster= FALSE)
p4=DimPlot(so,label=FALSE,group.by="Indivi_ID",reduction="UMAP_doublets",raster= FALSE)
ggsave(p1+p2+p3+p4,filename = file.path(PATH_O_fig,"0_2_ALS_extra_human_data_Umap_doublets12252025.png"),width =13, height = 8,dpi=100)

saveRDS(so, file = file.path(PATH_O_data_seq,"xxx.rds"),compress=F)







