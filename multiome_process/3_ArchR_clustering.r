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



input_path="xxx"
output_path="xxx"

allsamples = c("ALS1","ALS9","ALS11","ALS17","ALS3","ALS4","ALS8","ALS16")


projMulti3 <- loadArchRProject(path =input_path)
projMulti3 <- saveArchRProject(ArchRProj = projMulti3, outputDirectory = output_path, overwrite = TRUE, load = TRUE)

cell_count=table(projMulti3$Sample)
cell_count=cell_count[allsamples]

projMulti3 <- filterDoublets(projMulti3)
cell_count=table(projMulti3$Sample)
cell_count=cell_count[allsamples]
meta_data=data.frame(projMulti3@cellColData)





#######ATAC dimension reduction: LSI-ATAC
projMulti3 <- addIterativeLSI(
  ArchRProj = projMulti3,
  clusterParams = list(
    resolution = 0.2,
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures=25000,   #peak is the feature
  saveIterations = FALSE,
  useMatrix = "TileMatrix",
  depthCol = "nFrags",
  firstSelection = "Top",
  name = "LSI_ATAC_0.2",
  iterations=5,
dimsToUse=1:30,force = TRUE
)


projMulti3@reducedDims

#############LSI-RNA
#seRNA dimension reduction
projMulti3 <- addIterativeLSI(
  ArchRProj = projMulti3,
  clusterParams = list(
    resolution = 0.2,
    sampleCells = 10000,
    n.start = 10
  ),
    varFeatures = 5000,   #less than half of gene #
    saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix",
  depthCol = "Gex_nUMI",  #these values are used to minimize the related biases in the reduction related. For scATAC we recommend "nFrags" and for scRNA we recommend "Gex_nUMI".
  firstSelection = "Var",
  binarize = FALSE,
  name = "LSI_RNA_0.2",
iterations=5,
force = TRUE
)





###################batch effect correction using harmony
projMulti3 <- addHarmony(
    ArchRProj = projMulti3,
    reducedDims = "LSI_ATAC_0.2",
    name = "LSI_ATAC_0.2_Har",
    groupBy = "Sample"
)

projMulti3 <- addHarmony(
    ArchRProj = projMulti3,
    reducedDims = "LSI_RNA_0.2",
    name = "LSI_RNA_0.2_Har",
    groupBy = "Sample"
)



projMulti3 <- addCombinedDims(projMulti3, reducedDims = c("LSI_ATAC_0.2_Har", "LSI_RNA_0.2_Har"), name =  "LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined")  #combine ATAC_har and RNA_har

projMulti3 <- addUMAP(ArchRProj = projMulti3,reducedDims = "LSI_ATAC_0.2_Har",name = "UMAP_LSI_ATAC_0.2_Har", nNeighbors = 30, minDist = 0.5, metric = "cosine")
projMulti3 <- addUMAP(ArchRProj = projMulti3,reducedDims = "LSI_RNA_0.2_Har",name = "UMAP_LSI_RNA_0.2_Har", nNeighbors = 30, minDist = 0.5, metric = "cosine")
projMulti3 <- addUMAP(ArchRProj = projMulti3,reducedDims = "LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined",name = "UMAP_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined", nNeighbors = 30, minDist = 0.5, metric = "cosine")

projMulti3 <- addClusters(projMulti3, reducedDims = "LSI_ATAC_0.2_Har", name = "Clusters_ATAC_0.2_Har_reso0.1", resolution = 0.1, force = TRUE)
projMulti3 <- addClusters(projMulti3, reducedDims = "LSI_ATAC_0.2_Har", name = "Clusters_ATAC_0.2_Har_reso0.15", resolution = 0.15, force = TRUE)
projMulti3 <- addClusters(projMulti3, reducedDims = "LSI_ATAC_0.2_Har", name = "Clusters_ATAC_0.2_Har_reso0.2", resolution = 0.2, force = TRUE)
projMulti3 <- addClusters(projMulti3, reducedDims = "LSI_RNA_0.2_Har", name = "Clusters_RNA_0.2_Har_reso0.1", resolution = 0.1, force = TRUE)
projMulti3 <- addClusters(projMulti3, reducedDims = "LSI_RNA_0.2_Har", name = "Clusters_RNA_0.2_Har_reso0.15", resolution = 0.15, force = TRUE)
projMulti3 <- addClusters(projMulti3, reducedDims = "LSI_RNA_0.2_Har", name = "Clusters_RNA_0.2_Har_reso0.2", resolution = 0.2, force = TRUE)
projMulti3 <- addClusters(projMulti3, reducedDims = "LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined", name = "Clusters_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined_reso0.1", resolution = 0.1, force = TRUE)
projMulti3 <- addClusters(projMulti3, reducedDims = "LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined", name = "Clusters_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined_reso0.15", resolution = 0.15, force = TRUE)
projMulti3 <- addClusters(projMulti3, reducedDims = "LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined", name = "Clusters_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined_reso0.2", resolution = 0.2, force = TRUE)

projMulti3 <- saveArchRProject(ArchRProj = projMulti3, outputDirectory = getOutputDirectory(projMulti3), overwrite = TRUE, load = TRUE)
