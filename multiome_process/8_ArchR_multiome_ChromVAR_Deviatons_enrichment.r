library(ArchR)
library('Cairo')
datemark="240501"
addArchRGenome("hg38")
addArchRThreads(50)
addArchRLocking(locking = TRUE)
set.seed(1234)
library(ggplot2)
library(svglite)
library(ggpubr)
library(svglite)
library(BSgenome.Hsapiens.UCSC.hg38)


input_path="xxx"
output_path="xxx"


allsamples = c("ALS1","ALS9","ALS11","ALS17","ALS3","ALS4","ALS8","ALS16") #f
celltypes=c("ODC","Astro","OPC","MG","Neuron","Endo")
subcelltypes=c("ODC_1",   "ODC_2",   "ODC_3","ODC_4","ODC_5","Astro_1", "Astro_2","Astro_3","OPC",
"MG_1","MG_2","Neuron","Endo")
num_celltypes=length(celltypes)
num_subcelltypes=length(subcelltypes)
cluster_inte="Celltype"
subcluster_inte="SubCelltype"
cluster_inte="Clusters_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined_reso0.2"
umap_inte="UMAP_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined"
reduction_inte="LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined"
cluster_RNA="Clusters_RNA_0.2_Har_reso0.2"
umap_RNA="UMAP_LSI_RNA_0.2_Har"
cluster_atac="Clusters_ATAC_0.2_Har_reso0.2"
umap_ATAC="UMAP_LSI_ATAC_0.2_Har"


projMulti8 <- loadArchRProject(path =input_path)
projMulti8 <- saveArchRProject(ArchRProj = projMulti8, outputDirectory = output_path, overwrite = TRUE, load = TRUE)
if("Motif" %ni% names(projMulti8@peakAnnotation)){
    projMulti8 <- addMotifAnnotations(ArchRProj = projMulti8, motifSet = "cisbp", name = "Motif")
}
projMulti8 <- addBgdPeaks(projMulti8)
projMulti8 <- addDeviationsMatrix(
  ArchRProj = projMulti8,
  peakAnnotation = "Motif",    
  force = TRUE
)
projMulti8 <- saveArchRProject(ArchRProj = projMulti8, outputDirectory = getOutputDirectory(projMulti8), overwrite = TRUE, load = TRUE)








