library(ArchR)
library('Cairo')
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


allsamples = c("ALS1","ALS9","ALS11","ALS17","ALS3","ALS4","ALS8","ALS16")
celltypes=c("ODC","Astro","OPC","MG","Neuron","Endo")
subcelltypes=c("ODC_1",   "ODC_2",   "ODC_3","ODC_4","ODC_5","Astro_1", "Astro_2","Astro_3","OPC",
"MG_1","MG_2","Neuron","Endo")
num_celltypes=length(celltypes)
num_subcelltypes=length(subcelltypes)
pal_subcelltypes=c("#667FE1", "#ACD4EC", "#6F99AD","#64AAD2", "#8AB6E9FF","#FE9586" ,"#D87B8B","#D898B9","#EBCC78", "#5DA59E" ,"#56BC9B","#A48BCA", "#8C7F5F")
pal_celltypes=c("#8AB6E9FF","#D87B8B","#EBCC78", "#5DA59E" ,"#A48BCA", "#8C7F5F")
names(pal_celltypes)=celltypes
names(pal_subcelltypes)=subcelltypes
cluster_inte="Celltype"
subcluster_inte="SubCelltype"
cluster_inte="Clusters_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined_reso0.2"
umap_inte="UMAP_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined"
reduction_inte="LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined"
cluster_RNA="Clusters_RNA_0.2_Har_reso0.2"
umap_RNA="UMAP_LSI_RNA_0.2_Har"
cluster_atac="Clusters_ATAC_0.2_Har_reso0.2"
umap_ATAC="UMAP_LSI_ATAC_0.2_Har"


projMulti9 <- loadArchRProject(path =input_path)
projMulti9 <- saveArchRProject(ArchRProj = projMulti9, outputDirectory = output_path, overwrite = TRUE, load = TRUE)
motifPositions <- getPositions(projMulti9,name = "Motif")
ALS_genes <- c("C9orf72", "SOD1", "FUS", "TDP-43", "MATR3", "ATXN2","UBQLN2","TBK1","NEK1","OPTN")
motifs <- c("SPIB", "PITX2","NFIC","ZBTB3","SPI1","SOX13","SOX4")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

if(is.null(projMulti9@projectMetadata$GroupCoverages$SubCelltype)){
  projMulti9 <- addGroupCoverages(ArchRProj = projMulti9, groupBy = "SubCelltype")
}
seFoot <- getFootprints(
  ArchRProj = projMulti9,
  positions = motifPositions[markerMotifs],
  groupBy = "SubCelltype"
)

####way1:  Subtracting the Tn5 Bias
p=plotFootprints(
  seFoot = seFoot,
  ArchRProj = projMulti9,
  normMethod = "Subtract",
  plotName = "9_3_ALS_TF_Footprints-Subtract-Bias_Subtract",
  addDOC = FALSE,
  smoothWindow = 5,
  plot=TRUE,
  pal=pal_subcelltypes
)

###way2:  Dividing by the Tn5 Bias
plotFootprints(
  seFoot = seFoot,
  ArchRProj = projMulti9,
  normMethod = "Divide",
  plotName = "9_3_ALS_TF_Footprints-Divide-Bias_Divide",
  addDOC = FALSE,
  smoothWindow = 5,
    pal=pal_subcelltypes
)
seTSS <- getFootprints(
  ArchRProj = projMulti9,
  positions = GRangesList(TSS = getTSS(projMulti9)),
  groupBy = "SubCelltype",
  flank = 2000
)
plotFootprints(
  seFoot = seTSS,
  ArchRProj = projMulti9,
  normMethod = "None",
  plotName = "9_3_ALS_TSS_footprint_No-Normalization_SubCelltype",
  addDOC = FALSE,
  flank = 2000,
  flankNorm = 100,
      pal=pal_subcelltypes
)
projMulti9 <- saveArchRProject(ArchRProj = projMulti9, outputDirectory = getOutputDirectory(projMulti9), overwrite = TRUE, load = TRUE)

















