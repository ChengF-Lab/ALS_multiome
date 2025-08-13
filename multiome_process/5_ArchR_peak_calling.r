library(ArchR)
library('Cairo')
addArchRGenome("hg38")
addArchRThreads(50)
addArchRLocking(locking = TRUE)
set.seed(123)
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


cluster_inte="Celltype"
subcluster_inte="SubCelltype"
umap_inte="UMAP_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined"
reduction_inte="LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined_reso0.2"

cluster_RNA="Clusters_LSI_RNA_0.2_Har_reso0.2_celltype"
subcluster_RNA="Clusters_LSI_RNA_0.2_Har_reso0.2_subcelltype"
umap_RNA="UMAP_LSI_RNA_0.2_Har"

cluster_atac="Clusters_LSI_ATAC_0.2_Har_reso0.2_celltype"
cluster_atac="Clusters_LSI_ATAC_0.2_Har_reso0.2_subcelltype"
umap_ATAC="UMAP_LSI_ATAC_0.2_Har"



projMulti5 <- loadArchRProject(path ="/home/doul2/beegfs/doul2/Work/ALS/Archr_multiome_output240501/5_2_multiome_psedobulk_peakcaling_SubCelltype") #to avoid rerun add group, here we directly use previous archr object has finished add grou file
projMulti5 <- saveArchRProject(ArchRProj = projMulti5, outputDirectory = input_path, overwrite = TRUE, load = TRUE)

projMulti5$SubCelltype_Fibrinigen=paste0(projMulti5$SubCelltype,"_",projMulti5$Fibrinigen)
projMulti5 <- addGroupCoverages(ArchRProj = projMulti5, groupBy = "SubCelltype_Fibrinigen",
  minCells = 40,
  maxCells = 500,
  force = TRUE)

pathToMacs2 <- findMacs2()
projMulti5 <- addReproduciblePeakSet(
    ArchRProj = projMulti5,
    maxPeaks = 200000,
    groupBy = "SubCelltype_Fibrinigen",
    pathToMacs2 = pathToMacs2
)

projMulti5_PeakSet=getPeakSet(projMulti5)
    projMulti5_PeakSet$cluster <- names(projMulti5_PeakSet)
    projMulti5_PeakSet_df <- data.frame(projMulti5_PeakSet)
      nrow(projMulti5_PeakSet)
write.table(projMulti5_PeakSet_df,file.path(output_path,"5_3_ALS_peakset_SubCelltype_via_addGroupCoverages_SubCelltype_Fibrinigen240920_big_maxPeaks.tsv"),sep="\t",quote=FALSE)
saveRDS(projMulti5_PeakSet,file.path(output_path,"5_3_ALS_peakset_SubCelltype_via_addGroupCoverages_SubCelltype_Fibrinigen240920_big_maxPeaks.rds"))
projMulti5 <- addPeakMatrix(projMulti5)
######################
markersPeaks <- getMarkerFeatures(
    ArchRProj = projMulti5,
    useMatrix = "PeakMatrix",
    groupBy = "SubCelltype",
  bias = c("TSSEnrichment", "log10(nFrags)","log10(Gex_nUMI)"),
  testMethod = "wilcoxon"
)
markerList_peak<- getMarkers(markersPeaks,cutOff = "abs(Log2FC) >= 0",returnGR = TRUE)
saveRDS(markersPeaks, file =file.path(output_path,"5_3_ALS_denovomarker_peak240501_SubCelltype_via_addGroupCoverages_SubCelltype_Fibrinigen_big_maxPeaks.rds"))
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks,
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5",
  transpose = FALSE,
  nLabel=5
)
svglite(file.path(output_path,"5_3_ALS_denovomarker_peak_heatmap_240501_SubCelltype_via_addGroupCoverages_SubCelltype_Fibrinigen_big_maxPeaks.svg"), width = 4, height = 7)
ComplexHeatmap::draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

projMulti5 <- saveArchRProject(ArchRProj = projMulti5, outputDirectory = getOutputDirectory(projMulti5), overwrite = TRUE, load = TRUE)













