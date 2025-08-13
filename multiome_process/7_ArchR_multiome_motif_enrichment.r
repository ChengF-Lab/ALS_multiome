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

allsamples = c("ALS1","ALS9","ALS11","ALS17","ALS3","ALS4","ALS8","ALS16") #f
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


projMulti7 <- loadArchRProject(path =input_path)
projMulti7 <- saveArchRProject(ArchRProj = projMulti7, outputDirectory = output_path, overwrite = TRUE, load = TRUE)

projMulti7 <- addMotifAnnotations(ArchRProj = projMulti7, motifSet = "cisbp", name = "Motif")
pSet <- getPeakSet(ArchRProj = projMulti7)
pSet$name <- paste(seqnames(pSet), start(pSet), end(pSet), sep = "_")
matches <- getMatches(ArchRProj = projMulti7, name = "Motif")
rownames(matches) <- paste(seqnames(matches), start(matches), end(matches), sep = "_")
matches <- matches[pSet$name]

############Motif Enrichment in denovo-peaks
DARs_file=file.path(input_path,paste0("5_3_ALS_denovomarker_peak240501_Celltype_via_addGroupCoverages_SubCelltype_Fibrinigen_big_maxPeaks.rds"))
idata = readRDS(DARs_file)
enrichMotifs <- peakAnnoEnrichment(
    seMarker = idata,
    ArchRProj = projMulti7,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
motif_enrich_file=file.path(output_path,paste0("7_3_ALS_denovo_peaks_motif_enrichment_240501_Celltype_via_addGroupCoverages_SubCelltype_Fibrinigen_big_maxPeaks.rds"))
saveRDS(enrichMotifs, file =file.path(motif_enrich_file))
svglite(file.path(output_path,"7_3_ALS_denovo_peaks_motif_enrichments_heatmap_Celltype240920_via_addGroupCoverages_SubCelltype_Fibrinigen_big_maxPeaks.svg"), width = 5, height = 10)
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = FALSE)+ggtitle("Celltype-specific motif based on cs-bp")
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

projMulti7 <- saveArchRProject(ArchRProj = projMulti7, outputDirectory = getOutputDirectory(projMulti7), overwrite = TRUE, load = TRUE)












