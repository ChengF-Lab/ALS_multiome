library(ArchR)
library('Cairo')
datemark="240501"
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


projMulti4 <- loadArchRProject(path =input_multiome)
projMulti4 <- saveArchRProject(ArchRProj = projMulti4, outputDirectory = output_multiome, overwrite = TRUE, load = TRUE)


cluster_inte="Clusters_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined_reso0.2"
umap_inte="UMAP_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined"
cluster_RNA="Clusters_RNA_0.2_Har_reso0.2"
umap_RNA="UMAP_LSI_RNA_0.2_Har"
cluster_atac="Clusters_ATAC_0.2_Har_reso0.2"
umap_ATAC="UMAP_LSI_ATAC_0.2_Har"



MARKERS_new=c("PLP1", "MBP", "ST18",       # oligodendrocyte
"SLC1A2", "AQP4", "ATP1B2", "GFAP", "CD44","TNC",  #DAA# astrocyte
"VCAN", "CSPG4","OLIG1",      #"PDGRFA",#OPC
"SYT1", "RBFOX1", "CNTNAP2","MAP1B",   #"PDE1A",  #neuron
"CIITA", "CX3CR1", "CSF1R","P2RY12", # microglia
"CLDN5","ABCB1", "FLT1", # Endothelial
"DCN","PDGFRB" ,# , Fib
"SLC17A7", "GRIN1", "SNAP25", #eXNn
"GAD1","GAD2"  #eXNn, and InN
)



posi=match(MARKERS_new,rowData(denovo_markers_genescore0)$name)
heatmap_wellknown_markers_genescore0 <- plotMarkerHeatmap(
    seMarker = denovo_markers_genescore0,
    cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
    subsetMarkers= posi,transpose=TRUE,
    labelMarkers = MARKERS_new, labelRows=FALSE,binaryClusterRows=TRUE,clusterCols=TRUE
)


posi=match(MARKERS_new,rowData(denovo_markers_geneexpr0)$name)
heatmap_wellknown_markers_geneexpr0 <- plotMarkerHeatmap(
    seMarker = denovo_markers_geneexpr0,
    cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
    subsetMarkers= posi,transpose=TRUE,
    labelMarkers = MARKERS_new, labelRows=FALSE,binaryClusterRows=TRUE,clusterCols=TRUE
)
svglite(file.path(output_path,"ALS_wellknown_markergenes_gene_score0_240501.svg"), width =10, height = 5)
ComplexHeatmap::draw(heatmap_wellknown_markers_genescore0, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()
svglite(file.path(output_path,"ALS_wellknown_markergenes_gene_expr0_240501.svg"), width = 10, height = 5)
ComplexHeatmap::draw(heatmap_wellknown_markers_geneexpr0, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()




#rename the cell type name
cluster_inte="Clusters_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined_reso0.2"
umap_inte="UMAP_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined"
cluster_RNA="Clusters_RNA_0.2_Har_reso0.2"
umap_RNA="UMAP_LSI_RNA_0.2_Har"
cluster_atac="Clusters_ATAC_0.2_Har_reso0.2"
umap_ATAC="UMAP_ATAC_0.2_Har2"

cluster0_inte_list <-projMulti4$Clusters_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined_reso0.2
mapping_vector_subcelltype<- c(
"C1" = "MG_1",
"C2" = "Endo",
"C3" = "MG_2",
"C4" = "Astro_3",
"C5" = "ODC_1",
"C6" = "ODC_3",
"C7" = "ODC_4",
"C8" = "ODC_2",
"C9" = "ODC_5",
"C10"="Astro_2",
"C11" = "Astro_1",
"C12"="Neuron",
"C13" = "OPC"
)
projMulti4$Clusters_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined_reso0.2_subcelltype <- mapping_vector_subcelltype[cluster0_inte_list]


mapping_vector_celltype<- c(
"C1" = "MG",
"C2" = "Endo",
"C3" = "MG",
"C4" = "Astro",
"C5" = "ODC",
"C6" = "ODC",
"C7" = "ODC",
"C8" = "ODC",
"C9" = "ODC",
"C10"="Astro",
"C11" = "Astro",
"C12"="Neuron",
"C13" = "OPC"
)

projMulti4$Clusters_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined_reso0.2_celltype <- mapping_vector_celltype[cluster0_inte_list]

projMulti4$Celltype=projMulti4$Clusters_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined_reso0.2_celltype
projMulti4$SubCelltype=projMulti4$Clusters_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined_reso0.2_subcelltype

p1 <- plotEmbedding(projMulti4, name = "Celltype", embedding = "UMAP_LSI_ATAC_0.2_Har", size = 0.1, baseSize=14,labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(projMulti4, name = "Celltype", embedding = "UMAP_LSI_RNA_0.2_Har", size = 0.1, baseSize=14,labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(projMulti4, name = "Celltype", embedding = "UMAP_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined", size = 0.1, baseSize=14,labelAsFactors=F, labelMeans=F)
svglite(file.path(output_path,paste0("ALS_UMAP_after_harmony_Cluster240501.svg")), width =14, height = 5)
combined_plot <- ggarrange(p1, p2, p3, ncol = 3, nrow = 1)
print(combined_plot)
dev.off()

p1 <- plotEmbedding(projMulti4, name = "SubCelltype", embedding = "UMAP_LSI_ATAC_0.2_Har", size = 0.1, baseSize=14,labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(projMulti4, name = "SubCelltype", embedding = "UMAP_LSI_RNA_0.2_Har", size = 0.1, baseSize=14,labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(projMulti4, name = "SubCelltype", embedding = "UMAP_LSI_ATAC_0.2_Har_LSI_RNA_0.2_Har_Combined", size = 0.1, baseSize=14,labelAsFactors=F, labelMeans=F)
svglite(file.path(output_path,paste0("ALS_UMAP_after_harmony_Subcluster240501.svg")), width =14, height = 5)
combined_plot <- ggarrange(p1, p2, p3, ncol =3, nrow = 1)
print(combined_plot)
dev.off()

#######calculate denovo_markers
denovo_markers_geneexpr_celltype <- getMarkerFeatures(ArchRProj = projMulti4,
                        groupBy = "Celltype",  #
                        useMatrix = "GeneExpressionMatrix",
                        bias = c("TSSEnrichment", "log10(nFrags)","log10(Gex_nUMI)"))
saveRDS(denovo_markers_geneexpr_celltype, file =file.path(output_path,"ALS_celltype_denovo_marker_geneexpr_240501.rds"))

heatmap_denovo_markers_geneexpr_celltype <- plotMarkerHeatmap(
  seMarker = denovo_markers_geneexpr_celltype,
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  transpose = FALSE,
    nLabel = 10
)
svglite(file.path(output_path,"ALS_celltype_denovo_markergenes_gene_expr_240501.svg"), width = 5, height = 10)
ComplexHeatmap::draw(heatmap_denovo_markers_geneexpr_celltype, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()




denovo_markers_genescore_celltype <- getMarkerFeatures(ArchRProj = projMulti4,
                        groupBy = "Celltype",  #using the combined results from separately batch-corrected results
                        useMatrix = "GeneScoreMatrix",
                        bias = c("TSSEnrichment", "log10(nFrags)","log10(Gex_nUMI)"))

saveRDS(denovo_markers_genescore_celltype, file =file.path(output_path,"ALS_celltype_denovo_markers_genescore_240501.rds"))

heatmap_denovo_markers_genescore_celltype <- plotMarkerHeatmap(
  seMarker = denovo_markers_genescore_celltype,
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  transpose = FALSE,
  nLabel = 10
)
svglite(file.path(output_path,"ALS_celltype_denovo_markergenes_gene_score_240501.svg"), width = 5, height = 10)
ComplexHeatmap::draw(heatmap_denovo_markers_genescore_celltype, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


###add mete data to the ArchR object
sample_info<-read.csv(file.path("xxx","ALS_sequencing_data/Sample_info230920.tsv"),header=TRUE,sep="\t")
sample_info$IndiviID=sample_info$subject.ID
meta_data_multiome=data.frame(projMulti4@cellColData)
posi=match(meta_data_multiome$Sample,sample_info$subject.ID,)
projMulti4$Indivi_ID=meta_data_multiome$Sample
projMulti4$Sample_ID=meta_data_multiome$Sample
projMulti4$Tissue=as.character(sample_info$Tissue[posi])
projMulti4$Sex=as.character(sample_info$SEX[posi])
projMulti4$Fibrinigen=as.character(sample_info$Fibrinigen[posi])
projMulti4$Age_onset=as.character(sample_info$sage.at.onset[posi])
projMulti4$Age_death=as.character(sample_info$age.at.death[posi])
projMulti4$disease_furation_yrs=as.character(sample_info$disease.duration..yrs.[posi])
projMulti4$death_fixation_interval_hrs=as.character(sample_info$death.fixation.interval..hrs.[posi])
projMulti4$date_post_fixation_MRI=as.character(sample_info$date.of.post.fixation.MRI[posi])
projMulti4$ALS_group=as.character(sample_info$fALS.sALS[posi])
projMulti4$onset_location=as.character(sample_info$onset.location..E.distal.limb.L.proximal.limb.[posi])


allsamples = c("ALS1","ALS9","ALS11","ALS17","ALS3","ALS4","ALS8","ALS16") #first four: ALS; last forth: HC
celltypes=c("ODC","Astro","OPC","MG","Neuron","Endo")
subcelltypes=c("ODC_1",   "ODC_2",   "ODC_3","ODC_4","ODC_5","Astro_1", "Astro_2","Astro_3","OPC",
"MG_1","MG_2","Neuron","Endo")
projMulti4$Fibrinigen_new <- paste0(projMulti4$Fibrinigen," fibrinigen")
projMulti4$Fibrinigen=paste0(projMulti4$Fibrinigen,"_Fibrinigen")
projMulti4$Fibrinigen_new = gsub("LOW", "low", projMulti4$Fibrinigen_new)
projMulti4$Fibrinigen_new = gsub("HIGH", "high", projMulti4$Fibrinigen_new)

projMulti4$Celltype_Fibrinigen = paste0(projMulti4$Celltype,"_",projMulti4$Fibrinigen)
projMulti4$SubCelltype_Fibrinigen = paste0(projMulti4$SubCelltype,"_",projMulti4$Fibrinigen)
projMulti4 <- saveArchRProject(ArchRProj = projMulti4, outputDirectory = getOutputDirectory(projMulti4), overwrite = TRUE, load = TRUE)

