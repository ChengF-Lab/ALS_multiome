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

projMulti10 <- loadArchRProject(path =input_path)
projMulti10 <- saveArchRProject(ArchRProj = projMulti10, outputDirectory = output_path, overwrite = TRUE, load = TRUE)

##############################################################################################################################
##PART_I: CO-accessibility

projMulti10 <- addCoAccessibility(
    ArchRProj = projMulti10,
    reducedDims = reduction_inte
)
cA <- getCoAccessibility(
    ArchRProj = projMulti10,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = FALSE
)
saveRDS(cA,file=file.path(output_path,"10_3_ALS_co_accessibility_returnLoops_FALSE_corcut0.5_reso_1_240501_SubCelltype_big_maxPeaks.rds"))
cA <- getCoAccessibility(
    ArchRProj = projMulti10,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = TRUE
)

saveRDS(cA,file=file.path(output_path,"10_3_ALS_co_accessibility_returnLoops_true_corcut0.5_reso_1_240501_SubCelltype_big_maxPeaks.rds"))
cA[[1]]
##filter some strong results BASED ON fdr and other metrics
cALoops <- cA[[1]]
cALoops <- cALoops[cALoops$FDR < 10^-10]
cALoops <- cALoops[rowMins(cbind(cALoops$VarQuantile1,cALoops$VarQuantile2)) > 0.35]
cA <- getCoAccessibility(
    ArchRProj = projMulti10,
    corCutOff = 0.5,
    resolution = 1000,   #can change resolution
    returnLoops = TRUE
)
saveRDS(cA,file=file.path(output_path,"10_3_ALS_co_accessibility_returnLoops_true_corcut0.5_reso_1000_240501_Celltype_big_maxPeaks.rds"))


##############################################################################################################################
##PART_II: P2G
projMulti10 <- addPeak2GeneLinks(
    ArchRProj = projMulti10,
    reducedDims = reduction_inte,
    useMatrix="GeneExpressionMatrix"
)
p2g <- getPeak2GeneLinks(
    ArchRProj = projMulti10,
    corCutOff = 0.45,   ####
    resolution = 1,
    returnLoops = FALSE
)
p2g$geneName <- mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]
p2g$peakName <- (metadata(p2g)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2g$idxATAC]
p2g

write.table(p2g,file=file.path(output_path,"10_3_ALS_p2g_results_returnLoops_false_corcut_0.45_reso_1_240501_SubCelltype_big_maxPeaks.tsv"),sep="\t")
projMulti10 <- saveArchRProject(ArchRProj = projMulti10, outputDirectory = getOutputDirectory(projMulti10), overwrite = TRUE, load = TRUE)


