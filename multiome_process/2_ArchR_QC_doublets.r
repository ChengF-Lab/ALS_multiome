library(ArchR)
library('Cairo')
addArchRGenome("hg38")
addArchRThreads(50)
addArchRLocking(locking = TRUE)
set.seed(123)
library(ggplot2)
library(cowplot)


input_path="xxx"
output_path="xxx"


allsamples = c("ALS1","ALS9","ALS11","ALS17","ALS3","ALS4","ALS8","ALS16")
projMulti2 <- loadArchRProject(path =input_path)
projMulti2 <- saveArchRProject(ArchRProj = projMulti2, outputDirectory = output_path, overwrite = TRUE, load = TRUE)
meta_data=data.frame(projMulti2@cellColData)
meta_data=data.frame(projMulti2@cellColData)

QC_plot_details = function(data, features,mark,file_path) {
fig_list=list()
for (i in 1:length(features)){
feature=features[i]
feature=ensym(feature)
plot=ggplot(data, aes(x = Sample, y = !!feature)) +
  geom_violin()
fig_list[[i]]=plot
}
p=plot_grid(plotlist=fig_list,ncols=5)
ggsave(p,filename = file.path(output_path,paste0("ALS_QC_vlnplot_",mark,"231210.svg")),width =20,height =15,bg="white")
return(fig_list)
}

features=c("TSSEnrichment" ,   "ReadsInTSS" ,      "ReadsInPromoter" ,
"ReadsInBlacklist" ,"PromoterRatio"  ,   "NucleosomeRatio" ,
"nMultiFrags"   ,   "nMonoFrags",       "nFrags"     ,      "nDiFrags"   ,
"BlacklistRatio",   "Gex_nUMI"  ,       "Gex_nGenes"  ,     "Gex_MitoRatio"  ,
"Gex_RiboRatio" )

p=QC_plot_details(meta_data,features,"ALS_nuclei_quality_beforeQC_",output_path)


#############QC process
idxPass <- which(
projMulti2$TSSEnrichment < 20  &    #ATAC quality
projMulti2$BlacklistRatio<0.02 &
projMulti2$NucleosomeRatio<5  &
projMulti2$nFrags<5e4  &

projMulti2$Gex_nUMI >200 &     #RNA quality
projMulti2$Gex_nUMI < 1e4  &
projMulti2$Gex_nGenes >200 &
projMulti2$Gex_nGenes <4000 &
projMulti2$Gex_MitoRatio <0.1 &
projMulti2$Gex_RiboRatio <0.1
)


cellsPass <- projMulti2$cellNames[idxPass]
projMulti2=projMulti2[cellsPass, ]
meta_data=data.frame(projMulti2@cellColData)
p=QC_plot_details(meta_data,features,"ALS_nuclei_quality_afterQC_",output_path)
projMulti2 <- saveArchRProject(ArchRProj = projMulti2, outputDirectory = getOutputDirectory(projMulti2), overwrite = TRUE, load = TRUE)


#####################doublet removing###########
projMulti2 <- addDoubletScores(projMulti2, useMatrix = "TileMatrix", k=10, knnMethod = "UMAP", LSIMethod =1, force = FALSE,logFile = createLogFile("addDoubletScores"))
projMulti2 <- saveArchRProject(ArchRProj = projMulti2, outputDirectory = getOutputDirectory(projMulti2), overwrite = TRUE, load = TRUE)

