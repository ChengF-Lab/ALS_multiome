  library(Seurat)
  library(SeuratDisk)
  library(tidyverse)
  library(data.table)
  library(clusterProfiler)
  library(dplyr)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(enrichplot)
  library(ReactomePA)
  library(ArchR)
  library(rGREAT)
  library(bedr)
  library(svglite)

input_path="xxx"
output_path="xxx"

celltypes <- c("ODC", "Astro", "OPC", "MG", "Neuron", "Endo")
subcelltypes <- c("ODC_1", "ODC_2", "ODC_3", "ODC_4", "ODC_5", "Astro_1", "Astro_2", "Astro_3",
                  "OPC", "MG_1", "MG_2", "Neuron", "Endo")
num_celltypes <- length(celltypes)
num_subcelltypes <- length(subcelltypes)
OrgDb_mm="org.Mm.eg.db"
OrgDb_hs="org.Hs.eg.db"
allpeaks=read.csv("xxx.tsv",sep="\t")

allpeaks = allpeaks[,c(1,2,3)]
allpeaks.loc = paste0(allpeaks[,1],":",allpeaks[,2],"-",allpeaks[,3])
allpeaks.loc = bedr.sort.region(allpeaks.loc)
peak_denovo=readRDS("xxx.rds")
peakList <- getMarkers(peak_denovo, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")

for (celltype in celltypes[-6]){   #on Endo_specifc peaks, so delete this one
df=peakList[[celltype]]
df.loc = paste0(df[,1],":",df[,3],"-",df[,4])
 df.loc = bedr.sort.region(df.loc)
 df.intersect = bedr.join.region(df.loc,allpeaks.loc)
 df = unique(df.intersect[,2:4])
 colnames(df) = colnames(allpeaks)
 df = df[df$seqnames != ".",]
 df$start = as.numeric(df$start)
 df$end = as.numeric(df$end)
 job = submitGreatJob(df,bg=allpeaks,species="hg38")
 tb = getEnrichmentTables(job)
availableCategories(job)
df=peakList[[celltype]]
tab <- getEnrichmentTable(job,catogery="GO")
tab <- getEnrichmentTable(res)
saveRDS(tb,file.path(output_path,paste0("xxx_rGREAT.rds")))
cat(celltype,"is done")
}




