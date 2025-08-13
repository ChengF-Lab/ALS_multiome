library(Seurat)
library(SeuratDisk)
library(tidyverse)
set.seed(123)
library(svglite)
rm(list = ls())
options(stringsAsFactors = F)
library(data.table)
library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
library(enrichplot)
library(ReactomePA)


input_path="xxx"
output_path="xxx"

allsamples <- c("ALS1", "ALS9", "ALS11", "ALS17", "ALS3", "ALS4", "ALS8", "ALS16")
celltypes <- c("ODC", "Astro", "OPC", "MG", "Neuron", "Endo")
subcelltypes <- c("ODC_1", "ODC_2", "ODC_3", "ODC_4", "ODC_5", "Astro_1", "Astro_2", "Astro_3",
                  "OPC", "MG_1", "MG_2", "Neuron", "Endo")
num_celltypes <- length(celltypes)
num_subcelltypes <- length(subcelltypes)
groups <- c("LOW_Fibrinigen", "HIGH_Fibrinigen")
compare_group1 <- c("LOW_Fibrinigen")
compare_group2 <- c("HIGH_Fibrinigen")
compare_groups <- paste0(compare_group2, "_vs_", compare_group1)
num_groups <- length(groups)
num_com_groups <- length(compare_group2)

OrgDb_mm="org.Mm.eg.db"
OrgDb_hs="org.Hs.eg.db"

log2FC_cutoff=0.25
FDR_cutoff=0.05

DEGs_data=read.csv("xxx.tsv",sep="\t")
DEGs_data$compare_group=compare_groups
DEGs_data$gene=DEGs_data$name

for (celltype in subcelltypes) {
   for (i in 1:length(compare_group1)) {
     group1 <- compare_group1[i]
     group2 <- compare_group2[i]
     comparison=paste0(group2,"_vs_",group1)
     mark=paste0(celltype,"_",comparison)
     data=DEGs_data[DEGs_data$Celltype == celltype,   ]
     data=data[data$compare_group == comparison, ]
       data=data[,c("gene","Log2FC")]
       rownames(data)=1:(dim(data)[1])
        geneList<- data$Log2FC
        names(geneList)= data$gene
        geneList=sort(geneList,decreasing = T)



 if (length(geneList)>20) {
 cat(celltype,"have enough DEGs")
 gse_results <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Hs.eg.db,
  ont          = "ALL",
  keyType      = "SYMBOL",
  nPerm        = 1000,
  minGSSize    = 2,
  maxGSSize    = 500,
   pvalueCutoff = 0.1,
   pAdjustMethod = "BH",
  verbose      = FALSE,
  seed = 6, by = "fgsea")
 if (dim(gse_results@result)[1]>0){
require(DOSE)
show_term=min(10,dim(gse_results@result)[1])
 p=dotplot(gse_results,
  color = "NES",
  split=".sign",title=mark, showCategory=show_term)+facet_grid(~.sign)#点状图
 fig_filename=file.path(output_path,paste0("12_ALS_subcluster_DEG_",mark,"_enrichemnt_clusterprofiler_goGSEA_dotplot_241020.svg"))
 svglite(fig_filename)
 print(p)
 dev.off()
 data_file=file.path(PATH_O_data_enrich_SubCluster,paste0("12_ALS_subcluster_DEG_",mark,"_enrichemnt_clusterprofiler_goGSEA_241020.rds"))
 saveRDS(gse_results,file=data_file)
 }
}
 cat(mark, "gsego enriched termes #(p.adj<0.1)",dim(gse_results@result)[1],"\n")

 #######A: KEGG ANALYSIS
       gene.df <- bitr(data$gene,
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
#                 OrgDb =org.Mm.eg.db)   #mouse data
                 OrgDb = org.Hs.eg.db)
       gene.df = gene.df[!duplicated(gene.df[c("ENTREZID")]),]
       data2=merge(data,gene.df,by.x="gene",by.y="SYMBOL")
       geneList<-data2$Log2FC
        names(geneList)=data2$ENTREZID
        geneList=sort(geneList,decreasing = T)
 if (length(geneList)>20) {
 gse_results <- gseKEGG(geneList, organism = "hsa",
 			  pvalueCutoff = 0.1,
 			  pAdjustMethod = "BH",
#               nPerm = 1000,
              minGSSize = 2,
              maxGSSize = 500,
              verbose = TRUE, seed = FALSE,
  				keyType  = "ncbi-geneid",
              by = "fgsea")

 if (dim(gse_results@result)[1]>0){
show_term=min(10,dim(gse_results@result)[1])
 p=dotplot(gse_results,
 color = "NES",
 split=".sign",title=mark,showCategory=show_term)+facet_grid(~.sign)#点状图
fig_filename=file.path(output_path,paste0("12_ALS_subcluster_DEG_",mark,"_enrichemnt_clusterprofiler_keggGSEA_dotplot_241020.svg"))
svglite(fig_filename)
print(p)
dev.off()
 gse_results = setReadable(gse_results, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
 data_file=file.path(output_path,paste0("12_ALS_subcluster_DEG_",mark,"_enrichemnt_clusterprofiler_keggGSEA_241020.rds"))
 saveRDS(gse_results,file=data_file)
 }
}
 cat(mark, "gsekegg enriched terms #(p<0.1)",dim(gse_results@result)[1],"\n")
   }
 }







