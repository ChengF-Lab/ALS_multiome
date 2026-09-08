#DEGs enrichment analysis |


library(Seurat)
library(SeuratDisk)
library(DoubletFinder)
library(harmony)
library(Matrix)
library(matrixStats)
library(data.table)
library(dplyr)
library(tidyr)
library(edgeR)
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggrepel)
library(ggpubr)
library(viridis)
library(paletteer)
library(scCustomize)
library(clustree)
library(tidydr)
library(svglite)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)

library(org.Hs.eg.db)
library(org.Mm.eg.db)


file_path <- "xxx"
PATH_O <- file.path(file_path, "ALS_review_comments_output11052025/1_Cellfraction_DEGs_enrichments_out12302025")
PATH_O_data <- file.path(PATH_O,"Data")
PATH_O_data_seq <- file.path(PATH_O_data,"data_seq")
PATH_O_fig <- file.path(PATH_O,"Figures")
PATH_O_fig_polish <- file.path(PATH_O_fig,"Figure_polished")
PATH_O_fig_polish_data <- file.path(PATH_O_fig_polish,"Data")


allsamples = c("ALS1","ALS9","ALS11","ALS17", "ALS7","ALS12",      "ALS3","ALS4","ALS8","ALS16","ALS10","ALS14") #low; high
celltypes=c("ODC","Astro","OPC","MG","Neuron","Endo")
subcelltypes=c("ODC_1",   "ODC_2",   "ODC_3","ODC_4","ODC_5","Astro_1", "Astro_2","Astro_3","OPC",
"MG_1","MG_2","Neuron","Endo")
num_celltypes=length(celltypes)
num_subcelltypes=length(subcelltypes)



subcelltypes_merged=c("ODC_1",   "ODC_2",   "ODC_3","ODC_4","ODC_5","Astro_1", "Astro_2","Astro_3","OPC",
"MG_1","MG_2","ExN","InN","Endo","Pericytes","Fibroblast")
pal_subcelltypes_merged=c("#667FE1", "#ACD4EC", "#6F99AD","#64AAD2", "#8AB6E9FF","#FE9586" ,"#D87B8B","#D898B9","#EBCC78", "#5DA59E" ,"#56BC9B","#9e6f8a", "#A3939D" ,  "#8C7F5F", "#C17E73", "#476D87")   ##E95C59",
subcelltypes_merged_need=c("ODC_1",   "ODC_2",   "ODC_3","ODC_4","ODC_5","Astro_1", "Astro_2","Astro_3","OPC",
"MG_1","MG_2","ExN","InN","Endo")
pal_subcelltypes_merged_need=c("#667FE1", "#ACD4EC", "#6F99AD","#64AAD2", "#8AB6E9FF","#FE9586" ,"#D87B8B","#D898B9","#EBCC78", "#5DA59E" ,"#56BC9B","#9e6f8a", "#A3939D" ,  "#8C7F5F")   ##E95C59",
SubCelltype_merge_scANVI_Celltype=c("ODC","Astro","OPC","MG","Neuron","Vascular")
pal_SubCelltype_merge_scANVI_Celltype=c("#8AB6E9FF","#D87B8B","#EBCC78", "#5DA59E" ,"#A48BCA", "#8C7F5F")
names(pal_SubCelltype_merge_scANVI_Celltype)=SubCelltype_merge_scANVI_Celltype
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)
groups= c("LOW", "HIGH")
compare_group1=c("LOW")
compare_group2=c("HIGH")
num_groups=length(groups)
num_com_groups=length(compare_group2)
pal_groups=c("#A7C9DF","#4880B8")
names(pal_groups)=groups

DEGS_psedobulk=file.path(PATH_O_data,"DEGs_psedobulk_edgeR",paste0("1.tsv"))
DEGS_psedobulk=read.csv(DEGS_psedobulk,sep="\t")
DEGS_psedobulk$celltype=DEGS_psedobulk$Celltype


p_cutoff=0.05
logFC_cutoff=0.2



############################################psedo-bulk based DEGs--> enrich########################################
DEGs_data=DEGS_psedobulk
DEGs_data$effect_size=DEGs_data$logFC
DEGs_data$compare_group="HIGH_vs_LOW"

for (celltype in subcelltypes_merged_need) {
   for (i in 1:length(compare_group1)) {
     group1 <- compare_group1[i]
     group2 <- compare_group2[i]
     comparison=paste0(group2,"_vs_",group1)
     mark=paste0(celltype,"_",comparison)
     data=DEGs_data[DEGs_data$celltype == celltype,   ]
     data=data[data$compare_group == comparison, ]
       data=data[,c("gene","effect_size")]
       rownames(data)=1:(dim(data)[1])
        geneList<- data$effect_size
        names(geneList)= data$gene
        geneList=sort(geneList,decreasing = T)

####### A: GO ANALYSIS
if (length(geneList) > 20) {
  cat(celltype, "have enough DEGs")
  gse_results <- gseGO(
    geneList      = geneList,
    OrgDb         = org.Hs.eg.db,
    ont           = "ALL",
    keyType       = "SYMBOL",
    nPerm         = 1000,
    minGSSize     = 2,
    maxGSSize     = 500,
    pvalueCutoff  = 0.1,
    pAdjustMethod = "none",
    verbose       = FALSE,
    seed          = 6,
    by            = "fgsea"
  )
}

require(DOSE)
show_term=min(10,dim(gse_results@result)[1])
 p=dotplot(gse_results,
  color = "NES",
  split=".sign",title=mark, showCategory=show_term)+facet_grid(~.sign)#点状图
 fig_filename=file.path(PATH_O_fig,"Enrichments_clusterProfiler/DEGs_psedobulk_edgeR",paste0("3_ALS_integratedALS_subcluster_DEG_edgeR_",mark,"_enrichemnt_clusterprofiler_goGSEA_dotplot_01052026.svg"))
 svglite(fig_filename)
 print(p)
 dev.off()
 data_file=file.path(PATH_O_data,"Enrichments_clusterProfiler/DEGs_psedobulk_edgeR",paste0("3_ALS_integratedALS_subcluster_DEG_edgeR_",mark,"_enrichemnt_clusterprofiler_goGSEA_01052026.rds"))
 saveRDS(gse_results,file=data_file)
 }
}

 cat(mark, "gsego enriched termes #(p.adj<0.1)",dim(gse_results@result)[1],"\n")
       gene.df <- bitr(data$gene,
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                 OrgDb = org.Hs.eg.db)
       gene.df = gene.df[!duplicated(gene.df[c("ENTREZID")]),]
       data2=merge(data,gene.df,by.x="gene",by.y="SYMBOL")
       geneList<-data2$effect_size
        names(geneList)=data2$ENTREZID
        geneList=sort(geneList,decreasing = T)
 if (length(geneList)>20) {
 gse_results <- gseKEGG(geneList, organism = "hsa",
 			  pvalueCutoff = 0.1,
 			  pAdjustMethod = "none",
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
 fig_filename=file.path(PATH_O_fig,"Enrichments_clusterProfiler/DEGs_psedobulk_edgeR",paste0("3_ALS_integratedALS_subcluster_DEG_edgeR_",mark,"_enrichemnt_clusterprofiler_keggGSEA_dotplot_01052026.svg"))
svglite(fig_filename)
print(p)
dev.off()
 gse_results = setReadable(gse_results, OrgDb = org.Hs.eg.db,keyType="ENTREZID")
 data_file=file.path(PATH_O_data,"Enrichments_clusterProfiler/DEGs_psedobulk_edgeR",paste0("3_ALS_integratedALS_subcluster_DEG_edgeR_",mark,"_enrichemnt_clusterprofiler_keggGSEA_01052026.rds"))
 saveRDS(gse_results,file=data_file)
 }
}
 cat(mark, "gsekegg enriched terms #(p<0.1)",dim(gse_results@result)[1],"\n")
   }
 }



#merge all enriched pathways together
gsea_enrich_kegg_all=data.frame()
p_cut=0.05
NES_cutoff=1
 for (celltype in subcelltypes_merged_need) {
   for (i in 1:length(compare_group1)) {
     group1 <- compare_group1[i]
     group2 <- compare_group2[i]
     comparison=paste0(group2,"_vs_",group1)
     mark=paste0(celltype,"_",comparison)
      data_file=file.path(PATH_O_data,"Enrichments_clusterProfiler/DEGs_psedobulk_edgeR",paste0("3_ALS_integratedALS_subcluster_DEG_edgeR_",mark,"_enrichemnt_clusterprofiler_keggGSEA_01052026.rds"))
      data=readRDS(data_file)
      data=data@result
      data$celltype=celltype
      data$comparison=comparison
 if (nrow(gsea_enrich_kegg_all) > 0) {
   gsea_enrich_kegg_all <- rbind(gsea_enrich_kegg_all, data)
} else {
  gsea_enrich_kegg_all <- data
}
cat(celltype,dim(gsea_enrich_kegg_all)[1],"\n")
}
}

gsea_enrich_kegg_all$abs_NES=abs(gsea_enrich_kegg_all$NES)
gsea_enrich_kegg_all$regulate=ifelse (gsea_enrich_kegg_all$NES>0,"UP","DOWN")
gsea_enrich_kegg_all=gsea_enrich_kegg_all[order(gsea_enrich_kegg_all$abs_NES,decreasing=TRUE),]

file_path_all=file.path(PATH_O_data,"Enrichments_clusterProfiler/DEGs_psedobulk_edgeR",paste0("3_ALS_integratedALS_subcluster_DEG_edgeR_enrichemnt_clusterprofiler_keggGSEA_01052026_ALL.tsv"))
write.table(gsea_enrich_kegg_all,file_path_all,sep="\t",quote=FALSE)

gsea_enrich_kegg_all_sig=gsea_enrich_kegg_all[gsea_enrich_kegg_all$pvalue<p_cut & abs(gsea_enrich_kegg_all$NES)>NES_cutoff,]
file_path_all=file.path(PATH_O_data,"Enrichments_clusterProfiler/DEGs_psedobulk_edgeR",paste0("3_ALS_integratedALS_subcluster_DEG_edgeR_enrichemnt_clusterprofiler_keggGSEA_01052026_ALL_sig.tsv"))
write.table(gsea_enrich_kegg_all_sig,file_path_all,sep="\t",quote=FALSE)


