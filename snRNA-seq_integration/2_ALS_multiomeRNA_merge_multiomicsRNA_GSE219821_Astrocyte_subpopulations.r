##merge RNA data with public data within motor cortex
library(Seurat)
library(harmony)
library(cowplot)
library(ggplot2)
library(svglite)
library(scCustomize)

#清空环境变量
rm(list = ls())
options(stringsAsFactors = F)
#加载相关的包
library(data.table)
library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)#加载所需物种基因组注释信息，可有20个物种选择，http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
library(org.Mm.eg.db)
library(ggplot2)
library(enrichplot)
library(ReactomePA)




in_path="xxx.rds"
out_path="xxx"


so=readRDS(file.path(out_path_data_seq,"xxx.rds"))
so$ALS_Fibrinogen=paste(so$ALS_diagnosis,so$Fibrinigen,sep="_")
so_astro=subset(so,subset=Celltype %in% c("Astro"))

metadata=so_astro@meta.data
counts=so_astro@assays$RNA@counts
so_astro <- CreateSeuratObject(counts = counts,  min.cells = 3, min.features = 200,meta.data=metadata)


##redo cluster
dim.usage=30
so_astro <- NormalizeData(so_astro)
so_astro <- FindVariableFeatures(so_astro, selection.method = "vst", nfeatures = 2000)
so_astro <- ScaleData(so_astro)
so_astro <- RunPCA(so_astro,reduction.name = 'astro_PCA', reduction.key = 'astro_PCA_')
so_astro <- RunHarmony(so_astro, reduction="astro_PCA",c("Indivi_ID","data_source"),reduction.save = "astro_harmony_indiviID_datasource")
so_astro <- FindNeighbors(so_astro,reduction="astro_harmony_indiviID_datasource",k.param = 30,  dims = 1:dim.usage,graph.name="astro_harmony_indiviID_datasource_knn")
so_astro <- RunUMAP(so_astro,  dims = 1:dim.usage,reduction = "astro_harmony_indiviID_datasource",reduction.name = 'UMAP_astro_harmony_indiviID_datasource')
so_astro <- FindClusters(so_astro, resolution = c(0.1,0.05), graph.name="astro_harmony_indiviID_datasource_knn",) #

saveRDS(so_astro,file.path(out_path,"xxx.rds"))


p2=FeaturePlot(so_astro,features=c("C3","AQP4","AQP1","CD44","GFAP"),raster=TRUE)
p21=DimPlot(so_astro,reduction="UMAP_astro_harmony_indiviID_datasource",group.by="astro_harmony_indiviID_datasource_knn_res.0.05",label=TRUE,raster=TRUE)
p=plot_grid(p2+p21, ncol = 1)
ggsave(p,filename = file.path(out_path,"xxx.svg"),width =6,height = 10,dpi=100)



so$Astro_subcelltype=so$astro_harmony_indiviID_datasource_knn_res.0.1
saveRDS(so,file.path(out_path,"xxx.rds"))
