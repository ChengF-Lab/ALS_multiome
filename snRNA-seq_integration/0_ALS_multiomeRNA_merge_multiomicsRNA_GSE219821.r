##merge RNA data with public data within motor cortex
library(Seurat)
library(harmony)
library(cowplot)
library(ggplot2)
library(svglite)



muliome_ALS_RNA_file="xxx.rds"
muliomics_ALS_RNA_file="xxx.rds"
out_path="xxx"

so1=readRDS(muliome_ALS_RNA_file)
so2=readRDS(muliomics_ALS_RNA_file)

meta1=so1@meta.data


meta2=so2@meta.data
meta=rbind(meta1,meta2)


gene_matrix1=so1@assays$RNA@counts
gene_matrix2=so2@assays$RNA@counts

overlapped_genes=intersect(rownames(gene_matrix1),rownames(gene_matrix2))
gene_matrix1=gene_matrix1[overlapped_genes,]
gene_matrix2=gene_matrix2[overlapped_genes,]
gene_matrix_merge=cbind(gene_matrix1,gene_matrix2)

meta_merge=rbind(meta1,meta2)
identical(rownames(meta_merge),colnames(gene_matrix_merge))

so_merge <- CreateSeuratObject(counts = gene_matrix_merge,
                                 meta.data = meta_merge,
                                )
saveRDS(so_merge,file.path(out_path,"xxx.rds"),compress=FALSE)




# Normalization and scaling
so_merge <- NormalizeData(so_merge)
so_merge <- FindVariableFeatures(so_merge)
so_merge <- ScaleData(so_merge)

dim.usage=30
so_merge <- RunPCA(so_merge,reduction.name = 'PCA', reduction.key = 'PCA_')
so_merge <- RunHarmony(so_merge, reduction="PCA",c("Indivi_ID","data_source"),reduction.save = "harmony_IndiviID_datasource")
so_merge <- RunUMAP(so_merge,  dims = 1:dim.usage,reduction = "harmony_IndiviID_datasource",reduction.name = 'UMAP_harmony_IndiviID_datasource')
so_merge <- FindNeighbors(so_merge,reduction="harmony_IndiviID_datasource",graph.name="harmony_IndiviID_datasource_knn",k.param = 30,  dims = 1:dim.usage)
so_merge <- FindClusters(so_merge, graph.name="harmony_IndiviID_datasource_knn",resolution = c(0.1))

saveRDS(so_merge,file.path(out_path,"xxx.rds"),compress=FALSE)


##QC metrics
so=readRDS(file.path(out_path_data_seq,"xxx.rds"))
mito.genes <- grep(pattern = "^MT", rownames(so))  #For human
mito.genes <- rownames(so)[grep(pattern = "^MT-", rownames(so),ignore.case = TRUE)]  #For mouse
rb.genes <- rownames(so)[grep(pattern = "^RP[SL]", rownames(so),ignore.case = TRUE)]  #For mouse
hb.genes <- rownames(so)[grep(pattern = "^HB[^(P)]", rownames(so),ignore.case = TRUE)]  #For mouse
so[["percent.mt"]] <- PercentageFeatureSet(so, features =mito.genes)
so[["percent.ribo"]] <- PercentageFeatureSet(so, features =rb.genes)
so[["percent.hb"]] <- PercentageFeatureSet(so, features =hb.genes)

P=VlnPlot(so,features=c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","percent.hb"),pt.size=0, group.by="Indivi_ID",raster=TRUE)
svglite(file.path(out_path,"0_ALS_RNA_merge_before_QC_Vlnplot_241115.svg"),width=15,height=5)
print(P)
dev.off()


##do QC for merged data
so=subset(so, subset=nCount_RNA<=3e4 & nFeature_RNA<=6e3)  #1e4 for our data; #4000 for our data
so=subset(so, subset=nCount_RNA>=200)  #1e4 for our data
so=subset(so, subset=nFeature_RNA>=200)   #4000 for our data

P=VlnPlot(so,features=c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","percent.hb"),pt.size=0, group.by="Indivi_ID",raster=TRUE)
svglite(file.path(out_path,"0_ALS_RNA_merge_after_QC_Vlnplot_241115.svg"),width=10,height=5)
print(P)
dev.off()
saveRDS(so,file.path(out_path,"xxx.rds"),compress=FALSE)


