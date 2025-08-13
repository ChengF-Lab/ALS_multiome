library(ArchR)
library('Cairo')
addArchRGenome("hg38")
addArchRThreads(50)
addArchRLocking(locking = TRUE)
set.seed(123)
library(EnsDb.Hsapiens.v86)


input_path="xxx"
output_path="xxx"

input_arrows <- list.files(path = output_atac_arrow, pattern = "arrow")
names <- unlist(lapply(input_arrow , function(x){
  name <- unlist(strsplit(x, ".")[[1]])
}))

names_arrows <- sub("\\.arrow", "", input_arrows)
input_arrows <- paste(output_atac_arrow,"/",input_arrows,sep="")
names(input_arrows) <- names_arrows
allsamples = c("ALS1","ALS9","ALS11","ALS17","ALS3","ALS4","ALS8","ALS16")
input_arrows=input_arrows[allsamples]


projMulti1 <- ArchRProject(ArrowFiles = input_path,outputDirectory = output_path,copyArrows = TRUE)  #create archr object use ouw atac data
projMulti1 <- saveArchRProject(ArchRProj = projMulti1, outputDirectory = getOutputDirectory(projMulti1), overwrite = TRUE, load = TRUE)


input_rna <- list.files(path = cellranger_path, pattern = "h5$")
names <- unlist(lapply(input_rna , function(x){
  name <- unlist(strsplit(x, "[_]")[[1]])
  name <- name[1]
}))
input_rna <- unlist(lapply(input_rna , function(x){
  x <- file.path(cellranger_path, x)
}))
names(input_rna ) <- names
input_rnas=input_rna[allsamples]

seRNA <- import10xFeatureMatrix(
 input = input_rna,
  names = names(input_rna),
  strictMatch = TRUE,   #remove genes whose metadata is mis-matched across samples, you can set strictMatch = TRUE
  features = genes(EnsDb.Hsapiens.v86)   #all of the genes encoded on the mitochondrial genome are rescued
)


counts_RNA=seRNA@assays@data$data
all_RNA_genes=rownames(seRNA)
rows_to_keep <- which(rowSums(counts_RNA) > 10)
RNA_genes_keep <- all_RNA_genes[rows_to_keep]
seRNA_keep=seRNA[RNA_genes_keep,]
counts_RNA=seRNA_keep@assays@data$data
rows_to_keep <- which(rowSums(counts_RNA) > 10)
identical(rownames(seRNA_keep),RNA_genes_keep)




projMulti1 <- loadArchRProject(path =input_path)
cellsToKeep <- which(getCellNames(projMulti1) %in% colnames(seRNA_keep))
projMulti2 <- subsetArchRProject(ArchRProj = projMulti1, cells = getCellNames(projMulti1)[cellsToKeep], outputDirectory = output_multiome,force=TRUE)
projMulti2 <- addGeneExpressionMatrix(input = projMulti2, seRNA = seRNA_keep, scaleTo = 10000,strictMatch = TRUE, force = TRUE)
projMulti2 <- saveArchRProject(ArchRProj = projMulti2, outputDirectory = getOutputDirectory(projMulti2), overwrite = TRUE, load = TRUE)
saveRDS(seRNA,file.path(output_multiome_supp_data,"xxx.rds"))

