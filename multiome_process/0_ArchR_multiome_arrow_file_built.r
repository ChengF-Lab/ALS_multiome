library(ArchR)
library('Cairo')
addArchRGenome("hg38")
addArchRThreads(50)
addArchRLocking(locking = TRUE)
set.seed(123)


cellranger_path="xxx"
output_path="xxx"


input_atac <- list.files(path = cellranger_path, pattern = "gz$")
names <- unlist(lapply(input_atac , function(x){
  name <- unlist(strsplit(x, "[_]")[[1]])  #use first split of filenames as indivi name
  name <- name[1]
}))
input_atac <- unlist(lapply(input_atac, function(x){
  x <- file.path(cellranger_path, x)
}))
names(input_atac) <- names
allsamples = c("ALS1","ALS9","ALS11","ALS17","ALS3","ALS4","ALS8","ALS16")

input_atac=input_atac[allsamples]  #for each data with related data path
sample_num=length(allsamples)
ArrowFiles <- createArrowFiles(
  inputFiles = input_atac,
  sampleNames = names(input_atac),
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,  #A boolean value indicating whether to add a "Tile Matrix" to each ArrowFile. A Tile Matrix is a counts matrix that, instead of using peaks, uses a fixed-width sliding window of bins across the whole genome. This matrix can be used in many downstream ArchR operations.
  addGeneScoreMat = TRUE,   #A boolean value indicating whether to add a Gene-Score Matrix to each ArrowFile. A Gene-Score Matrix uses ATAC-seq signal proximal to the TSS to estimate gene activity.
  force=TRUE
)

ArrowFiles=paste0(allsamples,".arrow")
file.copy(from = ArrowFiles, to = output_path)



