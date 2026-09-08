#prepare transcriptomic data from multiome and additional data (snRNA-seq) for integration |

library(Seurat)
library(SeuratDisk)
library(DoubletFinder)
library(harmony)
library(tidyverse)
library(Matrix)
library(matrixStats)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(paletteer)
library(clustree)
library(tidydr)
set.seed(123)

file_path <- "xxx"
PATH_O <- file.path(file_path, "ALS_review_comments_output11052025/0_Output_Integration_snRNAseq_together12122025")
PATH_O_data <- file.path(PATH_O,"Data")
PATH_O_data_seq <- file.path(PATH_O_data,"data_seq")
PATH_O_fig <- file.path(PATH_O,"Figures")
PATH_O_fig_polish <- file.path(PATH_O_fig,"Figure_polished")


muliome_ALS_RNA_file="/home/doul2/beegfs/doul2/Work/ALS/Archr_multiome_output240501/Data_seurat_others/data_seq/6_5_ALS_multiome_RNA_expression_log1p_merge_seurat_count_DEGs_cal240920.rds"



Indivi_IDs=c("ALS07","ALS10","ALS12","ALS14")
Sample_IDs=c("ALS-7","ALS-10","ALS-12","ALS-14")
num_samp=dim(Sample_IDs)[1]


pal_celltypes=c("#8AB6E9FF","#D87B8B","#EBCC78", "#5DA59E" ,"#A48BCA", "#8C7F5F")
pal_celltypes_light = c("#ACD4EC","#FE9586","#B3E0A6",  "#FFCE54", "#AC92EC", "#C68D71", "#AAB2BD", "#FFC685")
pal_subcelltypes=c("#667FE1", "#ACD4EC", "#6F99AD","#64AAD2", "#8AB6E9FF","#FE9586" ,"#D87B8B","#D898B9","#EBCC78", "#5DA59E" ,"#56BC9B","#A48BCA", "#8C7F5F")
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)


###only perform label transfer for ODC, astrocytes, microglia, others are keep same as previously annotated results

so_merged = readRDS(file.path(PATH_O_data_seq,"2_2_ALS_multiomeRNA_merge_extraRNA_with_clustering_1225205.rds"))
#KEEP NEURONS SUBCELLTYPE AS InN and ExN
so_merged$SubCelltype<- ifelse(
  so_merged$SubCelltype %in% c("InN_LHX6","InN_VIP_RELN") ,   #keep as neurons
"InN",
so_merged$SubCelltype
)

so_merged$SubCelltype_scanvi <- ifelse(
  so_merged$SubCelltype %in% c("ODC","Astro","MG") & so_merged$data_source =="extra_snRNA-seq" ,   #keep as neurons
"Unknown",
so_merged$SubCelltype
)

so_merged$SubCelltype_scanvi <- ifelse(
  so_merged$SubCelltype_scanvi %in% c("Neuron") & so_merged$data_source =="multiome_new" ,   #have subcluster for neuron
"Unknown",
so_merged$SubCelltype_scanvi
)

meta=so_merged@meta.data
gene_matrix=so_merged@assays$RNA@counts
so_merged_all_nuclei_for_scanvi_input <- CreateSeuratObject(counts = gene_matrix,
                                 meta.data = meta,
                           )
SaveH5Seurat(so_merged_all_nuclei_for_scanvi_input, file.path(PATH_O_data_seq,"4_label_transfer_scANVI","4_1_ALS_multiomeRNA_merge_extraRNA_scanvi_input_12252025_all_nuclei_mapping_ODCAstrocytesMGNeurons.h5Seurat"))
Convert( file.path(PATH_O_data_seq,"4_label_transfer_scANVI","4_1_ALS_multiomeRNA_merge_extraRNA_scanvi_input_12252025_all_nuclei_mapping_ODCAstrocytesMGNeurons.h5Seurat"), dest = "h5ad")














































































