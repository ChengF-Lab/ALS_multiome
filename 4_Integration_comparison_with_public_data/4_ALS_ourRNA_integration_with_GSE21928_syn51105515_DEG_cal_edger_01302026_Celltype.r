#calculate DEGs

import os
import sys
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import harmonypy as hm
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
import umap
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap


np.random.seed(1234)





file_path = "xxx"
PATH_O = os.path.join(file_path, "ALS_review_comments_output11052025/6_Integration_snRNA_with_public_GSE_01052026_syn51105515_results_01302026")
PATH_O_data = os.path.join(PATH_O, "Data")
PATH_O_data_seq = os.path.join(PATH_O_data, "data_seq")
PATH_O_fig = os.path.join(PATH_O, "Figures")
PATH_O_fig_polish = os.path.join(PATH_O_fig, "Figure_polished")
PATH_O_fig_polish_data = os.path.join(PATH_O_fig_polish, "Data")
PATH_syn51105515 = "xxx/0_Data_collection_integration_results/0_data_download_preparation_results/ALS_FTD/Syn51105515/Data/data_seq"
PATH_syn51105515_meta="xxx/0_Data_collection_integration_results/0_data_download_preparation_results/ALS_FTD/Syn51105515/Data/data_seq/Syn51105515_metadata_08152025.tsv"
Celltypes_updated_va=['ExN', 'ODC', 'Astro', 'OPC', 'MG', 'InN', 'T cells', 'Other vascular', 'Endo']
Celltypes_updated_va_need=['ExN', 'ODC', 'Astro', 'OPC', 'MG', 'InN', 'Endo']








Integrated_Cellclass = c("ODC", "Astro", "OPC", "MG", "ExN", "InN", "Endo", "Vascular","T cells")
Integrated_Cellclass_need = c("ODC", "Astro", "OPC", "MG", "ExN", "InN", "Endo", "Vascular")
num_Cellclass = length(Integrated_Cellclass_need)





####################################################################################################
####################################################################################################
####FibHigh_vs_CtrlGSE219281_syn51105515


Celltype_updated_vas = c("ODC", "Astro", "OPC", "MG", "ExN", "InN", "Endo")

Region_need="PMC"

for (celltype in Celltype_updated_vas){
celltype_safe <- gsub("/", "_", celltype)
celltype_safe <- gsub(" ", "_", celltype_safe)
pseudobulk_with_meta = file.path(PATH_O_data_DEGs,"Psedobulk",
                            paste0("3_integrated_ALS_DEGs_Step1_pseudobulk_counts_Celltype_updated_vas_",celltype_safe,"_01302026.tsv"))

pseudobulk_with_meta <- read.csv(pseudobulk_with_meta, sep="\t",row.names = 1)
pseudobulk_with_meta$Fibrinogenlevel_ALS_type=paste(pseudobulk_with_meta$Fibrinogen_level, pseudobulk_with_meta$ALS_type,sep="_")
pseudobulk_with_meta_need <- pseudobulk_with_meta[pseudobulk_with_meta$Region == Region_need, ]

pseudobulk_with_meta_need <- pseudobulk_with_meta_need[pseudobulk_with_meta_need$Fibrinogenlevel_ALS_type %in% c("Unknown_Control","HIGH_ALS"), ]
pseudobulk_with_meta_need$Fibrinogenlevel_ALS_type <-
  factor(pseudobulk_with_meta_need$Fibrinogenlevel_ALS_type,
         levels = c("Unknown_Control","HIGH_ALS"))  ##then we have AD vs. Normal

pseudobulk_with_meta_need <- pseudobulk_with_meta_need[,
  c("Fibrinogenlevel_ALS_type", setdiff(colnames(pseudobulk_with_meta_need), "Fibrinogenlevel_ALS_type"))
]


pseudobulk_meta=pseudobulk_with_meta_need[,1:18]
pseudobulk_count=pseudobulk_with_meta_need[,19:ncol(pseudobulk_with_meta_need)]

##remove samples with NA
covariates_to_check <- c("Age", "Sex", "PMI",
                           "Fibrinogenlevel_ALS_type","ALS_group_sf","data_source")

# 1) 仅报告这三个协变量的缺失情况（不删除）
na_counts <- sapply(pseudobulk_meta[, covariates_to_check, drop=FALSE], function(x) sum(is.na(x)))
cat("\nMissing values (Age/Sex/PMI):\n")
print(na_counts)

auto_fix_covariates_age_sex_pmi <- function(meta, age_col="Age", sex_col="Sex", pmi_col="PMI") {
  # --- AGE ---
  use_age_missing <- any(is.na(meta[[age_col]]))
  if (use_age_missing) {
    meta$Age_missing <- is.na(meta[[age_col]])
    meta$Age_filled  <- meta[[age_col]]
    meta$Age_filled[is.na(meta$Age_filled)] <- median(meta$Age_filled, na.rm=TRUE)
    age_terms <- "Age_filled + Age_missing"
  } else {
    age_terms <- age_col
  }
  # --- PMI ---
  use_pmi_missing <- any(is.na(meta[[pmi_col]]))
  if (use_pmi_missing) {
    meta$PMI_missing <- is.na(meta[[pmi_col]])
    meta$PMI_filled  <- meta[[pmi_col]]
    meta$PMI_filled[is.na(meta$PMI_filled)] <- median(meta$PMI_filled, na.rm=TRUE)
    pmi_terms <- "PMI_filled + PMI_missing"
  } else {
    pmi_terms <- pmi_col
  }
  # --- SEX ---
  # keep as factor; if missing, set to "Unknown"
  meta[[sex_col]] <- as.character(meta[[sex_col]])
  if (any(is.na(meta[[sex_col]]) | meta[[sex_col]] == "")) {
    meta[[sex_col]][is.na(meta[[sex_col]]) | meta[[sex_col]] == ""] <- "Unknown"
  }
  meta[[sex_col]] <- factor(meta[[sex_col]])

  # report
  cat("\n[Auto-check] Missingness summary:\n")
  cat("  Age missing:", sum(is.na(meta[[age_col]])), "\n")
  cat("  PMI missing:", sum(is.na(meta[[pmi_col]])), "\n")
  cat("  Sex missing/blank fixed to 'Unknown':",
      sum(meta[[sex_col]] == "Unknown", na.rm=TRUE), "\n")

  list(meta=meta, age_terms=age_terms, pmi_terms=pmi_terms, sex_col=sex_col)
}


res <- auto_fix_covariates_age_sex_pmi(pseudobulk_meta)
pseudobulk_meta <- res$meta
age_terms_for_edger=res$age_terms
pmi_terms_for_edger=res$pmi_terms
sex_col_for_edger=res$sex_col
stopifnot(identical(rownames(pseudobulk_meta), rownames(pseudobulk_count)))  #
counts_clean <- t(pseudobulk_count)

dge <- DGEList(counts = counts_clean)

keep <- filterByExpr(dge, group = pseudobulk_meta$Fibrinogenlevel_ALS_type)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

form_str <- paste0(
  "~ factor(Fibrinogenlevel_ALS_type) + ",
  age_terms_for_edger, " + ",
  sex_col_for_edger, " + ",
  pmi_terms_for_edger
)


design <- model.matrix(as.formula(form_str), data = pseudobulk_meta)



# Fit model
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)


qlf <- glmQLFTest(fit, coef = "factor(Fibrinogenlevel_ALS_type)HIGH_ALS")  # change to one you need

degs <- topTags(qlf, n = Inf)$table
degs$Region=Region_need    ####
degs$celltype=celltype      ####
degs$comparison = "FibHigh vs. Ctrl(GSE219281_syn51105515)"   ####
degs$gene = rownames(degs)   ####
write.table(degs, file = file.path(PATH_O_data_DEGs,"Edger_results","FibHigh_vs_CtrlGSE219281_syn51105515",
        paste0("4_ALS_DEGs_FibHigh_vs_CtrlGSE219281_syn51105515_", Region_need,"_", celltype_safe, "_edgeR_covariates_01302026.tsv")),sep="\t",quote=FALSE)



PValue_cutoff=0.5
FDR_cutoff=0.5
logFC_cutoff=0.15


df <- degs
df <- df %>%
  mutate(
    sig = case_when(
      PValue < PValue_cutoff & logFC >  logFC_cutoff ~ "Up",
      PValue < PValue_cutoff & logFC < -logFC_cutoff ~ "Down",
      TRUE ~ "NS"
    )
  )

df_top <- df %>%
  filter(sig != "NS") %>%
  arrange(PValue) %>%
  slice_head(n = 10)

p <- ggplot(df, aes(x = logFC, y = -log10(PValue))) +
  geom_point(aes(color = sig), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("Up"="#d73027", "Down"="#4575b4", "NS"="grey70")) +
  geom_hline(yintercept = -log10(PValue_cutoff), linetype="dashed") +
  geom_text_repel(
    data = df_top,
    aes(label = gene),
    size = 3,
    max.overlaps = 20
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = paste0(celltype," log(Fold Change)"),
    y = "-log10(PValue)",
    color = "Regulation"
  )
ggsave(
  filename = file.path(PATH_O_fig_DEGs,"Edger_results","FibHigh_vs_CtrlGSE219281_syn51105515",paste0("4_ALS_DEGs_FibHigh_vs_CtrlGSE219281_syn51105515_",Region_need,"_", celltype_safe, "_edgeR_covariates_Age_sex_pmi_12102025.png")),
  plot = p, width = 6.5, height = 5, dpi = 300, bg = "white"
)

cat(celltype,"is done")

png(file.path(PATH_O_fig_DEGs,"Edger_results","FibHigh_vs_CtrlGSE219281_syn51105515",paste0("4_ALS_DEGs_FibHigh_vs_CtrlGSE219281_syn51105515_",Region_need,"_", celltype_safe, "_edgeR_covariates_Age_sex_pmi_12102025_MA.png")), width=1800, height=1500, res=300)
plotMD(qlf, main="edgeR MA plot (plotMD)")
abline(h=c(-1,1), col="blue", lty=2)
dev.off()


png(
  file.path(PATH_O_fig_DEGs,"Edger_results","FibHigh_vs_CtrlGSE219281_syn51105515",
            paste0("4_ALS_DEGs_FibHigh_vs_CtrlGSE219281_syn51105515_",
                   Region_need,"_", celltype_safe, "_edgeR_covariates_Age_sex_pmi_12102025_MDS.png")),
  width=1800, height=1500, res=300
)
lcpm <- cpm(dge, log=TRUE, prior.count=2)
src <- factor(pseudobulk_meta$data_source)
grp <- factor(pseudobulk_meta$Fibrinogenlevel_ALS_type)
pchv <- ifelse(grp=="HIGH_ALS", 17, 16)
plotMDS(lcpm, col=as.numeric(src), pch=pchv, main="MDS: color=source, shape=group")
legend("topleft", legend=levels(src), col=seq_along(levels(src)), pch=16, bty="n")
legend("topright", legend=levels(grp), pch=c(16,17), bty="n")
dev.off()
}


####################################################################################################
####################################################################################################
####C9ALS_vs_CtrlGSE219281
Celltype_updated_vas = c("ODC", "Astro", "OPC", "MG", "ExN", "InN", "Endo")
Region_need="PMC"
for (celltype in Celltype_updated_vas){
celltype_safe <- gsub("/", "_", celltype)
celltype_safe <- gsub(" ", "_", celltype_safe)
pseudobulk_with_meta = file.path(PATH_O_data_DEGs,"Psedobulk",
                            paste0("3_integrated_ALS_DEGs_Step1_pseudobulk_counts_Celltype_updated_vas_",celltype_safe,"_01302026.tsv"))
pseudobulk_with_meta <- read.csv(pseudobulk_with_meta, sep="\t",row.names = 1)
pseudobulk_with_meta$Fibrinogenlevel_ALS_type=paste(pseudobulk_with_meta$Fibrinogen_level, pseudobulk_with_meta$ALS_type,sep="_")
pseudobulk_with_meta_need <- pseudobulk_with_meta[pseudobulk_with_meta$Region == Region_need & pseudobulk_with_meta$data_source =="GSE219281", ]

pseudobulk_with_meta_need <- pseudobulk_with_meta_need[pseudobulk_with_meta_need$Fibrinogenlevel_ALS_type %in% c("Unknown_Control","Unknown_C9-ALS"), ]
pseudobulk_with_meta_need$Fibrinogenlevel_ALS_type <-
  factor(pseudobulk_with_meta_need$Fibrinogenlevel_ALS_type,
         levels = c("Unknown_Control","Unknown_C9-ALS"))

pseudobulk_with_meta_need <- pseudobulk_with_meta_need[,
  c("Fibrinogenlevel_ALS_type", setdiff(colnames(pseudobulk_with_meta_need), "Fibrinogenlevel_ALS_type"))
]

pseudobulk_meta=pseudobulk_with_meta_need[,1:18]
pseudobulk_count=pseudobulk_with_meta_need[,19:ncol(pseudobulk_with_meta_need)]

covariates_to_check <- c("Age", "Sex", "PMI",
                           "Fibrinogenlevel_ALS_type","ALS_group_sf","data_source")
na_counts <- sapply(pseudobulk_meta[, covariates_to_check, drop=FALSE], function(x) sum(is.na(x)))
cat("\nMissing values (Age/Sex/PMI):\n")
print(na_counts)

auto_fix_covariates_age_sex_pmi <- function(meta, age_col="Age", sex_col="Sex", pmi_col="PMI") {
  use_age_missing <- any(is.na(meta[[age_col]]))
  if (use_age_missing) {
    meta$Age_missing <- is.na(meta[[age_col]])
    meta$Age_filled  <- meta[[age_col]]
    meta$Age_filled[is.na(meta$Age_filled)] <- median(meta$Age_filled, na.rm=TRUE)
    age_terms <- "Age_filled + Age_missing"
  } else {
    age_terms <- age_col
  }
  use_pmi_missing <- any(is.na(meta[[pmi_col]]))
  if (use_pmi_missing) {
    meta$PMI_missing <- is.na(meta[[pmi_col]])
    meta$PMI_filled  <- meta[[pmi_col]]
    meta$PMI_filled[is.na(meta$PMI_filled)] <- median(meta$PMI_filled, na.rm=TRUE)
    pmi_terms <- "PMI_filled + PMI_missing"
  } else {
    pmi_terms <- pmi_col
  }
  meta[[sex_col]] <- as.character(meta[[sex_col]])
  if (any(is.na(meta[[sex_col]]) | meta[[sex_col]] == "")) {
    meta[[sex_col]][is.na(meta[[sex_col]]) | meta[[sex_col]] == ""] <- "Unknown"
  }
  meta[[sex_col]] <- factor(meta[[sex_col]])
  cat("\n[Auto-check] Missingness summary:\n")
  cat("  Age missing:", sum(is.na(meta[[age_col]])), "\n")
  cat("  PMI missing:", sum(is.na(meta[[pmi_col]])), "\n")
  cat("  Sex missing/blank fixed to 'Unknown':",
      sum(meta[[sex_col]] == "Unknown", na.rm=TRUE), "\n")
  list(meta=meta, age_terms=age_terms, pmi_terms=pmi_terms, sex_col=sex_col)
}


res <- auto_fix_covariates_age_sex_pmi(pseudobulk_meta)
pseudobulk_meta <- res$meta
age_terms_for_edger=res$age_terms
pmi_terms_for_edger=res$pmi_terms
sex_col_for_edger=res$sex_col
stopifnot(identical(rownames(pseudobulk_meta), rownames(pseudobulk_count)))  #

counts_clean <- t(pseudobulk_count)
dge <- DGEList(counts = counts_clean)
keep <- filterByExpr(dge, group = pseudobulk_meta$Fibrinogenlevel_ALS_type)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

form_str <- paste0(
  "~ factor(Fibrinogenlevel_ALS_type) + ",
  age_terms_for_edger, " + ",
  sex_col_for_edger, " + ",
  pmi_terms_for_edger
)


design <- model.matrix(as.formula(form_str), data = pseudobulk_meta)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = "factor(Fibrinogenlevel_ALS_type)Unknown_C9-ALS")  # change to one you need

degs <- topTags(qlf, n = Inf)$table
degs$Region=Region_need
degs$celltype=celltype
degs$comparison = "C9ALS vs. Ctrl(GSE219281)"
degs$gene = rownames(degs)
write.table(degs, file = file.path(PATH_O_data_DEGs,"Edger_results","C9ALS_vs_CtrlGSE2192815",
        paste0("4_ALS_DEGs_C9ALS_vs_CtrlGSE219281_", Region_need,"_", celltype_safe, "_edgeR_covariates_01302026.tsv")),sep="\t",quote=FALSE)

PValue_cutoff=0.5
FDR_cutoff=0.5
logFC_cutoff=0.15

df <- degs
df <- df %>%
  mutate(
    sig = case_when(
      PValue < PValue_cutoff & logFC >  logFC_cutoff ~ "Up",
      PValue < PValue_cutoff & logFC < -logFC_cutoff ~ "Down",
      TRUE ~ "NS"
    )
  )

df_top <- df %>%
  filter(sig != "NS") %>%
  arrange(PValue) %>%
  slice_head(n = 10)
p <- ggplot(df, aes(x = logFC, y = -log10(PValue))) +
  geom_point(aes(color = sig), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("Up"="#d73027", "Down"="#4575b4", "NS"="grey70")) +
  geom_hline(yintercept = -log10(PValue_cutoff), linetype="dashed") +
  geom_text_repel(
    data = df_top,
    aes(label = gene),
    size = 3,
    max.overlaps = 20
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = paste0(celltype," log(Fold Change)"),
    y = "-log10(PValue)",
    color = "Regulation"
  )

ggsave(
  filename = file.path(PATH_O_fig_DEGs,"Edger_results","C9ALS_vs_CtrlGSE2192815",paste0("4_ALS_DEGs_C9ALS_vs_CtrlGSE219281_",Region_need,"_", celltype_safe, "_edgeR_covariates_Age_sex_pmi_12102025.png")),
  plot = p, width = 6.5, height = 5, dpi = 300, bg = "white"
)
cat(celltype,"is done")

png(file.path(PATH_O_fig_DEGs,"Edger_results","C9ALS_vs_CtrlGSE2192815",paste0("4_ALS_DEGs_C9ALS_vs_CtrlGSE219281_",Region_need,"_", celltype_safe, "_edgeR_covariates_Age_sex_pmi_12102025_MA.png")), width=1800, height=1500, res=300)
plotMD(qlf, main="edgeR MA plot (plotMD)")
abline(h=c(-1,1), col="blue", lty=2)
dev.off()

png(
  file.path(PATH_O_fig_DEGs,"Edger_results","C9ALS_vs_CtrlGSE2192815",
            paste0("4_ALS_DEGs_C9ALS_vs_CtrlGSE219281_",
                   Region_need,"_", celltype_safe, "_edgeR_covariates_Age_sex_pmi_12102025_MDS.png")),
  width=1800, height=1500, res=300
)
lcpm <- cpm(dge, log=TRUE, prior.count=2)
src <- factor(pseudobulk_meta$data_source)
grp <- factor(pseudobulk_meta$Fibrinogenlevel_ALS_type)
pchv <- ifelse(grp=="Unknown_C9-ALS", 17, 16)
plotMDS(lcpm, col=as.numeric(src), pch=pchv, main="MDS: color=source, shape=group")
legend("topleft", legend=levels(src), col=seq_along(levels(src)), pch=16, bty="n")
legend("topright", legend=levels(grp), pch=c(16,17), bty="n")
dev.off()
}



####################################################################################################
####################################################################################################
####SALS_vs_CtrlSYN51105515
Celltype_updated_vas = c("ODC", "Astro", "OPC", "MG", "ExN", "InN", "Endo")
Region_need="PMC"
for (celltype in Celltype_updated_vas){
celltype_safe <- gsub("/", "_", celltype)
celltype_safe <- gsub(" ", "_", celltype_safe)
pseudobulk_with_meta = file.path(PATH_O_data_DEGs,"Psedobulk",
                            paste0("3_integrated_ALS_DEGs_Step1_pseudobulk_counts_Celltype_updated_vas_",celltype_safe,"_01302026.tsv"))

pseudobulk_with_meta <- read.csv(pseudobulk_with_meta, sep="\t",row.names = 1)
pseudobulk_with_meta$Fibrinogenlevel_ALS_type=paste(pseudobulk_with_meta$Fibrinogen_level, pseudobulk_with_meta$ALS_type,sep="_")
pseudobulk_with_meta_need <- pseudobulk_with_meta[pseudobulk_with_meta$Region == Region_need & pseudobulk_with_meta$data_source =="Syn51105515", ]

pseudobulk_with_meta_need <- pseudobulk_with_meta_need[pseudobulk_with_meta_need$Fibrinogenlevel_ALS_type %in% c("Unknown_Control","Unknown_ALS"), ]
pseudobulk_with_meta_need$Fibrinogenlevel_ALS_type <-
  factor(pseudobulk_with_meta_need$Fibrinogenlevel_ALS_type,
         levels = c("Unknown_Control","Unknown_ALS"))  ##then we have AD vs. Normal
pseudobulk_with_meta_need <- pseudobulk_with_meta_need[,
  c("Fibrinogenlevel_ALS_type", setdiff(colnames(pseudobulk_with_meta_need), "Fibrinogenlevel_ALS_type"))
]

pseudobulk_meta=pseudobulk_with_meta_need[,1:18]
pseudobulk_count=pseudobulk_with_meta_need[,19:ncol(pseudobulk_with_meta_need)]

covariates_to_check <- c("Age", "Sex", "PMI",
                           "Fibrinogenlevel_ALS_type","ALS_group_sf","data_source")
na_counts <- sapply(pseudobulk_meta[, covariates_to_check, drop=FALSE], function(x) sum(is.na(x)))
cat("\nMissing values (Age/Sex/PMI):\n")
print(na_counts)


auto_fix_covariates_age_sex_pmi <- function(meta, age_col="Age", sex_col="Sex", pmi_col="PMI") {
  use_age_missing <- any(is.na(meta[[age_col]]))
  if (use_age_missing) {
    meta$Age_missing <- is.na(meta[[age_col]])
    meta$Age_filled  <- meta[[age_col]]
    meta$Age_filled[is.na(meta$Age_filled)] <- median(meta$Age_filled, na.rm=TRUE)
    age_terms <- "Age_filled + Age_missing"
  } else {
    age_terms <- age_col
  }
  use_pmi_missing <- any(is.na(meta[[pmi_col]]))
  if (use_pmi_missing) {
    meta$PMI_missing <- is.na(meta[[pmi_col]])
    meta$PMI_filled  <- meta[[pmi_col]]
    meta$PMI_filled[is.na(meta$PMI_filled)] <- median(meta$PMI_filled, na.rm=TRUE)
    pmi_terms <- "PMI_filled + PMI_missing"
  } else {
    pmi_terms <- pmi_col
  }
  meta[[sex_col]] <- as.character(meta[[sex_col]])
  if (any(is.na(meta[[sex_col]]) | meta[[sex_col]] == "")) {
    meta[[sex_col]][is.na(meta[[sex_col]]) | meta[[sex_col]] == ""] <- "Unknown"
  }
  meta[[sex_col]] <- factor(meta[[sex_col]])
  cat("\n[Auto-check] Missingness summary:\n")
  cat("  Age missing:", sum(is.na(meta[[age_col]])), "\n")
  cat("  PMI missing:", sum(is.na(meta[[pmi_col]])), "\n")
  cat("  Sex missing/blank fixed to 'Unknown':",
      sum(meta[[sex_col]] == "Unknown", na.rm=TRUE), "\n")
  list(meta=meta, age_terms=age_terms, pmi_terms=pmi_terms, sex_col=sex_col)
}


res <- auto_fix_covariates_age_sex_pmi(pseudobulk_meta)
pseudobulk_meta <- res$meta
age_terms_for_edger=res$age_terms
pmi_terms_for_edger=res$pmi_terms
sex_col_for_edger=res$sex_col
stopifnot(identical(rownames(pseudobulk_meta), rownames(pseudobulk_count)))  #

counts_clean <- t(pseudobulk_count)
dge <- DGEList(counts = counts_clean)
keep <- filterByExpr(dge, group = pseudobulk_meta$Fibrinogenlevel_ALS_type)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
form_str <- paste0(
  "~ factor(Fibrinogenlevel_ALS_type) + ",
  age_terms_for_edger, " + ",
  sex_col_for_edger, " + ",
  pmi_terms_for_edger
)


design <- model.matrix(as.formula(form_str), data = pseudobulk_meta)



dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

qlf <- glmQLFTest(fit, coef = "factor(Fibrinogenlevel_ALS_type)Unknown_ALS")  # change to one you need

degs <- topTags(qlf, n = Inf)$table
degs$Region=Region_need
degs$celltype=celltype
degs$comparison = "sALS vs. Ctrl(Syn51105515)"
degs$gene = rownames(degs)
write.table(degs, file = file.path(PATH_O_data_DEGs,"Edger_results","sALS_vs_CtrlSYN51105515",
        paste0("4_ALS_DEGs_sALS_vs_CtrlSYN51105515_", Region_need,"_", celltype_safe, "_edgeR_covariates_01302026.tsv")),sep="\t",quote=FALSE)


PValue_cutoff=0.5
logFC_cutoff=0.15
df <- degs
df <- df %>%
  mutate(
    sig = case_when(
      PValue < PValue_cutoff & logFC >  logFC_cutoff ~ "Up",
      PValue < PValue_cutoff & logFC < -logFC_cutoff ~ "Down",
      TRUE ~ "NS"
    )
  )

df_top <- df %>%
  filter(sig != "NS") %>%
  arrange(PValue) %>%
  slice_head(n = 10)

p <- ggplot(df, aes(x = logFC, y = -log10(PValue))) +
  geom_point(aes(color = sig), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("Up"="#d73027", "Down"="#4575b4", "NS"="grey70")) +

  geom_hline(yintercept = -log10(PValue_cutoff), linetype="dashed") +
  geom_text_repel(
    data = df_top,
    aes(label = gene),
    size = 3,
    max.overlaps = 20
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = paste0(celltype," log(Fold Change)"),
    y = "-log10(PValue)",
    color = "Regulation"
  )

ggsave(
  filename = file.path(PATH_O_fig_DEGs,"Edger_results","sALS_vs_CtrlSYN51105515",paste0("4_ALS_DEGs_sALS_vs_CtrlSYN51105515_",Region_need,"_", celltype_safe, "_edgeR_covariates_Age_sex_pmi_12102025.png")),
  plot = p, width = 6.5, height = 5, dpi = 300, bg = "white"
)

cat(celltype,"is done")


png(file.path(PATH_O_fig_DEGs,"Edger_results","sALS_vs_CtrlSYN51105515",paste0("4_ALS_DEGs_sALS_vs_CtrlSYN51105515_",Region_need,"_", celltype_safe, "_edgeR_covariates_Age_sex_pmi_12102025_MA.png")), width=1800, height=1500, res=300)
plotMD(qlf, main="edgeR MA plot (plotMD)")
abline(h=c(-1,1), col="blue", lty=2)
dev.off()



png(
  file.path(PATH_O_fig_DEGs,"Edger_results","sALS_vs_CtrlSYN51105515",
            paste0("4_ALS_DEGs_sALS_vs_CtrlSYN51105515_",
                   Region_need,"_", celltype_safe, "_edgeR_covariates_Age_sex_pmi_12102025_MDS.png")),
  width=1800, height=1500, res=300
)
lcpm <- cpm(dge, log=TRUE, prior.count=2)
src <- factor(pseudobulk_meta$data_source)
grp <- factor(pseudobulk_meta$Fibrinogenlevel_ALS_type)
pchv <- ifelse(grp=="Unknown_C9-ALS", 17, 16)
plotMDS(lcpm, col=as.numeric(src), pch=pchv, main="MDS: color=source, shape=group")
legend("topleft", legend=levels(src), col=seq_along(levels(src)), pch=16, bty="n")
legend("topright", legend=levels(grp), pch=c(16,17), bty="n")
dev.off()
}









