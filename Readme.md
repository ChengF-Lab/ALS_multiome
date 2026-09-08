# 🧠 Single-nucleus multiome analyses of blood-brain barrier leakage, glial activation, angiogenesis and neuronal loss in Amyotrophic Lateral Sclerosis motor cortex

## 📖 Overview

This repository contains the analysis pipeline used to process, integrate, and analyze single-nucleus multiome (snRNA-seq + snATAC-seq) data from postmortem ALS and control motor cortex tissue, stratified by cortical atrophy/fibrinogen extravasation status. Analyses are organized into six sequential modules, in the order they are used in the manuscript:

| Order | Folder | Description                                                                            |
|---|---|----------------------------------------------------------------------------------------|
| 1 | [`0_multiome_process/`](#0_multiome_process) | Single-nucleus multiome (snRNA + snATAC) processing and comparison                     |
| 2 | [`1_snRNA_process/`](#1_snrna_process) | Standalone snRNA-seq processing for additional ALS patients)                           |
| 3 | [`2_Integration_ALS_transcriptomic_scANVI/`](#2_integration_als_transcriptomic_scanvi) | Integration of transcriptomic data from multiome + snRNA-seq datasets via scANVI       |
| 4 | [`3_Cell_cell_interaction/`](#3_cell_cell_interaction) | Cell–cell communication analysis (CellChat v2)                                         |
| 5 | [`4_Integration_comparison_with_public_data/`](#4_integration_comparison_with_public_data) | Integration/comparison with public ALS snRNA-seq datasets (GSE219281, and Syn51105515) |
| 6 | [`5_heritability_enrichment/`](#5_heritability_enrichment) | LDSC heritability enrichment analysis by implementing snATAC and GWAS data             |

---
## 🔧 Requirements

| Environment | Requirements                                                                                                                                                                                                                                 |
|---|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **R** | R 4.3.1; Seurat 4.3.0; SeuratObject 4.1.3; Signac 1.10.0; ArchR 1.0.3; harmony 0.1.1; ggplot2 3.4.2; SummarizedExperiment 1.30.2; SingleCellExperiment 1.22.0; clusterProfiler 4.8.1; rGREAT 2.5.7; CellChat 1.6.1; coloc 5.2.2; edgeR 5.2.2 |
| **Python** | Python; scvi-tools (scANVI); scanpy 1.9.5                                                                                                                                                                                                    |
| **External tools** | Cell Ranger ARC 2.0.0; LDSC 1.0.1; PLINK 1.90; bedtools 2.27.0                                                                                                                                                                               |


## 0_multiome_process

Processing of raw single-nucleus multiome (snRNA-seq + snATAC-seq) data generated with 10x Genomics Chromium and aligned using Cell Ranger ARC (v2.0.0), including:
- Cell Ranger ARC output import and QC filtering (nCount/nFeature thresholds per assay, doublet removal)
- Peak calling and chromatin accessibility quantification (ArchR)
- Gene activity score calculation
- RNA/ATAC joint dimensionality reduction and clustering
- Cell type annotation using canonical marker genes


## 1_snRNA_process

Processing of additional snRNA-seq data generated to increase statistical power for statistical analysis, including:
- QC filtering, normalization, and clustering (Seurat)
- Doublet detection and removal
- Cell type annotation using canonical markers
- Harmony batch integration across donors


## 2_Integration_ALS_transcriptomic_scANVI

Integration of multiome-derived and standalone snRNA-seq datasets using scANVI (scvi-tools) for batch correction and joint cell-type annotation, including:
- Reference/query label transfer across donors and sequencing batches
- Batch effect correction
- Joint UMAP embedding and cluster annotation
- Donor-aware pseudobulk differential expression (edgeR) between fibrinogen-high and fibrinogen-low groups

## 3_Cell_cell_interaction

Cell–cell communication analysis using CellChat to identify signaling interactions between glial, vascular, and neuronal populations (e.g., astrocyte–endothelial, oligodendrocyte–neuron), including:
- Ligand–receptor interaction inference per group (control, fibrinogen-low, fibrinogen-high)
- Differential interaction strength/number analysis between groups

## 4_Integration_comparison_with_public_data

Integration and comparative analysis with two publicly available ALS snRNA-seq datasets (GSE219281, Syn51105515) spanning control, sALS, and fALS samples, including:
- Cross-dataset integration and batch correction
- Differential expression comparison across ALS subtypes and cell types
- Marker gene/pathway overlap analysis against previously reported reactive astrocyte and endothelial signatures

## 5_heritability_enrichment

Heritability enrichment analysis linking cell-type-specific chromatin accessibility peaks (snATAC, from module 0) to ALS GWAS risk variants using stratified LD score regression (LDSC), including:
- Cell-type-specific peak calling (all ALS patients vs. fibrinogen-high subgroup)
- LDSC heritability partitioning across cell types
- Variant prioritization within cell-type-specific regulatory elements (e.g., rs631312 at the MOBP locus)
- Colocalization analysis (coloc) between GWAS signal and cell-type-specific accessible chromatin


## 📦 Data Availability

Raw sequencing data have been deposited in the NCBI Gene Expression Omnibus (GEO, https://www.ncbi.nlm.nih.gov/geo/) under accession code GSE301545 and are publicly available as of the date of publication. 

## ✉️ Contact

For questions, please contact doul2@ccf.org or doulijun777@gmail.com.