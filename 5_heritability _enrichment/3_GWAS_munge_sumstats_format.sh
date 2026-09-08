#!/bin/bash
conda activate ldsc


celltypes=("ODC" "Astro" "OPC" "MG" "Neuron") # "Endo")
subcelltypes=("ODC_1" "ODC_2" "ODC_3" "ODC_4" "ODC_5" "Astro_1" "Astro_2" "Astro_3" "OPC" "MG_1" "MG_2" "Neuron" "Endo")
celltypes_Fib=("ODC_HIGH_Fibrinigen" "Astro_HIGH_Fibrinigen" "OPC_HIGH_Fibrinigen" "MG_HIGH_Fibrinigen" "Neuron_HIGH_Fibrinigen" "Endo_HIGH_Fibrinigen" \
"ODC_LOW_Fibrinigen" "Astro_LOW_Fibrinigen" "OPC_LOW_Fibrinigen" "MG_LOW_Fibrinigen" "Neuron_LOW_Fibrinigen" "Endo_LOW_Fibrinigen")

bed_path_LOW_Celltype="xxx/1_peak_bed_format/"
annot_path_LOW_Celltype="/xxx/2_peak_annotation/"

snplist=xxx/LDSC/w_hm3.snplist     #must have columns SNP, A1, A2.' (so this one can't be used)
xxx/ldsc/munge_sumstats.py \
    --sumstats "${GWAS_data_path}/ALS_GCST005647_Nicolas2018_hg19_and_hg38.tsv" \
    --merge-alleles $snplist \
    --snp rsid \
    --a1 effect_allele \
    --a2 other_allele \
    --N 80610 \
    --p P \
    --chunksize 500000 \
    --out "${munge_sumstats_path}/4_ALS_GCST005647_Nicolas2018"

echo "4_ALS_GCST005647_Nicolas2018 is done"