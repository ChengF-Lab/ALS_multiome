#!/bin/bash


celltypes=("ODC" "Astro" "OPC" "MG" "Neuron") # "Endo")
subcelltypes=("ODC_1" "ODC_2" "ODC_3" "ODC_4" "ODC_5" "Astro_1" "Astro_2" "Astro_3" "OPC" "MG_1" "MG_2" "Neuron" "Endo")
celltypes_Fib=("ODC_HIGH_Fibrinigen" "Astro_HIGH_Fibrinigen" "OPC_HIGH_Fibrinigen" "MG_HIGH_Fibrinigen" "Neuron_HIGH_Fibrinigen" "Endo_HIGH_Fibrinigen" \
"ODC_LOW_Fibrinigen" "Astro_LOW_Fibrinigen" "OPC_LOW_Fibrinigen" "MG_LOW_Fibrinigen" "Neuron_LOW_Fibrinigen" "Endo_LOW_Fibrinigen")

bed_path_LOW_Celltype="xxx/1_peak_bed_format/"
annot_path_LOW_Celltype="/xxx/2_peak_annotation/"


for celltype in "${celltypes[@]}"
do
     python  xxx/ldsc/ldsc.py \
       --h2 "${GWAS_sumstats_path}/4_ALS_GCST005647_Nicolas2018.sumstats.gz" \
       --ref-ld-chr "${annot_path_LOW_Celltype}/${celltype}/Add_ALL1_only_baselineLD_snp." \
       --w-ld-chr xxx/LDSC/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
      --frqfile-chr xxx/LDSC/1000G_Phase3_frq/1000G.EUR.QC. \
       --overlap-annot \
        --n-blocks 1000 \
        --print-coefficients \
       --out "${ldsc_path_peak_GWAS_LOW_Celltype}/${celltype}/6_ALS_GCST005647_Nicolas2018_LDSC_partitioned_heritability_by_ALS_peak_annot"
    echo "${celltype} 4_ALS_GCST005647_Nicolas2018 is done by ALS peak anno"
done


