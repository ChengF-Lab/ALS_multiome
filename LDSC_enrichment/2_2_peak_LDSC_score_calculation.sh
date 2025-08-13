#!/bin/bash

celltypes=("ODC" "Astro" "OPC" "MG" "Neuron") # "Endo")
subcelltypes=("ODC_1" "ODC_2" "ODC_3" "ODC_4" "ODC_5" "Astro_1" "Astro_2" "Astro_3" "OPC" "MG_1" "MG_2" "Neuron" "Endo")
celltypes_Fib=("ODC_HIGH_Fibrinigen" "Astro_HIGH_Fibrinigen" "OPC_HIGH_Fibrinigen" "MG_HIGH_Fibrinigen" "Neuron_HIGH_Fibrinigen" "Endo_HIGH_Fibrinigen" \
"ODC_LOW_Fibrinigen" "Astro_LOW_Fibrinigen" "OPC_LOW_Fibrinigen" "MG_LOW_Fibrinigen" "Neuron_LOW_Fibrinigen" "Endo_LOW_Fibrinigen")

bed_path_LOW_Celltype="xxx/1_peak_bed_format/"
annot_path_LOW_Celltype="/xxx/2_peak_annotation/"



celltypes=("ODC" "Astro" "OPC")  # "MG" "Neuron") # "Endo")

export bed_path_LOW_Celltype annot_path_LOW_Celltype
for celltype in "${celltypes[@]}"; do
    export celltype  # Ensure celltype is also exported for access by parallel
    seq 1 22 | parallel -j 22 --env bed_path_LOW_Celltype,annot_path_LOW_Celltype,celltype --delay 0.2 \
    'python xxx/ldsc/ldsc.py --l2 \
    --bfile xxx/1000G_EUR_Phase3_plink/1000G.EUR.QC.{} \
    --ld-wind-cm 1 \
    --annot "${annot_path_LOW_Celltype}/${celltype}/Add_ALL1_only_baselineLD_snp.{}.annot.gz" \
    --thin-annot \
        --print-snps xxx/LDSC/baselineLD_v2.3_snp_list_split_by_dou/baselineLD_v2.3_snplist.{}.snp \
    --out "${annot_path_LOW_Celltype}/${celltype}/Add_ALL1_only_baselineLD_snp.{}"'
done




