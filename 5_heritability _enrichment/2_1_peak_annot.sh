#!/bin/bash

celltypes=("ODC" "Astro" "OPC" "MG" "Neuron") # "Endo")
subcelltypes=("ODC_1" "ODC_2" "ODC_3" "ODC_4" "ODC_5" "Astro_1" "Astro_2" "Astro_3" "OPC" "MG_1" "MG_2" "Neuron" "Endo")
celltypes_Fib=("ODC_HIGH_Fibrinigen" "Astro_HIGH_Fibrinigen" "OPC_HIGH_Fibrinigen" "MG_HIGH_Fibrinigen" "Neuron_HIGH_Fibrinigen" "Endo_HIGH_Fibrinigen" \
"ODC_LOW_Fibrinigen" "Astro_LOW_Fibrinigen" "OPC_LOW_Fibrinigen" "MG_LOW_Fibrinigen" "Neuron_LOW_Fibrinigen" "Endo_LOW_Fibrinigen")

bed_path_LOW_Celltype="xxx/1_peak_bed_format/"
annot_path_LOW_Celltype="/xxx/2_peak_annotation/"

for celltype in "${celltypes[@]}"
do
    for i in {11..22}
    do
        python /xxx/make_annot.py \
            --bed-file "${bed_path_LOW_Celltype}/1_ALS_denovo_peaks_LOW_Celltype_${celltype}_chr${i}_big_max_peaks250210_hg19.bed" \
            --bimfile "xxx/1000G_EUR_Phase3_plink/1000G.EUR.QC.${i}.bim" \
            --annot-file "${annot_path_LOW_Celltype}/${celltype}/${i}.annot.gz"

        # Add an "ALL" column where all values are set to 1
        gunzip -c "${annot_path_LOW_Celltype}/${celltype}/${i}.annot.gz" \
        | awk 'BEGIN{OFS="\t"}{if(NR==1){print $0,"ALL"}else{print $0,1}}' \
        | gzip -c > "${annot_path_LOW_Celltype}/${celltype}/Add_ALL1_only_baselineLD_snp.${i}.annot.gz"
    done
done