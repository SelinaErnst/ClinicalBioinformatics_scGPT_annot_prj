#!/bin/bash
#BSUB -J scGPT_pred_pbmc3k_quick
#BSUB -q long
#BSUB -n 20
#BSUB -R "rusage[mem=100GB] span[hosts=1]"
#BSUB -o /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project/experiment04-annotation/jobs/logs/scGPT_pred_pbmc3k_quick.out
#BSUB -e /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project/experiment04-annotation/jobs/logs/scGPT_pred_pbmc3k_quick.err
START_TIME=$(date +%s)
source ~/.bashrc
cd /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project/experiment04-annotation/scripts
module load GCC/14.1.0
conda activate scgpt_env
python 01_pbmc_annot_wo_fine_tuning.py
END_TIME=$(date +%s)
DURATION=$(($END_TIME - $START_TIME))
printf "\n"; printf "%.0s#" {1..58}; printf "\n"
printf "⏱️ Total Runtime: %02d:%02d:%02d\n" $(($DURATION/3600)) $(($DURATION%3600/60)) $(($DURATION%60))
printf "%.0s#" {1..58}; printf "\n\n"
