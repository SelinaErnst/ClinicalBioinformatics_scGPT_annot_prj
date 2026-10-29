#!/bin/bash
#BSUB -J scGPT_prediction_only
#BSUB -q verylong
#BSUB -n 20
#BSUB -R "rusage[mem=250GB]"
#BSUB -o /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project/experiment04-annotation/jobs/logs/scGPT_prediction_only.out
#BSUB -e /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project/experiment04-annotation/jobs/logs/scGPT_prediction_only.err

START_TIME=$(date +%s)

source ~/.bashrc
cd /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project/experiment04-annotation/scripts
module load GCC/14.1.0
conda activate scgpt_env
python 03_annot_using_fine_tuned.py

END_TIME=$(date +%s)
DURATION=$(($END_TIME - $START_TIME))

printf "\n"; printf "%.0s#" {1..58}; printf "\n"
printf "⏱️ Total Runtime: %02d:%02d:%02d\n" $(($DURATION/3600)) $(($DURATION%3600/60)) $(($DURATION%60))
printf "%.0s#" {1..58}; printf "\n\n"
