#!/bin/bash
#BSUB -J scGPT_prediction
#BSUB -q long
#BSUB -n 10
#BSUB -R "rusage[mem=10GB] span[hosts=1]"
#BSUB -W 05:00
#BSUB -o /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project/experiment03-scGPT/jobs/logs/scGPT_prediction.out
#BSUB -e /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project/experiment03-scGPT/jobs/logs/scGPT_prediction.err
START_TIME=$(date +%s)
source ~/.bashrc
cd /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project/experiment03-scGPT/scripts
module load GCC/14.1.0
conda activate scgpt_env
python 05_prediction_on_fine_tuned_model.py
END_TIME=$(date +%s)
DURATION=$(($END_TIME - $START_TIME))
printf "\n"; printf "%.0s#" {1..58}; printf "\n"
printf "⏱️ Total Runtime: %02d:%02d:%02d\n" $(($DURATION/3600)) $(($DURATION%3600/60)) $(($DURATION%60))
printf "%.0s#" {1..58}; printf "\n\n"
