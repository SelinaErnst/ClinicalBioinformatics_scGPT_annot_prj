#!/bin/bash
#BSUB -J downsample_20per_lsf
#BSUB -q long
#BSUB -n 1
#BSUB -R "rusage[mem=150GB] span[hosts=1]"
#BSUB -W 01:00
#BSUB -o /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project/experiment01-practice/jobs/logs/downsample_20per_lsf.out
#BSUB -e /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project/experiment01-practice/jobs/logs/downsample_20per_lsf.err
START_TIME=$(date +%s)
source ~/.bashrc
cd /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project/experiment01-practice/scripts
conda activate scanpy_env
python 03_downsample.py
END_TIME=$(date +%s)
DURATION=$(($END_TIME - $START_TIME))
printf "\n"; printf "%.0s#" {1..58}; printf "\n"
printf "⏱️ Total Runtime: %02d:%02d:%02d\n" $(($DURATION/3600)) $(($DURATION%3600/60)) $(($DURATION%60))
printf "%.0s#" {1..58}; printf "\n\n"
