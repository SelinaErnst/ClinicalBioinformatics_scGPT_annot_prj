#!/bin/bash
#BSUB -J test_lsf
#BSUB -q long
#BSUB -n 10
#BSUB -R "rusage[mem=150GB] span[hosts=1]"
#BSUB -W 02:00
#BSUB -o test_lsf.out
source ~/.bashrc
cd /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project/experiment01-practice/scripts/
conda activate scanpy_env
python 05_test_job.py