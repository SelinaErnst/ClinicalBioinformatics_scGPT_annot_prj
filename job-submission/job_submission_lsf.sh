#!/bin/bash
#BSUB -J test
#BSUB -R "rusage[mem=100GB]"
#BSUB -n 10
#BSUB -W 2:00
#BSUB -o out_%J.out
#BSUB -e out_%J.err

source ~/.bashrc
cd /omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/project/experiment01/scripts/
conda activate scanpy_env
python big_anndata.py