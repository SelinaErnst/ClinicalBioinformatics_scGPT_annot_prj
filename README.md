|Folder | Comments|
|-----|---------|
|experiment01-practice| Practicing with relevant python packages|
|experiment02-PBMC| Preparing PBMC3k dataset for further analysis|
|experiment03-scGPT| scGPT repository and tutorials|
|experiment04-annotation| annotation of datasets with scGPT|


# conda environments

### simple scanpy methods

The conda environment used for PBMC was `scanpy_env` as found in `/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/env`.  

```bash
my_dir='/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina'
conda env update -f $my_dir/project/environment.yml --prune
```

### scGPT analysis

The conda environment used for scGPT analysis is called `scgpt_env` and is located at `/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/env`.  

# adata

I work with scRNA-seq data. 
