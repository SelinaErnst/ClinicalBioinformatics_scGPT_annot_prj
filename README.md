|Folder | Comments|
|-----|---------|
|[experiment01-practice](experiment01-practice/README.md)| Practicing with relevant python packages|
|[experiment02-PBMC](experiment02-PBMC/README.md)| Preparing PBMC3k dataset for further analysis|
|[experiment03-scGPT](experiment03-scGPT/README.md)| scGPT repository and tutorials|
|[experiment04-annotation](experiment04-annotation/README.md)| annotation of datasets with scGPT|

# conda environments

The conda environments used for the analysis are found in `/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Selina/env`.  

- simple scanpy methods: [scanpy_env](experiment02-PBMC/environment.yml)
- scGPT analysis: scgpt_env or [scgpt_old_env](experiment03-scGPT/environment_old.yml)
- annotation using scGPT: scgpt_env or [annot_env](experiment04-annotation/environment.yml)

```bash
conda env update -f environment.yml --prune
```

# data

I work with scRNA-seq data. 

- `PBMC3k` : There is the dataset of 3k PBMCs from a Healthy Donor, used in a clustering tutorial on the [scanpy documentation](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering-2017.html) referencing Seurat's [guided clustering tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html).
- `ms` : Multiple Sclerosis dataset provided by [scGPT](https://github.com/bowang-lab/scGPT/tree/main/data)
- `scvi_pools` : dataset fusion-driven sarcoma