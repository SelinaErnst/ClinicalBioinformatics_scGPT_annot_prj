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

### 1. Preprocessing

```python
adata.var["mt"] = adata.var_names.str.startswith("MT-")
scanpy.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], inplace=True, log1p=False

scanpy.pp.normalize_total(adata,target_sum=1e4)
scanpy.pp.log1p(adata)
scanpy.pp.highly_variable_genes(adata,flavor='seurat',n_bins=20,n_top_genes=None)
```

### 2. Dimensionionality Reduction

```python
scanpy.tl.pca(adata)
scanpy.pp.neighbors(adata, n_neighbors=15, use_rep='X_PCA')
scanpy.tl.umap(adata)
```

### 3. Clustering

```python
scanpy.tl.leiden(adata)
```

[single-cell best practices](https://www.sc-best-practices.org/preprocessing_visualization/normalization.html)
