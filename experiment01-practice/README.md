# Practice with python packages

## Scripts
|Script | Comments|
|-----|---------|
|00_test.ipynb| test out packages and conda env|
|01_practice_anndata.ipynb| test out example adata file (SySa) |
|02_practice_scanpy.ipynb| scanpy legacy workflow (processing and clustering 3k PBMCs) |
|03_practice_scvi.ipynb| |
|04_downsample.py| downsample the big adata file for easier use |

## adata files
- big adata: `data/scvi_pools_2026_celltype_annotated.h5ad`.  
    - 481,411 (cells/barcodes) × 35,929 (genes/transcripts). 
    - Entity: RMA, EWING, URCS, SYSA, MYE, SFT, SEBA, ASPS, MLS, LGFMS, PLRMS, MYCHS, DDLS, LMS, PPMS, STET, cFS, DSRCT, TGCT, SARCNOS
    - Fusion: PAX7::FOXO1, PAX3::FOXO1, NCOA1::PAX3, FUS::ERG, YWHAE::NUTM2A, EWSR1::FLI1, EWSR1::ETV1, SS18::SSX1, EWSR1::POU5F1, NAB2::STAT6, SS18::SSX2, ASPSCR1::TFE3, DDIT3::FUS, CREB3L2::EWSR1, EWSR1::NR4A3, CREB3L2::FUS, NR4A3::SMARCA2, ETV6::NTRK3, EWSR1::WT1
- tumor example: `data/scvi_sysa_feb26.h5ad`.  
    - 112,903 (cells/barcodes) × 34,876 (genes/transcripts). 
    - Entity: SYSA
    - Fusion: SS18::SSX1, SS18::SSX2
- downsampled 20%: `data/scvi_pools_2026_downsampled_20per.h5ad`.  
    - 96,282 (cells/barcodes) × 35,929 (genes/transcripts). 


# Understanding scanpy and anndata

## Preprocessing

The files storing anndata objects (.h5ad) are pretty big and therefore hard to process. To actually work with them, it is sometimes necessary to work with job submissions. To run jobs on python scripts use: 
```bash
submit-job -m 250g -t 02:00:00 -c 10 -n scanpy_env -e experiment01 -j downsample_20per -s 04_downsample.py
```

#### Histogram of expression values

- downsampled 20%: `data/scvi_pools_2026_downsampled_20per.h5ad`.  

![Expression Histogram](images/adata_20per_hist_expression.png)
![Expression Histogram](images/adata_20per_hist_counts.png)

#### Quality Control via [scanpy.pp.calculate_qc_metrics](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html) and [scanpy.pp.filter_cells](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.filter_cells.html)
- GOAL: make it clean
- problems: low-quality cells,  contamination cell-free RNA, doublets

```python
# mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
scanpy.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], inplace=True, log1p=False
)
# example filter
adata = adata[adata.obs["pct_counts_mt"] < 10, :]
scanpy.pp.filter_cells(adata, min_genes=100) # cells with min deteced genes
scanpy.pp.filter_genes(adata, min_cells=3) # genes with min cell count 
```

- Filtering low quality cells
    - criteria for dying cell: 
        - low number of detected genes
        - low count depth
        - high fraction of mitochondrial counts (because of broken membrane)
            - mitochondrial counts: all captured RNA molecules in a single cell that belong to mitochondrial genes (start with MT)
        - very low fraction of ribosomal counts (no protein synthesis)
    - red blood cell contamination
        - high fraction of hemoglobin counts
- QC metrics
    - `obs` columns added: n_genes_by_counts, total_counts, total_counts_mt, pct_counts_mt
        - `n_genes_by_counts`: number of genes with (positive) counts in a cell 
        - `total_counts`: total number of counts for a cell (library size)
        - `total_counts_mt`:  total number of counts for a cell which are mitochondrial
        - `pct_counts_mt`:  proportion of total counts for a cell which are mitochondrial
    - `var` columns added: mt, total_counts, n_cells_by_count, mean_counts, ct_dropout_by_counts
        - `total_counts`: sum of counts for a gene
        - `mean_counts`: mean expression over all cells (for that gene)
        - `n_cells_by_count`: number of cells this expression is measured in (>0)
        - `ct_dropout_by_counts`: percentage of cells this feature does not appear in

#### Normalization (TP10K) via [scanpy.pp.normalize_total](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_total.html)
- GOAL: adjust to sequencing depth of a cell (count depth scaling)
```python
scanpy.pp.normalize_total(adata,target_sum=1e4)
```
- normalize each cell by set number like 10,000 reads (target_sum=1e4)
- $\text{Normalized Value}=(\frac{\text{Raw Count for Gene in Cell}}{\text{Total Counts in Cell}})\times {target\_sum}$
```python
scanpy.pp.normalize_total(adata,target_sum=None)
```
- the default target sum is set to be the median raw count depth in the dataset.
- $\text{Normalized Value}=(\frac{\text{Raw Count for Gene in Cell}}{\text{Total Counts in Cell}})\times \text{Median Raw Count Depth}$

#### Shifted Logarithm via [scanpy.pp.log1p](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.log1p.html)
- GOAL 1: reduce influence of (extreme) outliers (stabilizing variance)
- GOAL 2: highlight proportional shifts in gene expression scaling rather than raw differences 
```python
scanpy.pp.log1p(adata)
```

- log transformation: compresses large values, stretches small values
    - stabilizing variance: mean gene expression should not influence the variance
    - original problem: the mean increases -> the variance also increases (genes with generally high gene expression have higher variance)
- a small shift (pseudo-count = 1) is added so that zero values can be handled
- $X = log(X + 1)$

## feature selection
#### highly variable genes (HVGs) with [scanpy.pp.highly_variable_genes](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html)

 - HVG (informative genes): the gene expression of highly variable genese varies significantly throughout the cells 

```python
scanpy.pp.highly_variable_genes(adata,flavor='seurat',n_bins=20,n_top_genes=None)
```
- input: log1p normalized expression values 
- `seurat`
    - `dispersion`: measures how far apart or scattered data values are from a central value like the mean
        -  $d_{g} = log(\frac{\sigma_g^2}{\mu_g})$
            - $\mu_g = \frac{1}{N}\sum_{i=1}^{N}{y_{ig}}\quad, \quad\sigma_g^2 = \frac{1}{N-1}\sum_{i=1}^{N}{(y_{ig} - \mu_g)^2}$
        - computed in non-log space: the expr values $y$ are exponentiated beforehand (reverse of log1p transformation)
    - Expression binning
        - sorted by $log(\mu _{g})$
        - split into 20 discrete groups
        - only compare genes that have similar average expression levels
    - `Z-Score` for each gene: compare dispersions over genes in same bin
        - How many standard deviations away is this gene's variance compared to its peers with the exact same expression level?
        - $Z_g = \frac{D_g - \overline{D}_b}{SD_b}$
- highly variable: dispersion of gene expression among cells is greater than in other cells (high positive Z-score)
- example filtering and selection
    - `n_top_genes` : top N genes by Z-score
    - `cutoffs` (used if n_top_genes=None): 
        - `min_disp` = 0.5 (all genes with Z-score >= 0.5)
        - `max_disp` = inf
        - `min_mean` = 0.0125 (filter out unexpressed genes)
        - `max_mean` = 3 (filter out extreme housekeeping genes that dominate expr)

## Dimension reduction
#### PCA: dimension reduction with [scanpy.tl.pca](https://scanpy.readthedocs.io/en/1.9.x/generated/scanpy.tl.pca.html)
- GOAL: reduce dimensions, denoise data, find main axes of variation (which dimension retain most variation?)
```python
scanpy.tl.pca(adata)
scanpy.pp.pca(adata)
```

#### nearest neighbor graph with [scanpy.pp.neighbors](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.neighbors.html)

- GOAL: calculate distance between each cell
```python
scanpy.pp.neighbors(adata, n_neighbors=15, use_rep='X')
```
- distance for every gene $G$, between nodes/cells $p$ and $q$
    $$d(p,q) = \sqrt{\sum_{i=1}^G{(p_i-q_i)^2}}$$
```python
scanpy.pp.neighbors(adata, n_neighbors=15, use_rep='X_PCA')
```
- calculation for every gene is too large, also "Curse of Dimensionality" -> use PCA
- distance for every dimension $n$, between nodes/cells $p$ and $q$
    $$d(p,q) = \sqrt{\sum_{i=1}^n{(p_i-q_i)^2}}$$
- obsp
    - `distances`: euclidean distance
        - values range: 0 - inf
    - `connectivities`: weights
        - distances passed through exponential decay function 
        - values range: 0 - 1
        - decay function: $$W_{ij} = \exp{(-\frac{d(i,j)-p_i}{\sigma_i})}$$


#### UMAP: dimension reduction with [scanpy.tl.umap](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.umap.html)

- GOAL: take nearest neighbor graph and reduce dimensions without loosing much information about nearest neighbors

```python
scanpy.tl.umap(adata)
```

- tumor example: `data/scvi_sysa_feb26.h5ad`.  
![Expression Histogram](images/adata_sysa_umap_celltype.png)
- downsampled 20%: `data/scvi_pools_2026_downsampled_20per.h5ad`.  
![Expression Histogram](images/adata_20per_umap_celltype.png)
![Expression Histogram](images/adata_20per_umap_Entity.png)

## Clustering

#### Leiden: Clustering with [scanpy.tl.leiden](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html)

- GOAL: cluster cells based on nearest neighbor graph
- NOTE: improved version of Louvain algorithm 

```python
scanpy.tl.leiden(adata)
```
