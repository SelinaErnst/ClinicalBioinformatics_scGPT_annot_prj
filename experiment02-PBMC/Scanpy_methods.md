# Quality Control 
### QC metrics via [scanpy.pp.calculate_qc_metrics](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html) and [scanpy.pp.filter_cells](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.filter_cells.html)
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

# Normalization

### Normalization (TP10K) via [scanpy.pp.normalize_total](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.normalize_total.html)
- GOAL: adjust to sequencing depth of a cell (count depth scaling)
- TP10K: transcripts per 10k (total count of eacj cell sums up to 10k)
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

# Clustering

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
adata.layers["scaled"] = adata.X.toarray()
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"], layer="scaled")
sc.pp.scale(adata, max_value=10, layer="scaled")

# PCA
scanpy.tl.pca(adata_norm, n_comps=50, zero_center=True, mask_var="highly_variable", svd_solver='arpack')
scanpy.pp.pca(adata)
```
- reducing dimension by finding directions (new dimensions) that capture the greatest amount of variance in the data
    - principal component: direction in multi-dimensional space, vector of n-dimensions (n=number of genes), list of n-weights
        - from PC1, PC2, ... : every subsequent component captures next highest amoiunt of variance
        - each PC have to be uncorrelated / perpendicular to all previous components
    - PCA output: compressing cell expression
        - expression matrix: transformed and rotated (k x n -> k x m)
        - one cell in high-dimesnional space, value for gene gets multiplied by weight for gene, this happens for each PC creating m values (m = number of PCs)
- why use `highly_variable` ?
    - HVGs are more informative about cell types
    - genes are masked out if not highly variable: loading weight of exactly 0.0 across all PCs
- steps before PCA: 
    - [scanpy.pp.regress_out](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.regress_out.html)
        - regress out (mostly) unwanted sources of variation 
        - uses simple linear regression
        - prevents cells from clustering by read depth or mitochondrial content instead of actual cell types (don't treat dying vs healthy as dominant axes of variance)
    - [scanpy.pp.scale](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.scale.html)
        - shifts the expression of each gene so the mean is (0) and scales it so the variance is (1) (Z-score calculation)
        - clips extreme outliers at (10) standard deviations
        - don't let PCA prioritize genes with higher variance 

## Neighborhood graph

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


## Clustering

#### Leiden: Clustering with [scanpy.tl.leiden](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html)

- GOAL: cluster cells based on nearest neighbor graph
- NOTE: improved version of Louvain algorithm 

```python
scanpy.tl.leiden(adata, flavor='igraph',resolution=1,n_iterations=2,directed=False,use_weights=True)
scanpy.tl.leiden(adata, flavor='igraph',resolution=0.7,n_iterations=2,directed=False)
```

- parameter
    - `flavor`: which package's implementation to use (default: leidenalg, faster: igraph)
    - `resolution`: higher value leads to more clusters
    - `n_iterations`: How many iterations of the Leiden clustering algorithm to perform (default: -1). -1 algorithm runs until optimum, 2 is faster
    - `directed`: whether to treat the graph as directed or undirected
    - `use_weights`: If True, edge weights from the graph are used in the computation (placing more emphasis on stronger edges)


```python
sc.tl.paga(adata_norm, groups='leiden')
sc.tl.umap(adata, init_pos='paga')
```
- [scanpy.tl.paga](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.paga.html)
    - PAGA: partition-based graph abstraction
    - GOAL: spatial layout actually means something biologically (incorporate connection amongst clusters)
        - not simply differentiate between cell types but also look at cell trajectories (cells changing over time)
        - cells clusters that are somewhat related are spatially closer
    - requirement: sc.pp.neighbors, sc.tl.leiden
    - how it works
        - graph partitioning: provide cell clusters (e.g. leiden)
        - statistical connectivity scoring
        - PAGA confidence value: 0 - 1 (0 = completely isolated cell tyep, higher score suggests continuous biological transition)

# Annotation
## Marker Genes

- marker genes are used to distinguish between cell types
- when comparing clusters, there are genes that are only highly expressed in some

### Differential Expression with [scanpy.tl.rank_genes_groups](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html#scanpy-tl-rank-genes-groups)

- GOAL: find up-/down-regulated genes 
```python
scanpy.tl.rank_genes_groups(
    adata,
    groupby='leiden',
    mask_var='highly_variable',
    use_raw=False,
    groups='all',
    reference='rest',
    method='t-test'
    )

# helping functions
scanpy.get.rank_genes_groups_df(adata,group="0")
scanpy.pl.rank_genes_groups(adata_norm, n_genes=25, sharey=False)
scanpy.pl.rank_genes_groups_violin(adata, groups="0", n_genes=8)
```

- compare cluster to other cells (rest of the data or just subset of data) and find genes that are up- or down-regulated
    - compare: one cluster 0 to other cluster 1 (groups: ["0"], reference: "1")
    - compare: one cluster 0 to all other clusters (groups: "all", reference: "0")
- methods:
    - `wilcoxon`: Wilcoxon rank-sum (Mann-Whitney U)
        - Standard choice. Highly robust, non-parametric test that works on almost any single-cell dataset.
        - Slower on massive datasets 
    - `t-test`: Student's t-test
        - Quick exploratory analysis. Assumes a normal distribution.
        - Fast, but prone to false positives due to single-cell data sparsity.
    - `logreg`: Logistic Regression
        - Modern, highly discrete separation. Fits a multinomial model to classify cells.
        - Great at finding unique combinations of markers; slower to run.
