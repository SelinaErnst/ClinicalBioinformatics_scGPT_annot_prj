# Practice with python packages

## Scripts
|Script | Comments|
|-----|---------|
|00_test.ipynb| test out packages and conda env|
|01_practice_anndata.ipynb| test out example adata file (SySa) |
|02_practice_scanpy.ipynb| scanpy legacy workflow (processing and clustering 3k PBMCs) |
|03_downsample.py| downsample the big adata file for easier use |
|04_plot_umaps.py| umaps for scvi_pools datasets |

## adata files

The files are all originally located here: `/omics/odcf/analysis/OE0290_projects/heroes-aya_pools/AG_Thongjuea/Code/10x_multiomics/Katherine//integration/scvi/pansarcoma`

- big adata: `data/scvi_pools/scvi_pools_2026_celltype_annotated.h5ad`.  
    - 481,411 (cells/barcodes) × 35,929 (genes/transcripts). 
    - Entity: RMA, EWING, URCS, SYSA, MYE, SFT, SEBA, ASPS, MLS, LGFMS, PLRMS, MYCHS, DDLS, LMS, PPMS, STET, cFS, DSRCT, TGCT, SARCNOS
    - Fusion: PAX7::FOXO1, PAX3::FOXO1, NCOA1::PAX3, FUS::ERG, YWHAE::NUTM2A, EWSR1::FLI1, EWSR1::ETV1, SS18::SSX1, EWSR1::POU5F1, NAB2::STAT6, SS18::SSX2, ASPSCR1::TFE3, DDIT3::FUS, CREB3L2::EWSR1, EWSR1::NR4A3, CREB3L2::FUS, NR4A3::SMARCA2, ETV6::NTRK3, EWSR1::WT1
- tumor example: `data/scvi_pools/scvi_sysa_feb26.h5ad`.  
    - 112,903 (cells/barcodes) × 34,876 (genes/transcripts). 
    - Entity: SYSA
    - Fusion: SS18::SSX1, SS18::SSX2
- downsampled 20%: `data/scvi_pools/scvi_pools_2026_downsampled_20per.h5ad`.  
    - 96,282 (cells/barcodes) × 35,929 (genes/transcripts). 

#### Umpas

- tumor example: `scvi_sysa_feb26.h5ad`.  
![](images/adata_sysa_umap_celltype.png)
- scvi pools: `scvi_pools_2026_celltype_annotated.h5ad`.  
![](images/adata_umap_celltype.png)
- downsampled 20%: `scvi_pools_2026_downsampled_20per.h5ad`.  
![](images/adata_20per_umap_celltype.png)

## Downsampling

The files storing anndata objects (.h5ad) are pretty big and therefore hard to process. To actually work with them, it is helpful to have a smaller (downsampled) version which is created via job submissions. 

> To run jobs on python scripts use [provided](../job-submission/README.md) job submission scripts. 
```bash
submit-job -m 250g -t 02:00:00 -c 10 -n scanpy_env -e experiment01 -j downsample_20per -s 04_downsample.py
```

#### Distribution of expression values / counts

- downsampled 20%: `scvi_pools_2026_downsampled_20per.h5ad`.  
![Expression Histogram](images/adata_20per_hist_expression.png)
![Counts Histogram](images/adata_20per_hist_counts.png)


