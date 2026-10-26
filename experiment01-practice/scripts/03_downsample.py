import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import myutils

workspace,experiment=myutils.workspace()

# workspace = next(p for p in Path.cwd().resolve().parents if p.name.lower() == "selina")/"project"
# experiment = next(p for p in Path(__file__).resolve().parents if "experiment" in p.name.lower())

sys.path = [p for p in sys.path if "positron-python" not in p]
import anndata as ad
import scanpy as sc

# ================================================================================================ #
# READ ADATA                                                                                       #
# ================================================================================================ #

file='scvi_pools_2026_celltype_annotated.h5ad'
print(file)
data=workspace/"data/scvi_pools"/file
adata = ad.read_h5ad(data)
print(adata)

# ================================================================================================ #
# DOWNSAMPLE                                                                                       #
# ================================================================================================ #

file='scvi_pools_2026_downsampled_20per.h5ad'
print(file)
np.random.seed(42)
n_cells = int(adata.n_obs * 0.2)
selected_idx = np.random.choice(adata.n_obs, n_cells, replace=False)
adata_20per = adata[selected_idx].copy()
adata_20per.write_h5ad(workspace/"data/scvi_pools"/file)
print(adata_20per)

# ================================================================================================ #
# WRITE                                                                                            #
# ================================================================================================ #

numeric_cols = adata.obs.describe().columns
non_numeric_obs = adata.obs.loc[:, ~adata.obs.columns.isin(numeric_cols)]
adata_obs = non_numeric_obs.apply(lambda x: pd.Series({'number':len(x.dropna().unique()),'unique':list(x.dropna().unique())})).T
adata_obs = adata_obs.reset_index()
adata_obs = adata_obs.rename(columns = {'index': 'obs'})
adata_obs.to_csv(experiment/'results/scvi_pools_2026_obs_categorical.csv',index=False)
adata.obs.describe().to_csv(experiment/'results/scvi_pools_2026_obs_numerical.csv',index=True)

numeric_cols = adata.var.describe().columns
non_numeric_var = adata.var.loc[:, ~adata.var.columns.isin(numeric_cols)]
adata_var = non_numeric_var.apply(lambda x: pd.Series({'number':len(x.dropna().unique()),'unique':list(x.dropna().unique())})).T
adata_var = adata_var.reset_index()
adata_var = adata_var.rename(columns = {'index': 'var'})
adata_var.to_csv(experiment/'results/scvi_pools_2026_var_categorical.csv',index=False)
adata.var.describe().to_csv(experiment/'results/scvi_pools_2026_var_numerical.csv',index=True)
