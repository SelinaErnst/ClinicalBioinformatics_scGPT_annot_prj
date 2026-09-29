import sys

import myutils

workspace,experiment=myutils.workspace()

sys.path = [p for p in sys.path if "positron-python" not in p]
import anndata as ad
import scanpy as sc

obs_columns = ["Entity", "celltype", "Cohort"]

# ================================================================================================ #
# DOWNSAMPLED                                                                                      #
# ================================================================================================ #

file='scvi_pools_2026_downsampled_20per.h5ad'
data=workspace/"data"/file
adata = ad.read_h5ad(data)
print()
print(file)
print(adata)

for col in obs_columns:
    fig = sc.pl.umap(adata, color=col, show=False, return_fig=True)
    fig.savefig(experiment/f'images/adata_20per_umap_{col}.png', bbox_inches="tight")

# ================================================================================================ #
# BIG DATA                                                                                         #
# ================================================================================================ #

file='scvi_pools_2026_celltype_annotated.h5ad'
data=workspace/"data"/file
adata = ad.read_h5ad(data)
print()
print(file)
print(adata)


for col in obs_columns:
    fig = sc.pl.umap(adata, color=col, show=False, return_fig=True)
    fig.savefig(experiment/f'images/adata_umap_{col}.png', bbox_inches="tight")

# ================================================================================================ #
# SYSA                                                                                             #
# ================================================================================================ #

file='scvi_sysa_feb26.h5ad'
data=workspace/"data"/file
adata = ad.read_h5ad(data)
print()
print(file)
print(adata)

for col in obs_columns:
    fig = sc.pl.umap(adata, color=col, show=False, return_fig=True)
    fig.savefig(experiment/f'images/adata_sysa_umap_{col}.png', bbox_inches="tight")