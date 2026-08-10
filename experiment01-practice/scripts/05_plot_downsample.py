import sys

import utils

workspace,experiment=utils.workspace()

sys.path = [p for p in sys.path if "positron-python" not in p]
import anndata as ad
import scanpy as sc

# ================================================================================================ #
# READ ADATA                                                                                       #
# ================================================================================================ #
file='scvi_pools_2026_downsampled_20per.h5ad'
print(file)
data=workspace/"data"/file
adata = ad.read_h5ad(data)
print(adata)

fig = sc.pl.umap(adata, color="Entity", show=False, return_fig=True)
fig.savefig(experiment/'images/adata_20per_umap_Entity.png', bbox_inches="tight")
fig = sc.pl.umap(adata, color="celltype", show=False, return_fig=True)
fig.savefig(experiment/'images/adata_20per_umap_celltype.png', bbox_inches="tight")
