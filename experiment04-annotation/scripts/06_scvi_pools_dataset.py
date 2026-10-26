import sys

sys.path = [p for p in sys.path if "positron-python" not in p]

import anndata as ad

import myutils

workspace,experiment=myutils.workspace()

# ============================================================================ #
# ALL DATA                                                                     #
# ============================================================================ #

file='scvi_pools_2026_celltype_annotated.h5ad'
data=workspace/"data/scvi_pools"/file
adata = ad.read_h5ad(data)

# ============================================================================ #
# SUBSET                                                                       #
# ============================================================================ #

file='heroes_pool7_scRNA.h5ad'
data_file = workspace / 'data/scvi_pools/heroes_pool7_scRNA.h5ad'
adata_subset = ad.read_h5ad(data_file)

adata_subset_log1p = adata[adata.obs.index.isin(adata_subset.obs.index)]
adata_subset_log1p.write_h5ad(workspace/'data/scvi_pools/heroes_pool7_scRNA_positive.h5ad')
