import sys

import matplotlib.pyplot as plt

sys.path = [p for p in sys.path if "positron-python" not in p]

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

import myutils

workspace,experiment=myutils.workspace()

from scgpt.tokenizer.gene_tokenizer import GeneVocab

data_dir = experiment/ "data"
data_file = data_dir / "hereos_pool7_scRNA_batches.h5ad"
celltype_col = "consensus_celltype"
model_dir = workspace / "data/scGPT_pancancer"
vocab_file = model_dir / "vocab.json"


vocab = GeneVocab.from_file(vocab_file)
for s in ["<pad>", "<cls>", "<eoc>"]:
    if s not in vocab:
        vocab.append_token(s)
vocab.set_default_index(vocab["<pad>"])


# READ ----------------------------------------------------------------------- #
adata = sc.read(data_file)
# edit categorie types and index
adata.obs[celltype_col] = adata.obs[celltype_col].astype("category")
# batch ids to int
batch_id_labels = adata.obs["str_batch"].astype("category").cat.codes.values
adata.obs["batch_id"] = batch_id_labels
# cell types to int 
celltype_id_labels = adata.obs[celltype_col].astype("category").cat.codes.values
celltypes = adata.obs[celltype_col].unique()
num_types = len(np.unique(celltype_id_labels))
id2type = dict(enumerate(adata.obs[celltype_col].astype("category").cat.categories))
adata.obs["celltype_id"] = celltype_id_labels
# gene names 
adata.var["gene_name"] = adata.var.index.tolist()
adata.var.set_index(adata.var["gene_name"], inplace=True)
adata.var["id_in_vocab"] = [
    1 if gene in vocab else -1 for gene in adata.var["gene_name"]
]
gene_ids_in_vocab = np.array(adata.var["id_in_vocab"])
adata = adata[:, adata.var["id_in_vocab"] >= 0]
adata_test = adata[adata.obs["str_batch"] == "1"]
adata = adata[adata.obs["str_batch"] == "0"]
adata_test_raw = adata_test.copy()

from scgpt.preprocess import Preprocessor

preprocessor = Preprocessor(
    use_key="X",  # the key in adata.layers to use as raw data
    filter_gene_by_counts=False,  # step 1
    filter_cell_by_counts=False,  # step 2
    normalize_total=1e4,  # 3. whether to normalize the raw data and to what sum
    result_normed_key="X_normed",  # the key in adata.layers to store the normalized data
    log1p=False,  # 4. whether to log1p the normalized data
    result_log1p_key="X_log1p",
    subset_hvg=3000,  # 5. whether to subset the raw data to highly variable genes
    hvg_flavor="cell_ranger",
    binning=51,  # 6. whether to bin the raw data and to what number of bins
    result_binned_key="X_binned",  # the key in adata.layers to store the binned data
)

adata_copy = adata.copy()
adata_test_copy = adata_test.copy()
preprocessor(adata_copy, batch_key=None)
preprocessor(adata_test_copy, batch_key=None)