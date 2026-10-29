# ============================================================================ #
# IMPORTS                                                                      #
# ============================================================================ #

import os
import warnings

warnings.filterwarnings('ignore')

import json
import pickle
import sys
import time

import numpy as np
import scanpy as sc
import torch
import wandb

sys.path.insert(0, "../")
import scgpt as scg
from scgpt.model import TransformerModel
from scgpt.preprocess import Preprocessor
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.utils import set_seed

import myutils

# ============================================================================ #
# PATHS                                                                        #
# ============================================================================ #

workspace,experiment=myutils.workspace()
experiment = workspace / "experiment04-annotation"

result_dir = experiment/"results/pred-on-fine-tuned"

data_file =  experiment/ "data" / "heroes_pool8_scRNA_hvg.h5ad"  # ! Dataset
# data_file =  experiment/ "data" / "heroes_pool8_scRNA_positive.h5ad"  # ! Dataset
celltype_col = "celltype"

model_dir = experiment / "results/pred-with-fine-tuning" / 'pred_heroes-2026-10-25_172348' 
model_name = 'scGPT_pancancer'

model_file = model_dir/ 'model.pt'
model_config_file = model_dir / "args.json"
vocab_file = model_dir / "vocab.json"

# ============================================================================ #
# CONFIG                                                                       #
# ============================================================================ #

with open(model_config_file, "r") as f:
    model_configs = json.load(f)

os.environ["WANDB_MODE"] = "disabled" # ! Selina
run = wandb.init(
    config=model_configs,
    project="scGPT",
    reinit=True,
)
config = wandb.config
set_seed(config.seed)

# SCRIPT CONFIG -------------------------------------------------------------- #
config.eval_batch_size = config.batch_size
config.ECS = config.ecs_threshold > 0
config.explicit_zero_prob = config.MLM and config.include_zero_gene 
config.DAB_separate_optim = config.DAB > 1
config.hvg_flavor = "seurat_v3" if config.data_is_raw else "cell_ranger"
config.input_layer_key = {  # the values of this map coorespond to the keys in preprocessing
    "normed_raw": "X_normed",
    "log1p": "X_normed",
    "binned": "X_binned",
}[config.input_style]
if config.input_emb_style == "category":
    config.mask_value = config.n_bins + 1
    config.pad_value = config.n_bins  # for padding gene expr values
    config.n_input_bins = config.n_bins + 2
else:
    config.mask_value = -1
    config.pad_value = -2
    config.n_input_bins = config.n_bins

save_dir = result_dir/f"pred_{model_name}-{time.strftime('%Y-%m-%d_%H%M%S')}/"
save_dir.mkdir(parents=True, exist_ok=True)
print(f"CONFIG: save to {save_dir}\n")
# SCRIPT CONFIG -------------------------------------------------------------- #

# LOGGING -------------------------------------------------------------------- #
logger = scg.logger
scg.utils.add_file_handler(logger, save_dir / "run.log")
# LOGGING -------------------------------------------------------------------- #

# ============================================================================ #
# GENE VOCABULARY                                                              #
# ============================================================================ #

vocab = GeneVocab.from_file(vocab_file)
for s in config.special_tokens:
    if s not in vocab:
        vocab.append_token(s)
vocab.set_default_index(vocab["<pad>"])

# ============================================================================ #
# DATA PROCESSING                                                              #
# ============================================================================ #

# ! adata should have obs columns: celltype (celltype_col), str_batch
# ! adata var index should be gene names
# READ ----------------------------------------------------------------------- #
adata = sc.read(data_file)
# edit categorie types and index
adata.obs[celltype_col] = adata.obs[celltype_col].astype("category")
# batch ids to int
adata.obs["batch_id"]  = adata.obs["str_batch"] = "1"          
batch_id_labels = adata.obs["str_batch"].astype("category").cat.codes.values
adata.obs["batch_id"] = batch_id_labels
# cell types to int 
celltype_id_labels = adata.obs[celltype_col].astype("category").cat.codes.values
celltypes = adata.obs[celltype_col].unique()
num_types = len(np.unique(celltype_id_labels))


adata.obs["celltype_id"] = celltype_id_labels
# gene names 
adata.var["gene_name"] = adata.var.index.tolist()
adata.var.set_index(adata.var["gene_name"], inplace=True)
adata.var["id_in_vocab"] = [
    1 if gene in vocab else -1 for gene in adata.var["gene_name"]
]
gene_ids_in_vocab = np.array(adata.var["id_in_vocab"])
adata = adata[:, adata.var["id_in_vocab"] >= 0]

# READ ----------------------------------------------------------------------- #

# LOGGING -------------------------------------------------------------------- #
logger.info(
    f"match {np.sum(gene_ids_in_vocab >= 0)}/{len(gene_ids_in_vocab)} genes "
    f"in vocabulary of size {len(vocab)}."
)
# LOGGING -------------------------------------------------------------------- #

# PROCESSING ----------------------------------------------------------------- #
preprocessor = Preprocessor(
    use_key="X",  # the key in adata.layers to use as raw data
    filter_gene_by_counts=config.filter_gene_by_counts,  # step 1
    filter_cell_by_counts=config.filter_cell_by_counts,  # step 2
    normalize_total=1e4,  # 3. whether to normalize the raw data and to what sum
    result_normed_key="X_normed",  # the key in adata.layers to store the normalized data
    log1p=config.data_is_raw,  # 4. whether to log1p the normalized data
    result_log1p_key="X_log1p",
    subset_hvg=False,  # 5. whether to subset the raw data to highly variable genes
    hvg_flavor=config.hvg_flavor,
    binning=config.n_bins,  # 6. whether to bin the raw data and to what number of bins
    result_binned_key="X_binned",  # the key in adata.layers to store the binned data
)
preprocessor(adata, batch_key=None)
logger.info(f'\n{adata}\n')
# PROCESSING ----------------------------------------------------------------- #

# ============================================================================ #
# MODEL                                                                        #
# ============================================================================ #

# PREP ----------------------------------------------------------------------- #
genes = adata.var["gene_name"].tolist() 
gene_ids = np.array(vocab(genes), dtype=int) 
batch_ids = adata.obs["batch_id"].tolist() 
num_batch_types = len(set(batch_ids))
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# PREP ----------------------------------------------------------------------- #

# LOAD MODEL ----------------------------------------------------------------- #
model = TransformerModel(
    len(vocab),
    config.embsize,
    config.nheads,
    config.d_hid,
    config.nlayers,
    vocab=vocab,
    pad_value=config.pad_value,
    pad_token=config.pad_token,
    n_input_bins=config.n_input_bins,
    use_fast_transformer=config.fast_transformer,
    fast_transformer_backend=config.fast_transformer_backend,
    num_batch_labels=num_batch_types,
    use_batch_labels=config.INPUT_BATCH_LABELS,
    nlayers_cls=3,
    n_cls=num_types if config.CLS else 1,
    dropout=config.dropout,
    do_mvc=config.MVC,
    do_dab=config.DAB,
    domain_spec_batchnorm=config.DSBN,
    input_emb_style=config.input_emb_style,
    cell_emb_style=config.cell_emb_style,
    mvc_decoder_style=config.mvc_decoder_style,
    ecs_threshold=config.ecs_threshold,
    explicit_zero_prob=config.explicit_zero_prob,
    pre_norm=config.pre_norm,
)

try:
    model.load_state_dict(torch.load(model_file,map_location=torch.device('cpu')))
    logger.info(f"Loading all model params from {model_file}")
except:
    # only load params that are in the model and match the size
    model_dict = model.state_dict()
    pretrained_dict = torch.load(model_file,map_location=torch.device('cpu'))
    pretrained_dict = {
        k: v
        for k, v in pretrained_dict.items()
        if k in model_dict and v.shape == model_dict[k].shape
    }
    # for k, v in pretrained_dict.items():
        # logger.info(f"Loading params {k} with shape {v.shape}")
    model_dict.update(pretrained_dict)
    model.load_state_dict(model_dict)

model.to(device)
wandb.watch(model)
# LOAD MODEL ----------------------------------------------------------------- #


id_map_file = model_dir / "results.pkl"
with open(id_map_file, 'rb') as f:
    model_results = pickle.load(f)
id2type = model_results['id_maps']

# ============================================================================ #
# PREDICTION                                                                   #
# ============================================================================ #

predictions, labels, results = myutils.annotation.test( 
    model=model,
    adata=adata, 
    vocab=vocab, 
    device=device,
    gene_ids=gene_ids,
    config=config,
    logger=logger,
    )

save_dict = {
    "predictions": predictions,
    "results": results,
    "id_maps": id2type
}

# SAVE PKL ------------------------------------------------------------------- #
with open(save_dir / "results.pkl", "wb") as f:
    pickle.dump(save_dict, f)
# SAVE PKL ------------------------------------------------------------------- #

