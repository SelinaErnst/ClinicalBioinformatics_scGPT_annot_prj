# ============================================================================ #
# IMPORTS                                                                      #
# ============================================================================ #

import copy
import os
import warnings

warnings.filterwarnings('ignore')

import json
import pickle
import shutil
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import torch
import wandb
from scipy.sparse import issparse
from sklearn.metrics import (
        confusion_matrix,
)
from sklearn.model_selection import train_test_split
from torch import nn

sys.path.insert(0, "../")
import scgpt as scg
from scgpt.loss import (
        masked_mse_loss,
)
from scgpt.model import AdversarialDiscriminator, TransformerModel
from scgpt.preprocess import Preprocessor
from scgpt.tokenizer import tokenize_and_pad_batch
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.utils import set_seed

import myutils

# ============================================================================ #
# PATHS                                                                        #
# ============================================================================ #

dataset = "heroes"
model_name = "scGPT_pancancer"

avail_datasets = [
    "pbmc3k",
    "ms",
    "heroes",
]
avail_models = [
    "scGPT_human",
    "scGPT_pancancer",
    "scGPT_bc",
]
assert model_name in avail_models
assert dataset in avail_datasets

workspace,experiment=myutils.workspace()
experiment = workspace / "experiment04-annotation"

data_dir = experiment/ "data"
result_dir = experiment/"results/pred-with-fine-tuning"

if (dataset == "pbmc3k"):
    data_file = data_dir / "pbmc3k_scanpy_annotated_batches.h5ad"  # ! Dataset
    celltype_col = "celltype"
if (dataset == "ms"):
    data_file = workspace / "data/ms/ms_batches.h5ad"  # ! Dataset
    celltype_col = "celltype"
if (dataset == "heroes"):
    data_file = data_dir / "hereos_pool7_scRNA_batches.h5ad" # ! Dataset
    # data_file = data_dir / "hereos_pool7_scRNA_batches_hvg.h5ad"  # ! Dataset
    # data_file = data_dir / "hereos_pool7_scRNA_batches_small.h5ad" 
    # data_file = data_dir / "hereos_pool7_scRNA_batches_small_tumor.h5ad"
    celltype_col = "consensus_celltype"

if (model_name == "scGPT_human"):
    model_dir = workspace / "data/scGPT_human"
if (model_name == "scGPT_pancancer"):
    model_dir = workspace / "data/scGPT_pancancer"
if (model_name == "scGPT_bc"):
    model_dir = workspace / "data/scGPT_bc"

model_file = model_dir/ 'best_model.pt'
model_config_file = model_dir / "args.json"
vocab_file = model_dir / "vocab.json"

# ============================================================================ #
# CONFIG                                                                       #
# ============================================================================ #

# WORKFLOW CONFIG ------------------------------------------------------------ #
hyperparameter_defaults = { 

    'seed': 1,
    'dataset_name': dataset, # ! Dataset
    'load_model': model_dir, 
    'do_train': True, 
    
    'layer_size': 128,
    'nlayers': 4,  # number of nn.TransformerEncoderLayer in nn.TransformerEncoder
    'nhead': 4,  # number of heads in nn.MultiheadAttention
    'schedule_ratio': 0.9,  # ratio of epochs for learning rate schedule
    'pre_norm': False,
    'amp': True,  # Automatic Mixed Precision
    'freeze': False, #freeze
    'DSBN': False,  # Domain-spec batchnorm

    # settings for input and preprocessing
    'n_bins': 51,
    'pad_token': "<pad>",
    'special_tokens': ["<pad>", "<cls>", "<eoc>"],
    'mask_ratio': 0.0,
    'data_is_raw': False,
    'filter_gene_by_counts': False,
    'filter_cell_by_counts': False,
    'subset_hvg': 3000, # ! HVG
    'include_zero_gene': False,
    'max_seq_len': 3001,

    # input/output representation
    'input_style': 'binned',

    # settings for training
    'MLM': False,
    'CLS': True,
    'DAB': False, 
    "ADV" : False, 
    "CCE" : False, 
    'MVC': False, # Masked value prediction for cell embedding
    'INPUT_BATCH_LABELS': False,
    'input_emb_style': "continuous",
    'cell_emb_style': "cls",
    'adv_E_delay_epochs': 0,  # delay adversarial training on encoder for a few epochs 
    'adv_D_delay_epochs': 0, 
    'mvc_decoder_style': "inner product",
    'ecs_threshold': 0.0, # Elastic cell similarity objective, 0.0 to 1.0, 0.0 to disable
    'dab_weight': 0.0,
    'do_sample_in_train': False,
    'per_seq_batch_sample': False,

    # settings for optimizer
    'lr': 1e-4,
    'lr_ADV': 1e-3,
    'batch_size': 32, ##! BATCH SIZE
    'epochs': 10, ##! EPOCH
    'schedule_interval': 1,

    # settings for the model
    'fast_transformer': True,
    'fast_transformer_backend': "flash",
    'dropout': 0.2,  # dropout probability

    # logging
    'log_interval': 100,
    'save_eval_interval': 5,
    # 'do_eval_scib_metrics': True,
}
os.environ["WANDB_MODE"] = "disabled" # ! Selina
run = wandb.init(
    config=hyperparameter_defaults,
    project="scGPT",
    reinit=True,
)
config = wandb.config
set_seed(config.seed)
# WORKFLOW CONFIG ------------------------------------------------------------ #

# VALIDATION ----------------------------------------------------------------- #
assert config.input_style in ["normed_raw", "log1p", "binned"]
assert config.input_emb_style in ["category", "continuous", "scaling"]
if config.input_style == "binned" \
    and config.input_emb_style == "scaling":
        raise ValueError("input_emb_style `scaling` is not supported for binned input.")
elif (config.input_style == "log1p" or config.input_style == "normed_raw") \
    and (config.input_emb_style == "category"):
        raise ValueError("input_emb_style `category` is not supported for log1p or normed_raw input.")
if config.ADV and config.DAB:
    raise ValueError("ADV and DAB cannot be both True.")
# VALIDATION ----------------------------------------------------------------- #

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
save_dir = result_dir/f"pred_{config.dataset_name}-{time.strftime('%Y-%m-%d_%H%M%S')}/"
save_dir.mkdir(parents=True, exist_ok=True)
print(f"CONFIG: save to {save_dir}\n")
# SCRIPT CONFIG -------------------------------------------------------------- #

# LOGGING -------------------------------------------------------------------- #
logger = scg.logger
scg.utils.add_file_handler(logger, save_dir / "run.log")
logger.info(f"Model: {model_name}")
logger.info(f"Dataset: {dataset}")
logger.info(f"File: {data_file}")
# LOGGING -------------------------------------------------------------------- #

# MODEL CONFIG --------------------------------------------------------------- #
with open(model_config_file, "r") as f:
    model_configs = json.load(f)
# MODEL CONFIG --------------------------------------------------------------- #

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
    subset_hvg=config.subset_hvg,  # 5. whether to subset the raw data to highly variable genes
    hvg_flavor=config.hvg_flavor,
    binning=config.n_bins,  # 6. whether to bin the raw data and to what number of bins
    result_binned_key="X_binned",  # the key in adata.layers to store the binned data
)
preprocessor(adata, batch_key=None)
preprocessor(adata_test, batch_key=None)
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
batch_ids = np.array(batch_ids)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# PREP ----------------------------------------------------------------------- #

# LOGGING -------------------------------------------------------------------- #
logger.info(f"CUDA: {torch.cuda.is_available()}")
logger.info(f"TORCH DEVICE: {device}")
logger.info(f"PyTorch CPU threads: {torch.get_num_threads()}")
logger.info(f"PyTorch interop threads: {torch.get_num_interop_threads()}")
myutils.experiment_config(
    path=save_dir/'args.json',

    seed=config.seed,
    dataset_name=config.dataset_name,
    do_train=config.do_train,
    load_model=config.load_model,
    device=f"{device}",

    # settings for input and preprocessing
    data_is_raw=config.data_is_raw,
    filter_cell_by_counts=config.filter_cell_by_counts,
    filter_gene_by_counts=config.filter_gene_by_counts,
    subset_hvg=config.subset_hvg, 
    hvg_flavor=config.hvg_flavor,
    pad_token=config.pad_token,
    pad_value=config.pad_value,
    special_tokens=config.special_tokens,
    mask_value=config.mask_value,
    mask_ratio=config.mask_ratio,

    include_zero_gene=config.include_zero_gene,
    max_seq_len=config.max_seq_len,
    n_bins=config.n_bins,
    per_seq_batch_sample=config.per_seq_batch_sample,

    # input/output representation
    input_style=config.input_style,
    # output_style=config.output_style,

    # settings for training
    MLM=config.MLM, # whether to use masked language modeling, currently it is always on.
    CLS=config.CLS, # celltype classification objective
    ADV=config.ADV, # Adversarial training for batch correction
    CCE=config.CCE, # Contrastive cell embedding objective
    MVC=config.MVC, # Masked value prediction for cell embedding
    ECS=config.ECS,  # Elastic cell similarity objective
    DAB=config.DAB, # Domain adaptation by reverse backpropagation, set to 2 for separate optimizer
    DSBN=config.DSBN,
    INPUT_BATCH_LABELS=config.INPUT_BATCH_LABELS,
    input_emb_style=config.input_emb_style,
    cell_emb_style=config.cell_emb_style,
    adv_E_delay_epochs=config.adv_E_delay_epochs,
    adv_D_delay_epochs=config.adv_D_delay_epochs,
    mvc_decoder_style=config.mvc_decoder_style,
    ecs_threshold=config.ecs_threshold,
    dab_weight=config.dab_weight,

    # settings for optimizer
    lr=config.lr,
    lr_ADV=config.lr_ADV,
    batch_size=config.batch_size,
    eval_batch_size=config.eval_batch_size,
    epochs=config.epochs,
    schedule_interval=config.schedule_interval,
    schedule_ratio=config.schedule_ratio,

    # settings for the model
    fast_transformer=config.fast_transformer,
    fast_transformer_backend=config.fast_transformer_backend,
    layer_size=config.layer_size,
    embsize=model_configs["embsize"],
    d_hid=model_configs["d_hid"],
    nlayers=model_configs["nlayers"],
    nheads=model_configs["nheads"],
    n_layers_cls=model_configs["n_layers_cls"],
    dropout=config.dropout,
    pre_norm=config.pre_norm,
    ntokens=len(vocab),

    freeze=config.freeze,

    # logging
    log_interval=config.log_interval,

    amp=config.amp,
)
# LOGGING -------------------------------------------------------------------- #

# LOAD MODEL ----------------------------------------------------------------- #
model = TransformerModel(
    len(vocab),
    model_configs["embsize"],
    model_configs["nheads"],
    model_configs["d_hid"],
    model_configs["nlayers"],
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
    for k, v in pretrained_dict.items():
        logger.info(f"Loading params {k} with shape {v.shape}")
    model_dict.update(pretrained_dict)
    model.load_state_dict(model_dict)

pre_freeze_param_count = sum({p.data_ptr(): p.numel() for p in model.parameters() if p.requires_grad}.values())
# Freeze all pre-decoder weights
for name, para in model.named_parameters():
    if config.freeze and "encoder" in name and "transformer_encoder" not in name:
        para.requires_grad = False
post_freeze_param_count = sum({p.data_ptr(): p.numel() for p in model.parameters() if p.requires_grad}.values())
model.to(device)
wandb.watch(model)
# LOAD MODEL ----------------------------------------------------------------- #

# LOGGING -------------------------------------------------------------------- #
logger.info(f"Total Pre freeze Params {(pre_freeze_param_count )}")
logger.info(f"Total Post freeze Params {(post_freeze_param_count )}")
wandb.log(
        {
            "info/pre_freeze_param_count": pre_freeze_param_count,
            "info/post_freeze_param_count": post_freeze_param_count,
        },
)
# LOGGING -------------------------------------------------------------------- #

# ============================================================================ #
# FINE TUNING PREP                                                             #
# ============================================================================ #

def define_wandb_metrcis():
    wandb.define_metric("valid/mse", summary="min", step_metric="epoch")
    wandb.define_metric("valid/mre", summary="min", step_metric="epoch")
    wandb.define_metric("valid/dab", summary="min", step_metric="epoch")
    wandb.define_metric("valid/sum_mse_dab", summary="min", step_metric="epoch")
    wandb.define_metric("test/avg_bio", summary="max")

all_counts = (
    adata.layers[config.input_layer_key].toarray()
    if issparse(adata.layers[config.input_layer_key])
    else adata.layers[config.input_layer_key]
)

celltypes_labels = adata.obs["celltype_id"].tolist()  # make sure count from 0
celltypes_labels = np.array(celltypes_labels)

(
    train_data,
    valid_data,
    train_celltype_labels,
    valid_celltype_labels,
    train_batch_labels,
    valid_batch_labels,
) = train_test_split(
    all_counts, celltypes_labels, batch_ids, test_size=0.1, shuffle=True
)

tokenized_train = tokenize_and_pad_batch(
    train_data,
    gene_ids,
    max_len=config.max_seq_len,
    vocab=vocab,
    pad_token=config.pad_token,
    pad_value=config.pad_value,
    append_cls=True,  # append <cls> token at the beginning
    include_zero_gene=config.include_zero_gene,
)
tokenized_valid = tokenize_and_pad_batch(
    valid_data,
    gene_ids,
    max_len=config.max_seq_len,
    vocab=vocab,
    pad_token=config.pad_token,
    pad_value=config.pad_value,
    append_cls=True,
    include_zero_gene=config.include_zero_gene,
)

logger.info(
    f"train set number of samples: {tokenized_train['genes'].shape[0]}, "
    f"\n\t feature length: {tokenized_train['genes'].shape[1]}"
)
logger.info(
    f"valid set number of samples: {tokenized_valid['genes'].shape[0]}, "
    f"\n\t feature length: {tokenized_valid['genes'].shape[1]}"
)

if config.ADV:
    discriminator = AdversarialDiscriminator(
        d_model=model_configs["embsize"],
        n_cls=num_batch_types,
    ).to(device)
else:
    discriminator = None

criterion = masked_mse_loss
optimizer = torch.optim.Adam(
    model.parameters(), lr=config.lr, eps=1e-4 if config.amp else 1e-8
)
scheduler = torch.optim.lr_scheduler.StepLR(
    optimizer, config.schedule_interval, gamma=config.schedule_ratio
)

if config.DAB_separate_optim:
    optimizer_dab = torch.optim.Adam(model.parameters(), lr=config.lr)
    scheduler_dab = torch.optim.lr_scheduler.StepLR(
        optimizer_dab, config.schedule_interval, gamma=config.schedule_ratio
    )
if config.ADV:
    criterion_adv = nn.CrossEntropyLoss()  # consider using label smoothing
    optimizer_E = torch.optim.Adam(model.parameters(), lr=config.lr_ADV)
    scheduler_E = torch.optim.lr_scheduler.StepLR(
        optimizer_E, config.schedule_interval, gamma=config.schedule_ratio
    )
    optimizer_D = torch.optim.Adam(discriminator.parameters(), lr=config.lr_ADV)
    scheduler_D = torch.optim.lr_scheduler.StepLR(
        optimizer_D, config.schedule_interval, gamma=config.schedule_ratio
    )
else:
    criterion_adv, optimizer_E, optimizer_D = None, None, None

scaler = torch.cuda.amp.GradScaler(enabled=config.amp)

# ============================================================================ #
# FINE TUNING                                                                  #
# ============================================================================ #

best_val_loss = float("inf")
best_avg_bio = 0.0
best_model = None
define_wandb_metrcis()

for epoch in range(1, config.epochs + 1):
    epoch_start_time = time.time()
    train_data_pt, valid_data_pt = myutils.annotation.prepare_data(
        tokenized_train=tokenized_train,
        train_batch_labels=train_batch_labels,
        train_celltype_labels=train_celltype_labels,
        tokenized_valid=tokenized_valid,
        valid_batch_labels=valid_batch_labels,
        valid_celltype_labels=valid_celltype_labels,
        sort_seq_batch=config.per_seq_batch_sample,
        epoch=epoch,
        config=config,
        )
    train_loader = myutils.annotation.prepare_dataloader(
        data_pt=train_data_pt,
        batch_size=config.batch_size,
        shuffle=False,
        intra_domain_shuffle=True,
        drop_last=False,
        config=config,
    )
    valid_loader = myutils.annotation.prepare_dataloader(
        data_pt=valid_data_pt,
        batch_size=config.eval_batch_size,
        shuffle=False,
        intra_domain_shuffle=False,
        drop_last=False,
        config=config,
    )

    if config.do_train:
        myutils.annotation.train(
            model=model,
            loader=train_loader,
            vocab=vocab,
            epoch=epoch,
            scheduler=scheduler,
            optimizer=optimizer,
            optimizer_D=optimizer_D,
            optimizer_E=optimizer_E,
            discriminator=discriminator,
            device=device,
            config=config,
            logger=logger,
        )
    val_loss, val_err = myutils.annotation.evaluate(
        model=model,
        loader=valid_loader,
        vocab=vocab,
        device=device,
        epoch=epoch,
        config=config,
    )
    elapsed = time.time() - epoch_start_time
    logger.info("-" * 89)
    logger.info(
        f"| end of epoch {epoch:3d} | time: {elapsed:5.2f}s | "
        f"valid loss/mse {val_loss:5.4f} | err {val_err:5.4f}"
    )
    logger.info("-" * 89)

    if val_loss < best_val_loss:
        best_val_loss = val_loss
        best_model = copy.deepcopy(model)
        best_model_epoch = epoch
        logger.info(f"Best model with score {best_val_loss:5.4f}")

    scheduler.step()
    if config.DAB_separate_optim:
        scheduler_dab.step()
    if config.ADV:
        scheduler_D.step()
        scheduler_E.step()


# ============================================================================ #
# PREDICTION                                                                   #
# ============================================================================ #

# SAVE MODEL ----------------------------------------------------------------- #
shutil.copy(vocab_file, save_dir / "vocab.json")
torch.save(best_model.state_dict(), save_dir / "model.pt")
# SAVE MODEL ----------------------------------------------------------------- #


predictions, labels, results = myutils.annotation.test( 
    model=best_model,
    adata=adata_test, 
    vocab=vocab, 
    device=device,
    gene_ids=gene_ids,
    config=config,
    logger=logger,
    )

adata_test_raw.obs["predictions"] = [id2type[p] for p in predictions]
save_dict = {
    "predictions": predictions,
    "labels": labels,
    "results": results,
    "id_maps": id2type
}

# SAVE PKL ------------------------------------------------------------------- #
# adata_test_raw.write(save_dir / f"pred_{config.dataset_name}.h5ad")
with open(save_dir / "results.pkl", "wb") as f:
    pickle.dump(save_dict, f)
# SAVE PKL ------------------------------------------------------------------- #

# UMAP ----------------------------------------------------------------------- #
palette_ = plt.rcParams["axes.prop_cycle"].by_key()["color"] 
palette_ = plt.rcParams["axes.prop_cycle"].by_key()["color"] + plt.rcParams["axes.prop_cycle"].by_key()["color"] + plt.rcParams["axes.prop_cycle"].by_key()["color"]
palette_ = {c: palette_[i] for i, c in enumerate(celltypes)}

if 'X_umap' in adata_test_raw.obsm or 'umap' in adata_test_raw.obsm:

    with plt.rc_context({"figure.figsize": (6, 4), "figure.dpi": (300)}):
        sc.pl.umap(
            adata_test_raw,
            color=[celltype_col, "predictions"],
            palette=palette_,
            show=False,
            ncols=1,
        )
        plt.savefig(save_dir / "results.png", dpi=300, bbox_inches="tight")

    results["test/cell_umap"] = wandb.Image(
        str(save_dir / "results.png"),
        caption=f"predictions macro f1 {results['test/macro_f1']:.3f}",
    )
# UMAP ----------------------------------------------------------------------- #

# CONFUSION MATRIX ----------------------------------------------------------- #
celltypes = [id2type[i] for i in range(len(celltypes))]
for i in {id2type[p] for p in predictions}:
    if i not in celltypes:
        celltypes.remove(i)
cm = confusion_matrix(labels, predictions)
cm = cm.astype("float") / cm.sum(axis=1)[:, np.newaxis]
cm = pd.DataFrame(cm, index=celltypes[:cm.shape[0]], columns=celltypes[:cm.shape[1]])
plt.figure(figsize=(10, 10))
plt.tight_layout()
sns.heatmap(cm, annot=True, fmt=".1f", cmap="Blues")
plt.tight_layout()
plt.savefig(save_dir / "confusion_matrix.png", dpi=300, bbox_inches="tight")

results["test/confusion_matrix"] = wandb.Image(
    str(save_dir / "confusion_matrix.png"),
    caption="confusion matrix",
)
# CONFUSION MATRIX ----------------------------------------------------------- #

wandb.log(results)

