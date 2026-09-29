# ============================================================================ #
# IMPORTS                                                                      #
# ============================================================================ #

import os
import warnings

warnings.filterwarnings('ignore')

import json
import pickle
import shutil
import sys
import time
import traceback
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scgpt as scg
import scvi
import seaborn as sns
import torch
import wandb
from anndata import AnnData
from scgpt.model import TransformerModel
from scgpt.preprocess import Preprocessor
from scgpt.tokenizer import random_mask_value, tokenize_and_pad_batch
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.utils import set_seed
from scipy.sparse import issparse
from sklearn.metrics import (
    adjusted_rand_score,
    confusion_matrix,
    normalized_mutual_info_score,
)
from sklearn.model_selection import train_test_split
from torch import nn
from torch.nn import functional as F
from torch.utils.data import DataLoader, Dataset
from torchtext._torchtext import (
    Vocab as VocabPybind,
)
from torchtext.vocab import Vocab

sys.path.insert(0, "../")
import scgpt as scg
from scgpt import SubsetsBatchSampler
from scgpt.loss import (
    criterion_neg_log_bernoulli,
    masked_mse_loss,
    masked_relative_error,
)
from scgpt.model import AdversarialDiscriminator, TransformerModel
from scgpt.preprocess import Preprocessor
from scgpt.tokenizer import random_mask_value, tokenize_and_pad_batch
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.utils import category_str2int, eval_scib_metrics, set_seed

import myutils

# ============================================================================ #
# PATHS                                                                        #
# ============================================================================ #

workspace,experiment=myutils.workspace()
experiment = workspace / "experiment04-annotation"

data_dir = experiment/ "data"
model_dir = workspace / "experiment03-scGPT/scGPT/save/scGPT_human"

model_file = model_dir/ 'best_model.pt'
model_config_file = model_dir / "args.json"
vocab_file = model_dir / "vocab.json"

# ============================================================================ #
# CONFIG                                                                       #
# ============================================================================ #

hyperparameter_defaults = { 
    'seed': 0,
    'dataset_name': "pbmc3k",
    'do_train': True,
    'n_bins': 51,
    'layer_size': 128,
    'nlayers': 4,  # number of nn.TransformerEncoderLayer in nn.TransformerEncoder
    'nhead': 4,  # number of heads in nn.MultiheadAttention
    'dropout': 0.2,  # dropout probability
    'schedule_ratio': 0.9,  # ratio of epochs for learning rate schedule
    'save_eval_interval': 5,
    'pre_norm': False,
    'amp': True,  # Automatic Mixed Precision
    'freeze': False, #freeze
    'DSBN': False,  # Domain-spec batchnorm
    'load_model': experiment / "scGPT/save/scGPT_human", 

    # !SELINA -------------------------------------------------------------------- #
    # settings for input and preprocessing
    'pad_token': "<pad>",
    'special_tokens': ["<pad>", "<cls>", "<eoc>"],
    'mask_ratio': 0.0,
    'data_is_raw': False,
    'filter_gene_by_counts': False,
    'filter_cell_by_counts': False,
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
    'adv_E_delay_epochs': 0,  # delay adversarial training on encoder for a few epochs # !
    'adv_D_delay_epochs': 0, 
    'mvc_decoder_style': "inner product",
    'ecs_threshold': 0.0, # Elastic cell similarity objective, 0.0 to 1.0, 0.0 to disable
    'dab_weight': 0.0,
    'do_sample_in_train': False,
    'per_seq_batch_sample': False,
    # settings for optimizer
    'lr': 1e-4,
    'lr_ADV': 1e-3,
    'batch_size': 32, 
    'eval_batch_size': 32, 
    'epochs': 10,
    'schedule_interval': 1,

    # settings for the model
    'fast_transformer': True,
    'fast_transformer_backend': "flash",

    # !SELINA -------------------------------------------------------------------- #
}

with open(model_config_file, "r") as f:
    model_configs = json.load(f)
os.environ["WANDB_MODE"] = "disabled" # ! Selina
run = wandb.init(
    config=hyperparameter_defaults,
    project="scGPT",
    reinit=True,
)
config = wandb.config
ECS = config.ecs_threshold > 0
explicit_zero_prob = config.MLM and config.include_zero_gene 


# ADDITIONAL ----------------------------------------------------------------- #
input_layer_key = {  # the values of this map coorespond to the keys in preprocessing
    "normed_raw": "X_normed",
    "log1p": "X_normed",
    "binned": "X_binned",
}[config.input_style]

n_bins = config.n_bins
if config.input_emb_style == "category":
    mask_value = n_bins + 1
    pad_value = n_bins  # for padding gene expr values
    n_input_bins = n_bins + 2
else:
    mask_value = -1
    pad_value = -2
    n_input_bins = n_bins

save_dir = experiment/f"results/pred-wo-fine-tuning/pred_{config.dataset_name}-{time.strftime('%b%d-%H-%M')}/"
# ADDITIONAL ----------------------------------------------------------------- #

set_seed(0)

save_dir.mkdir(parents=True, exist_ok=True)
print(f"CONFIG: save to {save_dir}\n")
logger = scg.logger
scg.utils.add_file_handler(logger, save_dir / "run.log")

# ============================================================================ #
# GENE VOCABULARY                                                              #
# ============================================================================ #

vocab = GeneVocab.from_file(vocab_file)

for s in config.special_tokens:
    if s not in vocab:
        vocab.append_token(s)


# ============================================================================ #
# DATA PROCESSING                                                              #
# ============================================================================ #

# READ ----------------------------------------------------------------------- #
adata = sc.read(data_dir / "pbmc3k_scanpy_annotated_batches.h5ad")
# edit categorie types and index
adata.obs["celltype"] = adata.obs["celltype"].astype("category")

# batch ids to int
batch_id_labels = adata.obs["str_batch"].astype("category").cat.codes.values
adata.obs["batch_id"] = batch_id_labels
# cell types to int 
celltype_id_labels = adata.obs["celltype"].astype("category").cat.codes.values
celltypes = adata.obs["celltype"].unique()
num_types = len(np.unique(celltype_id_labels))
id2type = dict(enumerate(adata.obs["celltype"].astype("category").cat.categories))
adata.obs["celltype_id"] = celltype_id_labels
# gene names as column
adata.var["gene_name"] = adata.var.index.tolist()

adata.var["id_in_vocab"] = [
    1 if gene in vocab else -1 for gene in adata.var["gene_name"]
]
gene_ids_in_vocab = np.array(adata.var["id_in_vocab"])
adata = adata[:, adata.var["id_in_vocab"] >= 0]


adata_test = adata[adata.obs["str_batch"] == "1"]
adata = adata[adata.obs["str_batch"] == "0"]
adata_test_raw = adata_test.copy()
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
    hvg_flavor="seurat_v3" if config.data_is_raw else "cell_ranger",
    binning=config.n_bins,  # 6. whether to bin the raw data and to what number of bins
    result_binned_key="X_binned",  # the key in adata.layers to store the binned data
)

preprocessor(adata, batch_key=None)

genes = adata.var["gene_name"].tolist() 
gene_ids = np.array(vocab(genes), dtype=int) 
batch_ids = adata.obs["batch_id"].tolist() 
num_batch_types = len(set(batch_ids))

# ============================================================================ #
# MODEL                                                                        #
# ============================================================================ #

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# utils.experiment_config()
logger.info(f"CUDA: {torch.cuda.is_available()}")
logger.info(f"TORCH DEVICE: {device}")
logger.info(f"CONFIG: {config}")
logger.info(f"PyTorch CPU threads: {torch.get_num_threads()}")
logger.info(f"PyTorch interop threads: {torch.get_num_interop_threads()}")

model = TransformerModel(
    len(vocab),
    model_configs["embsize"],
    model_configs["nheads"],
    model_configs["d_hid"],
    model_configs["nlayers"],
    vocab=vocab,
    pad_value=pad_value,
    pad_token=config.pad_token,
    n_input_bins=n_input_bins,
    use_fast_transformer=config.fast_transformer,
    fast_transformer_backend=config.fast_transformer_backend,
    num_batch_labels=len(set(batch_ids)),
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
    explicit_zero_prob=config.MLM and config.include_zero_gene,
    pre_norm=config.pre_norm,
)

try:
    model.load_state_dict(torch.load(model_file,map_location=torch.device('cpu')))
    print(f"Loading all model params from {model_file}")
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
        # print(f"Loading params {k} with shape {v.shape}")
        model_dict.update(pretrained_dict)
        model.load_state_dict(model_dict)

model.to(device)


# ============================================================================ #
# FUNCTIONS AND CLASSES                                                        #
# ============================================================================ #

criterion_cls = nn.CrossEntropyLoss()
criterion_dab = nn.CrossEntropyLoss()

# dataset
class SeqDataset(Dataset):
    def __init__(self, data: Dict[str, torch.Tensor]):
        self.data = data

    def __len__(self):
        return self.data["gene_ids"].shape[0]

    def __getitem__(self, idx):
        return {k: v[idx] for k, v in self.data.items()}

# EVALUATION ----------------------------------------------------------------- #

def evaluate(model: nn.Module, loader: DataLoader, config, return_raw: bool = False) -> float:
    """
    Evaluate the model on the evaluation data.
    """
    model.eval()
    total_loss = 0.0
    total_error = 0.0
    total_dab = 0.0
    total_num = 0
    predictions = []
    with torch.no_grad():
        for batch_data in loader:
            input_gene_ids = batch_data["gene_ids"].to(device)
            input_values = batch_data["values"].to(device)
            target_values = batch_data["target_values"].to(device)
            batch_labels = batch_data["batch_labels"].to(device)
            celltype_labels = batch_data["celltype_labels"].to(device)

            src_key_padding_mask = input_gene_ids.eq(vocab[config.pad_token])
            with torch.cuda.amp.autocast(enabled=config.amp):
                output_dict = model(
                    input_gene_ids,
                    input_values,
                    src_key_padding_mask=src_key_padding_mask,
                    batch_labels=batch_labels if config.INPUT_BATCH_LABELS or config.DSBN else None,
                    CLS=config.CLS,  # evaluation does not need CLS or CCE
                    CCE=False,
                    MVC=False,
                    ECS=False,
                    do_sample=False,
                    #generative_training = False,
                )
                output_values = output_dict["cls_output"]
                loss = criterion_cls(output_values, celltype_labels)

                if config.DAB:
                    loss_dab = criterion_dab(output_dict["dab_output"], batch_labels)

            total_loss += loss.item() * len(input_gene_ids)
            accuracy = (output_values.argmax(1) == celltype_labels).sum().item()
            total_error += (1 - accuracy / len(input_gene_ids)) * len(input_gene_ids)
            total_dab += loss_dab.item() * len(input_gene_ids) if config.DAB else 0.0
            total_num += len(input_gene_ids)
            preds = output_values.argmax(1).cpu().numpy()
            predictions.append(preds)

    wandb.log(
        {
            "valid/mse": total_loss / total_num,
            "valid/err": total_error / total_num,
            "valid/dab": total_dab / total_num,
            "valid/sum_mse_dab": (total_loss + config.dab_weight * total_dab) / total_num,
        },
    )

    if return_raw:
        return np.concatenate(predictions, axis=0)

    return total_loss / total_num, total_error / total_num

# TESTING -------------------------------------------------------------------- #

def test(model: nn.Module, adata: DataLoader, config) -> float:
    all_counts = (
        adata.layers[input_layer_key].toarray()
        if issparse(adata.layers[input_layer_key])
        else adata.layers[input_layer_key]
    )

    celltypes_labels = adata.obs["celltype_id"].tolist()  # make sure count from 0
    celltypes_labels = np.array(celltypes_labels)

    batch_ids = adata.obs["batch_id"].tolist()
    batch_ids = np.array(batch_ids)

    tokenized_test = tokenize_and_pad_batch(
        all_counts,
        gene_ids,
        max_len=config.max_seq_len,
        vocab=vocab,
        pad_token=config.pad_token,
        pad_value=pad_value,
        append_cls=True,  # append <cls> token at the beginning
        include_zero_gene=config.include_zero_gene,
    )

    input_values_test = random_mask_value(
        tokenized_test["values"],
        mask_ratio=config.mask_ratio,
        mask_value=mask_value,
        pad_value=pad_value,
    )

    test_data_pt = {
        "gene_ids": tokenized_test["genes"],
        "values": input_values_test,
        "target_values": tokenized_test["values"],
        "batch_labels": torch.from_numpy(batch_ids).long(),
        "celltype_labels": torch.from_numpy(celltypes_labels).long(),
    }

    test_loader = DataLoader(
        dataset=SeqDataset(test_data_pt),
        batch_size=config.eval_batch_size,
        shuffle=False,
        drop_last=False,
        num_workers=min(len(os.sched_getaffinity(0)), config.eval_batch_size // 2),
        pin_memory=True,
    )

    model.eval()
    predictions = evaluate(
        model,
        loader=test_loader,
        config=config,
        return_raw=True,
    )

    # compute accuracy, precision, recall, f1
    from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score

    accuracy = accuracy_score(celltypes_labels, predictions)
    precision = precision_score(celltypes_labels, predictions, average="macro")
    recall = recall_score(celltypes_labels, predictions, average="macro")
    macro_f1 = f1_score(celltypes_labels, predictions, average="macro")

    logger.info(
        f"Accuracy: {accuracy:.3f}, Precision: {precision:.3f}, Recall: {recall:.3f}, "
        f"Macro F1: {macro_f1:.3f}"
    )

    results = {
        "test/accuracy": accuracy,
        "test/precision": precision,
        "test/recall": recall,
        "test/macro_f1": macro_f1,
    }

    return predictions, celltypes_labels, results


# ============================================================================ #
# PREDICTION                                                                   #
# ============================================================================ #

predictions, labels, results = test(model, adata_test, config)

adata_test_raw.obs["predictions"] = [id2type[p] for p in predictions]
adata_test_raw.write(save_dir / f"pred_{config.dataset_name}.h5ad")

save_dict = {
    "predictions": predictions,
    "labels": labels,
    "results": results,
    "id_maps": id2type
}

with open(save_dir / "results.pkl", "wb") as f:
    pickle.dump(save_dict, f)

# plot
palette_ = plt.rcParams["axes.prop_cycle"].by_key()["color"] 
palette_ = plt.rcParams["axes.prop_cycle"].by_key()["color"] + plt.rcParams["axes.prop_cycle"].by_key()["color"] + plt.rcParams["axes.prop_cycle"].by_key()["color"]
palette_ = {c: palette_[i] for i, c in enumerate(celltypes)}

with plt.rc_context({"figure.figsize": (6, 4), "figure.dpi": (300)}):
    sc.pl.umap(
        adata_test_raw,
        color=["celltype", "predictions"],
        palette=palette_,
        show=False,
        ncols=1,
    )
    plt.savefig(save_dir / "results.png", dpi=300, bbox_inches="tight")

results["test/cell_umap"] = wandb.Image(
    str(save_dir / "results.png"),
    caption=f"predictions macro f1 {results['test/macro_f1']:.3f}",
)

from sklearn.metrics import confusion_matrix

celltypes = list(celltypes)
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

wandb.log(results)
