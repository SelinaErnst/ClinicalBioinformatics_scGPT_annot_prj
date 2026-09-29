# ============================================================================ #
# IMPORTS                                                                      #
# ============================================================================ #

import os
import warnings
warnings.filterwarnings('ignore')
import numpy as np
from anndata import AnnData
import anndata as ad
import scanpy as sc
import myutils
import json
from scipy.sparse import issparse
from typing import List, Tuple, Dict, Union, Optional
import timeff

import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader

import scgpt as scg
from scgpt.model import TransformerModel
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.utils import set_seed
from scgpt.preprocess import Preprocessor
from scgpt.tokenizer import tokenize_and_pad_batch, random_mask_value

import wandb

# ============================================================================ #
# PATHS                                                                        #
# ============================================================================ #

workspace,experiment=myutils.workspace()
experiment = workspace / "experiment03-scGPT"

run = 'ms-Sep03-14-24'

model_dir = experiment/ f'save/dev_{run}'
model_file = model_dir/ 'model.pt'
model_config_file = model_dir / "args.json"
vocab_file = model_dir / "vocab.json"
data_dir = experiment/ "scGPT/data/ms"

# ============================================================================ #
# CONFIG                                                                       #
# ============================================================================ #

with open(model_config_file, "r") as f:
    model_configs = json.load(f)
model_configs
os.environ["WANDB_MODE"] = "disabled" # ! Selina
run = wandb.init(
    config=model_configs,
    project="scGPT",
    reinit=True,
)
config = wandb.config

input_layer_key = {  # the values of this map coorespond to the keys in preprocessing
    "normed_raw": "X_normed",
    "log1p": "X_normed",
    "binned": "X_binned",
}[config.input_style]


set_seed(config.seed)

save_dir = experiment/f"predictions/fine-tuned/{run}"
result_file = save_dir/f"pred_{config.dataset_name}.h5ad"
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
adata = sc.read(data_dir / "c_data.h5ad")
adata_test = sc.read(data_dir / "filtered_ms_adata.h5ad")

# edit categorie types and index
adata.obs["celltype"] = adata.obs["Factor Value[inferred cell type - authors labels]"].astype("category")
adata_test.obs["celltype"] = adata_test.obs["Factor Value[inferred cell type - authors labels]"].astype("category")
adata.obs["batch_id"]  = adata.obs["str_batch"] = "0"
adata_test.obs["batch_id"]  = adata_test.obs["str_batch"] = "1"          
adata.var.set_index(adata.var["gene_name"], inplace=True)
adata_test.var.set_index(adata.var["gene_name"], inplace=True)
adata_test_raw = adata_test.copy()
adata_batch0_raw = adata.copy()
adata = adata.concatenate(adata_test, batch_key="str_batch")

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


genes = adata.var["gene_name"].tolist()
gene_ids = np.array(vocab(genes), dtype=int)
# celltypes_labels = adata.obs["celltype_id"].tolist()  # make sure count from 0
# celltypes_labels = np.array(celltypes_labels)
batch_ids = adata.obs["batch_id"].tolist()
# batch_ids = np.array(batch_ids)

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

preprocessor(adata_test, batch_key=None)

# ============================================================================ #
# MODEL                                                                        #
# ============================================================================ #

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

logger.info(f"CUDA: {torch.cuda.is_available()}")
logger.info(f"TORCH DEVICE: {device}")
logger.info(f"CONFIG: {config}")
logger.info(f"PyTorch CPU threads: {torch.get_num_threads()}")
logger.info(f"PyTorch interop threads: {torch.get_num_interop_threads()}")

model = TransformerModel(
    config.ntokens,
    config.embsize,
    config.nheads,
    config.d_hid,
    config.nlayers,
    vocab=vocab,
    pad_value=config.pad_value,
    pad_token=config.pad_token,
    n_input_bins=config.n_bins + 2 if (config.input_emb_style == "category")  else config.n_bins,
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
        pad_value=config.pad_value,
        append_cls=True,  # append <cls> token at the beginning
        include_zero_gene=config.include_zero_gene,
    )

    input_values_test = random_mask_value(
        tokenized_test["values"],
        mask_ratio=config.mask_ratio,
        mask_value=config.mask_value,
        pad_value=config.pad_value,
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
adata_test_raw.write(result_file)