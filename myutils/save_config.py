
import json
from pathlib import Path
from typing import Any

import numpy as np


def make_json_serializable(obj: Any):
    """Convert common Python/NumPy/PyTorch objects into JSON-safe values."""
    if isinstance(obj, Path):
        return str(obj)

    if isinstance(obj, (np.integer,)):
        return int(obj)

    if isinstance(obj, (np.floating,)):
        return float(obj)

    if isinstance(obj, (np.bool_,)):
        return bool(obj)

    if isinstance(obj, dict):
        return {
            str(k): make_json_serializable(v)
            for k, v in obj.items()
        }

    if isinstance(obj, (list, tuple)):
        return [make_json_serializable(v) for v in obj]

    if isinstance(obj, set):
        return [make_json_serializable(v) for v in obj]

    return obj


def save_json(data: dict, path: Path):
    with open(path, "w") as f:
        json.dump(
            make_json_serializable(data),
            f,
            indent=4,
            sort_keys=True,
        )


# ---------------------------------------------------------------------------- #
# ! ADV: script [bool]
# ! adv_E_delay_epochs: script [int]
# ! adv_D_delay_epochs: script [int]
# ! amp: *wandb [bool]
# ! batch_size: *wandb, args, script [int]
# ! CCE: *script [bool]
# ! cell_emb_style: *script [str]
# ! CLS: script [bool]
# data_source: args [str]
# ! DAB: *script [bool]
## DAB_separate_optim: script [bool]
# ! DSBN: *wandb [bool]
# ! data_is_raw: script [bool]
# ! dataset_name: *wandb, script [str]
# ! dab_weight: *wandb, script [float]
# ! device: *script [str]
# ! d_hid: *args, script [int]
# dist_backend: args [str]
# dist_url: args [str]
# distributed: args [bool]
# do_eval_scib_metrics: script [bool]
## do_sample_in_train: script [bool]
# ! do_train: *wandb [bool]
# ! dropout: *wandb, args, script [float]
# ! ECS: script [bool]
# ! ecs_thres: *wandb [float]
## ecs_threshold: script [float] - ecs_thres
# ! embsize: *args, script [int]
# ! epochs: *wandb, args, script [int]
# ! eval_batch_size: args, *script [int]
## explicit_zero_prob: script [bool]
# ! fast_transformer: *wandb, args, script [bool]
# ! fast_transformer_backend: *script [str]                      ["linear", "flash"]
# ! filter_cell_by_counts: *script [bool]
# ! filter_gene_by_counts: *script [bool]
# fp16: args [bool]
# ! freeze: *wandb [bool]
# gpu: args [int]
# grad_accu_steps: args [int]
# ! hvg_flavor: *script [str]                                       ["seurat_v3", "cell_ranger"]
# ! INPUT_BATCH_LABELS: *script [bool]
# ! input_emb_style: args, *script [str]                          ["category", "continuous", "scaling"]
# ! include_zero_gene: *wandb, script [bool]
# ! input_style: args, *script [str]                             ["normed_raw", "log1p", "binned"]
# ! layer_size: *wandb [int]
# !? load_model: *wandb [Path], args []
# local_rank: args [int]
# ! log_interval: args, *script [int]
# ! lr: *wandb, args, script [float] 
# ! lr_ADV: script [float]
# !? mask_ratio: *wandb [float], args [list], script [float]
# ! mask_value: args, *script [int]
# ! max_seq_len: args, *script [int]
# ! MLM: *script [bool]
## model_dir: script [Path] - load_model
# ! MVC: *wandb, args, script [bool]
# ! mvc_decoder_style: *script [str]
# ! n_bins: *wandb, args, script [int]
# n_hvg: args []
## n_input_bins: *script [int]
#! n_layers_cls: *args, script [int]
# ? nhead: wandb, script, *args [int] - nheads
# !? nheads: *args [int]
# ! nlayers: wandb, *args, script [int]
# no_cls: args [bool]
# no_cce: args [bool]
# ! num_batch_types: *script [int]
# ! ntokens: *script [int]
# ! output_style: *script [str]                                    ["normed_raw", "log1p", "binned"]
# ! pad_token: wandb, args, *script [str]
# ! pad_value: args, *script [int]
# ! per_seq_batch_sample: *script [bool]
# ! pre_norm: *wandb [bool]
# rank: args [int]
# save_dir: args, *script [str]
# save_eval_interval: *wandb, script [int]
# save_interval: args [int]
# ! schedule_ratio: *wandb [float]
# ! schedule_interval: *script [int] 
# scheduler_interval: args [int]
# scheduler_factor: args [float]
# ! seed: *wandb [int]
# ! special_tokens: *script [list]
# training_tasks: args [str]
# trunc_by_sample: args [bool]
# USE_CCE: args [bool]
# USE_CLS: args [bool]
# USE_GENERATIVE_TRAINING: args [bool]
# valid_size_or_ratio: args [float]
## vocab_file: script [Path]
# vocab_path: args [str]
# warmup_ratio_or_step: args [float]
# world_size: args [int]
# ---------------------------------------------------------------------------- #

def experiment_config(path : Path | None =None, **kwargs):
    if path!=None and path.parent.exists(): save_json(kwargs,path)
    return kwargs