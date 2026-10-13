

import time
import warnings

import torch
import wandb
from scgpt.loss import (
    criterion_neg_log_bernoulli,
    masked_mse_loss,
)
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from torch import nn
from torch.utils.data import DataLoader

criterion = masked_mse_loss
criterion_cls = nn.CrossEntropyLoss()
criterion_dab = nn.CrossEntropyLoss()
criterion_adv = nn.CrossEntropyLoss()

from logging import Logger

from scgpt.model import AdversarialDiscriminator
from wandb.sdk.wandb_config import Config


def train(
    model: nn.Module, 
    loader: DataLoader, 
    vocab: GeneVocab, 
    epoch: int,
    scheduler: torch.optim.lr_scheduler.StepLR,
    optimizer: torch.optim.Adam,
    optimizer_D: torch.optim.Adam, 
    optimizer_E: torch.optim.Adam, 
    discriminator: AdversarialDiscriminator, 
    device: torch.device,
    config: Config,
    logger: Logger,
    ) -> None:
    """
    Train the model for one epoch.
    """

    scaler = torch.cuda.amp.GradScaler(enabled=config.amp)

    model.train()
    (
        total_loss,
        total_mse,
        total_cls,
        total_cce,
        total_mvc,
        total_ecs,
        total_dab,
        total_adv_E,
        total_adv_D,
        total_zero_log_prob,
        total_mvc_zero_log_prob,
    ) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    total_error = 0.0
    start_time = time.time()

    num_batches = len(loader)
    for batch, batch_data in enumerate(loader):

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
                CLS=config.CLS,
                CCE=config.CCE,
                MVC=config.MVC,
                ECS=config.ECS,
                do_sample=config.do_sample_in_train,
                #generative_training=False
            )
            masked_positions = input_values.eq(config.mask_value)  # the postions to predict
            loss = 0.0
            metrics_to_log = {}
            if config.MLM:
                loss_mse = criterion(
                    output_dict["mlm_output"], target_values, masked_positions
                )
                loss = loss + loss_mse
                metrics_to_log = {"train/mse": loss_mse.item()}
            if config.explicit_zero_prob:
                loss_zero_log_prob = criterion_neg_log_bernoulli(
                    output_dict["mlm_zero_probs"], target_values, masked_positions
                )
                loss = loss + loss_zero_log_prob
                metrics_to_log.update({"train/nzlp": loss_zero_log_prob.item()})
            if config.CLS:
                loss_cls = criterion_cls(output_dict["cls_output"], celltype_labels)
                loss = loss + loss_cls
                metrics_to_log.update({"train/cls": loss_cls.item()})

                error_rate = 1 - (
                    (output_dict["cls_output"].argmax(1) == celltype_labels)
                    .sum()
                    .item()
                ) / celltype_labels.size(0)
            if config.CCE:
                loss_cce = 10 * output_dict["loss_cce"]
                loss = loss + loss_cce
                metrics_to_log.update({"train/cce": loss_cce.item()})
            if config.MVC:
                loss_mvc = criterion(
                    output_dict["mvc_output"], target_values, masked_positions
                )
                loss = loss + loss_mvc
                metrics_to_log.update({"train/mvc": loss_mvc.item()})
            if config.MVC and config.explicit_zero_prob:
                loss_mvc_zero_log_prob = criterion_neg_log_bernoulli(
                    output_dict["mvc_zero_probs"], target_values, masked_positions
                )
                loss = loss + loss_mvc_zero_log_prob
                metrics_to_log.update({"train/mvc_nzlp": loss_mvc_zero_log_prob.item()})
            if config.ECS:
                loss_ecs = 10 * output_dict["loss_ecs"]
                loss = loss + loss_ecs
                metrics_to_log.update({"train/ecs": loss_ecs.item()})
            if config.DAB:
                # try weighting and separate optimizer
                loss_dab = criterion_dab(output_dict["dab_output"], batch_labels)
                loss = loss + config.dab_weight * loss_dab
                metrics_to_log.update({"train/dab": loss_dab.item()})

        model.zero_grad()
        scaler.scale(loss).backward()
        scaler.unscale_(optimizer)
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings("always")
            torch.nn.utils.clip_grad_norm_(
                model.parameters(),
                1.0,
                error_if_nonfinite=not scaler.is_enabled(),
            )
            if len(w) > 0:
                logger.warning(
                    f"Found infinite gradient. This may be caused by the gradient "
                    f"scaler. The current scale is {scaler.get_scale()}. This warning "
                    "can be ignored if no longer occurs after autoscaling of the scaler."
                )
        scaler.step(optimizer)
        scaler.update()

        if config.ADV:
            # rerun the model for adversarial training
            output_dict = model(
                input_gene_ids,
                input_values,
                src_key_padding_mask=src_key_padding_mask,
                batch_labels=batch_labels if config.INPUT_BATCH_LABELS or config.DSBN else None,
                CLS=config.CLS,
                CCE=config.CCE,
                MVC=config.MVC,
                ECS=config.ECS,
                do_sample=config.do_sample_in_train,
                #generative_training=False
            )

            # TRAINING DISCRIMINATOR
            loss_adv_D = criterion_adv(
                discriminator(output_dict["cell_emb"].detach()), batch_labels
            )
            if epoch > config.adv_D_delay_epochs:
                discriminator.zero_grad()
                loss_adv_D.backward()
                optimizer_D.step()

            # TRAINING ENCODER
            loss_adv_E = -criterion_adv(
                discriminator(output_dict["cell_emb"]), batch_labels
            )
            # NOTE: the loss is negative here because we want to maximize
            # the cross_entropy_loss, in other words, disguise against the discriminator
            if epoch > config.adv_E_delay_epochs:
                model.zero_grad()
                discriminator.zero_grad()
                loss_adv_E.backward()
                optimizer_E.step()

        wandb.log(metrics_to_log)

        total_loss += loss.item()
        total_mse += loss_mse.item() if config.MLM else 0.0
        total_cls += loss_cls.item() if config.CLS else 0.0
        total_cce += loss_cce.item() if config.CCE else 0.0
        total_mvc += loss_mvc.item() if config.MVC else 0.0
        total_ecs += loss_ecs.item() if config.ECS else 0.0
        total_dab += loss_dab.item() if config.DAB else 0.0
        total_adv_E += loss_adv_E.item() if config.ADV else 0.0
        total_adv_D += loss_adv_D.item() if config.ADV else 0.0
        total_zero_log_prob += loss_zero_log_prob.item() if config.explicit_zero_prob else 0.0
        total_mvc_zero_log_prob += (
            loss_mvc_zero_log_prob.item() if config.MVC and config.explicit_zero_prob else 0.0
        )
        total_error += error_rate
        if batch % config.log_interval == 0 and batch > 0:
            lr = scheduler.get_last_lr()[0]
            ms_per_batch = (time.time() - start_time) * 1000 / config.log_interval
            cur_loss = total_loss / config.log_interval
            cur_mse = total_mse / config.log_interval
            cur_cls = total_cls / config.log_interval if config.CLS else 0.0
            cur_cce = total_cce / config.log_interval if config.CCE else 0.0
            cur_mvc = total_mvc / config.log_interval if config.MVC else 0.0
            cur_ecs = total_ecs / config.log_interval if config.ECS else 0.0
            cur_dab = total_dab / config.log_interval if config.DAB else 0.0
            cur_adv_E = total_adv_E / config.log_interval if config.ADV else 0.0
            cur_adv_D = total_adv_D / config.log_interval if config.ADV else 0.0
            cur_zero_log_prob = (
                total_zero_log_prob / config.log_interval if config.explicit_zero_prob else 0.0
            )
            cur_mvc_zero_log_prob = (
                total_mvc_zero_log_prob / config.log_interval
                if config.MVC and config.explicit_zero_prob
                else 0.0
            )
            cur_error = total_error / config.log_interval
            logger.info(
                f"| epoch {epoch:3d} | {batch:3d}/{num_batches:3d} batches | "
                f"lr {lr:05.4f} | ms/batch {ms_per_batch:5.2f} | "
                f"loss {cur_loss:5.2f} | "
                + (f"mse {cur_mse:5.2f} | mre {cur_error:5.2f} |" if config.MLM else "")
                + (f"cls {cur_cls:5.2f} | " if config.CLS else "")
                + (f"err {cur_error:5.2f} | " if config.CLS else "")
                + (f"cce {cur_cce:5.2f} |" if config.CCE else "")
                + (f"mvc {cur_mvc:5.2f} |" if config.MVC else "")
                + (f"ecs {cur_ecs:5.2f} |" if config.ECS else "")
                + (f"dab {cur_dab:5.2f} |" if config.DAB else "")
                + (f"adv_E {cur_adv_E:5.2f} |" if config.ADV else "")
                + (f"adv_D {cur_adv_D:5.2f} |" if config.ADV else "")
                + (f"nzlp {cur_zero_log_prob:5.2f} |" if config.explicit_zero_prob else "")
                + (
                    f"mvc_nzlp {cur_mvc_zero_log_prob:5.2f} |"
                    if config.MVC and config.explicit_zero_prob
                    else ""
                )
            )
            total_loss = 0
            total_mse = 0
            total_cls = 0
            total_cce = 0
            total_mvc = 0
            total_ecs = 0
            total_dab = 0
            total_adv_E = 0
            total_adv_D = 0
            total_zero_log_prob = 0
            total_mvc_zero_log_prob = 0
            total_error = 0
            start_time = time.time()