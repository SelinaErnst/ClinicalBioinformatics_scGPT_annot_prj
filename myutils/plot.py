from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc


def inspect_prediction(adata, celltype: str, palette: dict, save_dir: Path | None = None):
    adata_plot = adata.copy()
    if celltype in adata_plot.obs["celltype"].cat.categories:
        prediction_col = f"{celltype} prediction"
        adata_plot.obs[prediction_col] = (
            adata_plot.obs["predictions"]
            .astype(str)
            .where(
                adata_plot.obs["celltype"] == celltype,
                "other",
            )
        )
        mask = adata_plot.obs["celltype"] == celltype
        # Add grey for background
        prediction_palette = palette.copy()
        prediction_palette["other"] = "lightgrey"
        fig, ax = plt.subplots(figsize=(7, 7))
        sc.pl.umap(
            adata_plot[~mask],
            color=prediction_col,
            palette={"other": "lightgrey"},
            size=15,
            ax=ax,
            show=False,
        )
        sc.pl.umap(
            adata_plot[mask],
            color=prediction_col,
            palette=prediction_palette,
            size=50,
            ax=ax,
            show=False,
        )
        if save_dir != None:
            plt.savefig(
                save_dir / f"umap_prediction_{celltype.replace(' ', '_')}.png",
                dpi=300,
                bbox_inches="tight",
            )



def select_celltype(adata, celltype: str, palette: dict | None = None, save_dir: Path | None = None):
    adata_plot = adata.copy()
    if (palette == None):
        if "celltype_colors" not in adata.uns: return
        palette = dict(
            zip(
                adata.obs["celltype"].cat.categories,
                adata.uns["celltype_colors"],
            )
        )
    if celltype in adata_plot.obs["celltype"].cat.categories:
        prediction_col = f"celltype {celltype}"
        adata_plot.obs[prediction_col] = (
            adata_plot.obs["celltype"]
            .astype(str)
            .where(
                adata_plot.obs["celltype"] == celltype,
                "other",
            )
        )
        mask = adata_plot.obs["celltype"] == celltype
        # Add grey for background
        prediction_palette = palette.copy()
        prediction_palette["other"] = "lightgrey"
        fig, ax = plt.subplots(figsize=(7, 7))
        sc.pl.umap(
            adata_plot[~mask],
            color=prediction_col,
            palette={"other": "lightgrey"},
            size=15,
            ax=ax,
            show=False,
        )
        sc.pl.umap(
            adata_plot[mask],
            color=prediction_col,
            palette=prediction_palette,
            size=50,
            ax=ax,
            show=False,
        )
        if save_dir != None:
            plt.savefig(
                save_dir / f"umap_prediction_{celltype.replace(' ', '_')}.png",
                dpi=300,
                bbox_inches="tight",
            )


def select_prediction(adata, celltype: str, palette: dict, save_dir: Path | None = None):
    adata_plot = adata.copy()
    if celltype in adata_plot.obs["predictions"].cat.categories:
        prediction_col = f"predicted {celltype}"
        adata_plot.obs[prediction_col] = (
            adata_plot.obs["predictions"]
            .astype(str)
            .where(
                adata_plot.obs["predictions"] == celltype,
                "other",
            )
        )
        mask = adata_plot.obs["predictions"] == celltype
        # Add grey for background
        prediction_palette = palette.copy()
        prediction_palette["other"] = "lightgrey"
        fig, ax = plt.subplots(figsize=(7, 7))
        sc.pl.umap(
            adata_plot[~mask],
            color=prediction_col,
            palette={"other": "lightgrey"},
            size=15,
            ax=ax,
            show=False,
        )
        sc.pl.umap(
            adata_plot[mask],
            color=prediction_col,
            palette=prediction_palette,
            size=50,
            ax=ax,
            show=False,
        )
        if save_dir != None:
            plt.savefig(
                save_dir / f"umap_prediction_{celltype.replace(' ', '_')}.png",
                dpi=300,
                bbox_inches="tight",
            )