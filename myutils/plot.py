from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc


def inspect_prediction(adata, celltype: str, palette: dict, save_dir: Path | None = None, celltype_col="celltype", prediction_col="predictions", show=True):
    adata_plot = adata.copy()

    if celltype_col not in adata_plot.obs.columns: return
    if prediction_col not in adata_plot.obs.columns: return
    
    adata_plot.obs[celltype_col] = adata_plot.obs[celltype_col].astype(str).astype('category')
    adata_plot.obs[prediction_col] = adata_plot.obs[prediction_col].astype(str).astype('category')
    
    if celltype in adata_plot.obs[celltype_col].cat.categories:
        plot_col = f"{celltype} prediction"
        
        adata_plot.obs[plot_col] = (
            adata_plot.obs[prediction_col]  
            .astype(str)
            .where(
                adata_plot.obs[celltype_col] == celltype,
                "other",
            )
        )
        
        mask = adata_plot.obs[celltype_col] == celltype
        
        adata_plot.obs[plot_col] = adata_plot.obs[plot_col].astype('category')
        
        prediction_palette = {cat: palette.get(cat, "#7f7f7f") for cat in adata_plot.obs[plot_col].cat.categories}
        prediction_palette["other"] = "lightgrey"
        
        fig, ax = plt.subplots(figsize=(7, 7))
        
        sc.pl.umap(
            adata_plot[~mask],
            color=plot_col,
            palette={"other": "lightgrey"},
            size=15,
            ax=ax,
            show=False,
            sort_order=False
        )
        
        sc.pl.umap(
            adata_plot[mask],
            color=plot_col,
            palette=prediction_palette,
            size=40,
            ax=ax,
            show=False,
            sort_order=False
        )
        
        if save_dir is not None:
            plt.savefig(
                save_dir / f"umap_prediction_{celltype.replace(' ', '_').replace('/','_')}.png",
                dpi=300,
                bbox_inches="tight",
            )
        if show: plt.show()


def select_celltype(adata, celltype: str, palette: dict | None = None, save_dir: Path | None = None, celltype_col="celltype"):
    adata_plot = adata.copy()
    if celltype_col not in adata_plot.obs.columns:
        return
        
    adata_plot.obs[celltype_col] = adata_plot.obs[celltype_col].astype(str).astype('category')

    if palette is None:
        if f"{celltype_col}_colors" not in adata_plot.uns:
            return
        palette = dict(zip(adata_plot.obs[celltype_col].cat.categories, adata_plot.uns[f"{celltype_col}_colors"]))

    if celltype in adata_plot.obs[celltype_col].cat.categories:
        plot_col = f"celltype {celltype}"
        
        adata_plot.obs[plot_col] = (
            adata_plot.obs[celltype_col]
            .astype(str)
            .where(adata_plot.obs[celltype_col] == celltype, "other")
        )
        adata_plot.obs[plot_col] = adata_plot.obs[plot_col].astype('category')
        mask = adata_plot.obs[celltype_col] == celltype

        prediction_palette = {cat: palette.get(cat, "#7f7f7f") for cat in adata_plot.obs[plot_col].cat.categories}
        prediction_palette["other"] = "lightgrey"

        fig, ax = plt.subplots(figsize=(7, 7))
        
        sc.pl.umap(
            adata_plot[~mask], color=plot_col, palette={"other": "lightgrey"}, 
            size=15, ax=ax, show=False, sort_order=False
        )
        sc.pl.umap(
            adata_plot[mask], color=plot_col, palette=prediction_palette, 
            size=50, ax=ax, show=False, sort_order=False
        )

        if save_dir is not None:
            safe_name = celltype.replace(' ', '_').replace('/', '_')
            plt.savefig(save_dir / f"umap_celltype_{safe_name}.png", dpi=300, bbox_inches="tight")
        
        plt.show()

def select_prediction(adata, celltype: str, palette: dict, save_dir: Path | None = None, prediction_col="predictions"):
    adata_plot = adata.copy()
    if prediction_col not in adata_plot.obs.columns:
        return
        
    adata_plot.obs[prediction_col] = adata_plot.obs[prediction_col].astype(str).astype('category')

    if celltype in adata_plot.obs[prediction_col].cat.categories:
        plot_col = f"predicted {celltype}"
        
        adata_plot.obs[plot_col] = (
            adata_plot.obs[prediction_col]
            .astype(str)
            .where(adata_plot.obs[prediction_col] == celltype, "other")
        )
        adata_plot.obs[plot_col] = adata_plot.obs[plot_col].astype('category')
        mask = adata_plot.obs[prediction_col] == celltype

        prediction_palette = {cat: palette.get(cat, "#7f7f7f") for cat in adata_plot.obs[plot_col].cat.categories}
        prediction_palette["other"] = "lightgrey"

        fig, ax = plt.subplots(figsize=(7, 7))
        
        sc.pl.umap(
            adata_plot[~mask], color=plot_col, palette={"other": "lightgrey"}, 
            size=15, ax=ax, show=False, sort_order=False
        )
        sc.pl.umap(
            adata_plot[mask], color=plot_col, palette=prediction_palette, 
            size=50, ax=ax, show=False, sort_order=False
        )

        if save_dir is not None:
            safe_name = celltype.replace(' ', '_').replace('/', '_')
            plt.savefig(save_dir / f"umap_prediction_{safe_name}.png", dpi=300, bbox_inches="tight")
            
        plt.show()

# def select_celltype(adata, celltype: str, palette: dict | None = None, save_dir: Path | None = None):
#     adata_plot = adata.copy()
#     if (palette == None):
#         if "celltype_colors" not in adata.uns: return
#         palette = dict(
#             zip(
#                 adata.obs["celltype"].cat.categories,
#                 adata.uns["celltype_colors"],
#             )
#         )
#     if celltype in adata_plot.obs["celltype"].cat.categories:
#         prediction_col = f"celltype {celltype}"
#         adata_plot.obs[prediction_col] = (
#             adata_plot.obs["celltype"]
#             .astype(str)
#             .where(
#                 adata_plot.obs["celltype"] == celltype,
#                 "other",
#             )
#         )
#         mask = adata_plot.obs["celltype"] == celltype
#         # Add grey for background
#         prediction_palette = palette.copy()
#         prediction_palette["other"] = "lightgrey"
#         fig, ax = plt.subplots(figsize=(7, 7))
#         sc.pl.umap(
#             adata_plot[~mask],
#             color=prediction_col,
#             palette={"other": "lightgrey"},
#             size=15,
#             ax=ax,
#             show=False,
#         )
#         sc.pl.umap(
#             adata_plot[mask],
#             color=prediction_col,
#             palette=prediction_palette,
#             size=50,
#             ax=ax,
#             show=False,
#         )
#         if save_dir != None:
#             plt.savefig(
#                 save_dir / f"umap_prediction_{celltype.replace(' ', '_')}.png",
#                 dpi=300,
#                 bbox_inches="tight",
#             )


# def select_prediction(adata, celltype: str, palette: dict, save_dir: Path | None = None):
#     adata_plot = adata.copy()
#     if celltype in adata_plot.obs["predictions"].cat.categories:
#         prediction_col = f"predicted {celltype}"
#         adata_plot.obs[prediction_col] = (
#             adata_plot.obs["predictions"]
#             .astype(str)
#             .where(
#                 adata_plot.obs["predictions"] == celltype,
#                 "other",
#             )
#         )
#         mask = adata_plot.obs["predictions"] == celltype
#         # Add grey for background
#         prediction_palette = palette.copy()
#         prediction_palette["other"] = "lightgrey"
#         fig, ax = plt.subplots(figsize=(7, 7))
#         sc.pl.umap(
#             adata_plot[~mask],
#             color=prediction_col,
#             palette={"other": "lightgrey"},
#             size=15,
#             ax=ax,
#             show=False,
#         )
#         sc.pl.umap(
#             adata_plot[mask],
#             color=prediction_col,
#             palette=prediction_palette,
#             size=50,
#             ax=ax,
#             show=False,
#         )
#         if save_dir != None:
#             plt.savefig(
#                 save_dir / f"umap_prediction_{celltype.replace(' ', '_')}.png",
#                 dpi=300,
#                 bbox_inches="tight",
#             )