"""Helpers used by sc_foundation_model_overview.qmd."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from napistu.utils import load_pickle
from napistu_torch.foundation_models.attention_patterns import LayerwiseAttentionInputs
from napistu_torch.foundation_models.constants import COMPARE_EMBEDDINGS_COMPARISONS
from napistu_torch.foundation_models.foundation_models import (
    FoundationModel,
    FoundationModelStore,
    FoundationModels,
)
from napistu_torch.foundation_models.gene_embeddings import GeneEmbeddingsSet

from sc_foundation_model_utils import get_model_label_maps, get_model_names, n_genes_in_dataset_embedding_summary, to_filename


def get_model_categories(
    model_name: str,
    model_outputs_dir: Path,
    embedding_dataset: str,
    ignore_categories_with: List[str],
) -> List[str]:
    store = FoundationModelStore(Path(model_outputs_dir) / model_name)
    return sorted([
        category
        for category in store.list_categories(embedding_dataset)
        if not any(substring in category for substring in ignore_categories_with)
    ])


def get_cache_path(model_name: str, category: str, cache_dir: Path) -> Path:
    return cache_dir / f"within_model_{to_filename(model_name)}_{to_filename(category)}.pkl"


def get_model_comparison_metadata(
    model_outputs_dir: Path,
    embedding_dataset: str,
    ignore_categories_with: List[str],
    ignored_models: List[str],
) -> Dict:
    model_names = get_model_names(ignored_models=ignored_models)
    models = FoundationModels.load_multiple(model_outputs_dir, model_names, verbose=False)

    model_metadata_summary = models.get_summary().rename(columns={
        "full_name": "model",
        "model": "model type",
        "variant": "model variant",
        "embed_dim": "# embedding dim",
        "n_layers": "# layers",
        "n_heads": "# heads",
        "parameter_count": "# parameters",
    })

    model_metadata_summary["# genes in dataset embedding"] = [
        n_genes_in_dataset_embedding_summary(
            model,
            model_outputs_dir=model_outputs_dir,
            embedding_dataset=embedding_dataset,
            ignore_categories_with=ignore_categories_with,
        )
        for model in models.models
    ]

    disk_name_by_full_name = {
        model.full_name: model.disk_name for model in models.models
    }

    return {
        "model_order": models.model_names,
        "disk_name_by_full_name": disk_name_by_full_name,
        "model_metadata_summary": model_metadata_summary,
    }


def model_facet_mosaic_layout(n: int) -> Tuple[List[List[str]], List[float]]:
    """Layout for one panel per model (~half the models per row)."""
    if n < 1:
        raise ValueError("n must be >= 1")
    n_cols = max(1, int(np.ceil(n / 2)))
    n_rows = int(np.ceil(n / n_cols))
    layout: List[List[str]] = []
    k = 0
    for _ in range(n_rows):
        row: List[str] = []
        for _col in range(n_cols):
            row.append(str(k) if k < n else ".")
            k += 1
        layout.append(row)
    width_ratios = [1.0] * n_cols
    return layout, width_ratios


def residual_corr_to_matrix(residual_corr: Dict[str, float], n_layers: int) -> pd.DataFrame:
    """Reshape pairwise layer correlation dict into an n_layers × n_layers DataFrame."""
    mat = np.full((n_layers, n_layers), np.nan)
    np.fill_diagonal(mat, 1.0)
    for key, rho in residual_corr.items():
        indices = [int(x) for x in re.findall(r"(?<=layer_)\d+", key)]
        if len(indices) == 2:
            i, j = indices
            mat[i, j] = rho
            mat[j, i] = rho
    return pd.DataFrame(mat)


def compute_within_model_layer_comparisons(
    model: FoundationModel,
    dataset_name: str,
    category: str,
    top_k: int,
    k_sensitivity_values: List[int],
    by_absolute_value: bool = False,
    ignore_self_attention: bool = True,
    verbose: bool = False,
) -> Dict:
    if top_k not in k_sensitivity_values:
        raise ValueError(f"top_k={top_k} must be in k_sensitivity_values={k_sensitivity_values}")

    residuals = model.load_category_residuals(dataset_name, category)
    lwa = LayerwiseAttentionInputs(residual_stream_embeddings=residuals, foundation_model=model)

    layer_rank_by_k = {}
    for k in k_sensitivity_values:
        layer_corr_k, layer_rank_k = lwa.compare_layer_attention_consistency(
            top_k=k,
            by_absolute_value=by_absolute_value,
            ignore_self_attention=ignore_self_attention,
            verbose=verbose,
        )
        layer_rank_by_k[k] = layer_rank_k
        if k == top_k:
            layer_corr = layer_corr_k
            layer_rank = layer_rank_k

    embedding_set = GeneEmbeddingsSet.from_gene_embeddings(
        [residuals[i] for i in sorted(residuals.keys())]
    )
    residual_stream_layer_corr = embedding_set.compare_embeddings()

    return {
        COMPARE_EMBEDDINGS_COMPARISONS.MODEL_LAYER_CORRELATIONS: layer_corr,
        COMPARE_EMBEDDINGS_COMPARISONS.MODEL_LAYER_RANK_AGREEMENT: layer_rank,
        "layer_rank_by_k": layer_rank_by_k,
        "residual_stream_layer_correlations": residual_stream_layer_corr,
    }


def aggregate_model_categories(
    model_name: str,
    k_sensitivity_values: List[int],
    cache_dir: Path,
    model_outputs_dir: Path,
    embedding_dataset: str,
    ignore_categories_with: List[str],
) -> Dict:
    categories = get_model_categories(
        model_name,
        model_outputs_dir=model_outputs_dir,
        embedding_dataset=embedding_dataset,
        ignore_categories_with=ignore_categories_with,
    )

    cat_layer_corrs: List[np.ndarray] = []
    cat_rank_dfs: List[pd.DataFrame] = []
    cat_rank_by_k: Dict[int, List[pd.DataFrame]] = {k: [] for k in k_sensitivity_values}
    cat_residual_pairs: Dict[str, List[float]] = {}

    for category in categories:
        data = load_pickle(get_cache_path(model_name, category, cache_dir=cache_dir))
        cat_layer_corrs.append(data[COMPARE_EMBEDDINGS_COMPARISONS.MODEL_LAYER_CORRELATIONS])
        cat_rank_dfs.append(data[COMPARE_EMBEDDINGS_COMPARISONS.MODEL_LAYER_RANK_AGREEMENT])
        for k in k_sensitivity_values:
            cat_rank_by_k[k].append(data["layer_rank_by_k"][k])
        for pair, rho in data["residual_stream_layer_correlations"].items():
            cat_residual_pairs.setdefault(pair, []).append(rho)

    return {
        COMPARE_EMBEDDINGS_COMPARISONS.MODEL_LAYER_CORRELATIONS: np.median(cat_layer_corrs, axis=0),
        COMPARE_EMBEDDINGS_COMPARISONS.MODEL_LAYER_RANK_AGREEMENT: (
            pd.concat(cat_rank_dfs)
            .groupby(["query_layer", "eval_layer"])["median_quantile"]
            .median()
            .reset_index()
        ),
        "layer_rank_by_k": {
            k: (
                pd.concat(dfs)
                .groupby(["query_layer", "eval_layer"])["median_quantile"]
                .median()
                .reset_index()
            )
            for k, dfs in cat_rank_by_k.items()
        },
        "residual_stream_layer_correlations": {
            pair: float(np.median(rhos)) for pair, rhos in cat_residual_pairs.items()
        },
    }


def plot_rank_agreement_by_layer(
    layer_rank_by_k,
    K_values,
    model_order,
    model_metadata_summary,
    figsize=(14, 4),
    group_gap=2,
    reference_value=0.5,
):
    k_colors = ["#E63946", "#2196F3", "#FF9800", "#4CAF50"]
    color_map = dict(zip(sorted(K_values), k_colors))

    model_type_map, model_variant_map = get_model_label_maps(model_metadata_summary)
    fig, ax = plt.subplots(figsize=figsize)

    x_pos_map = {}
    current_x = 0
    model_label_positions = {}

    for model in model_order:
        sample_k = sorted(K_values)[0]
        rank_df = layer_rank_by_k[model][sample_k]
        consecutive = rank_df[rank_df["eval_layer"] == rank_df["query_layer"] + 1]
        query_layers = sorted(consecutive["query_layer"].unique())

        model_start = current_x
        for layer in query_layers:
            x_pos_map[(model, layer)] = current_x
            current_x += 1
        model_label_positions[model] = (model_start + current_x - 1) / 2
        current_x += group_gap

    for model in model_order:
        for k in sorted(K_values):
            rank_df = layer_rank_by_k[model][k]
            consecutive = (
                rank_df[rank_df["eval_layer"] == rank_df["query_layer"] + 1]
                .sort_values("query_layer")
                .assign(median_quantile=lambda x: 1 - x["median_quantile"])
            )
            xs = [
                x_pos_map[(model, ql)]
                for ql in consecutive["query_layer"]
                if (model, ql) in x_pos_map
            ]
            ys = [
                row["median_quantile"]
                for _, row in consecutive.iterrows()
                if (model, row["query_layer"]) in x_pos_map
            ]
            ax.plot(
                xs,
                ys,
                marker="o",
                markersize=4,
                color=color_map[k],
                label=f"K = {k:,}",
                linewidth=1.5,
            )

    all_xticks = []
    all_xlabels = []
    for model in model_order:
        sample_k = sorted(K_values)[0]
        rank_df = layer_rank_by_k[model][sample_k]
        consecutive = rank_df[rank_df["eval_layer"] == rank_df["query_layer"] + 1]
        for layer in sorted(consecutive["query_layer"].unique()):
            if (model, layer) in x_pos_map:
                all_xticks.append(x_pos_map[(model, layer)])
                all_xlabels.append(str(layer))

    ax.set_xticks(all_xticks)
    ax.set_xticklabels(all_xlabels, fontsize=7)
    ax.axhline(reference_value, color="gray", linestyle=":", linewidth=1.5, zorder=0)
    ax.set_xlabel("Layer")
    ax.set_ylabel("Median quantile of top-K pairs\nin next layer")
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()

    y_max = ax.get_ylim()[1]
    y_type = y_max * 1.12
    y_variant = y_max * 1.04
    seen_types = {}

    for model, x_center in model_label_positions.items():
        model_type = model_type_map.get(model, model)
        variant = model_variant_map.get(model)
        if model_type not in seen_types:
            seen_types[model_type] = []
        seen_types[model_type].append(x_center)
        if pd.notna(variant):
            ax.text(
                x_center,
                y_variant,
                variant,
                ha="center",
                va="bottom",
                fontsize=8,
                color="#555555",
            )

    for model_type, positions in seen_types.items():
        ax.text(
            np.mean(positions),
            y_type,
            model_type,
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )

    handles = [
        plt.Line2D([0], [0], color=color_map[k], marker="o", label=f"K = {k:,}")
        for k in sorted(K_values)
    ]
    ax.legend(
        handles=handles,
        bbox_to_anchor=(1.01, 0.5),
        loc="center left",
        frameon=False,
        fontsize=9,
    )

    plt.suptitle(
        "Consecutive-layer rank consistency is stable across K",
        fontweight="bold",
        fontsize=11,
        y=1.15,
    )

    return fig, ax
