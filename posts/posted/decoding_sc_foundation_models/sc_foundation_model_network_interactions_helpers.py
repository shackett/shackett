"""Helpers used by sc_foundation_model_network_interactions.qmd."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import torch

from napistu.constants import BQB_DEFINING_ATTRS_LOOSE, SBML_DFS
from napistu.ontologies.constants import ONTOLOGIES
from napistu.network.constants import NAPISTU_GRAPH_EDGES, NAPISTU_GRAPH_VERTICES
from napistu_torch.foundation_models.attention_patterns import AttentionPatternsInputs
from napistu_torch.foundation_models.constants import FM_EDGELIST
from napistu_torch.foundation_models.foundation_models import FoundationModelStore, FoundationModels
from napistu_torch.utils.tensor_utils import compute_cosine_distances_torch

from sc_foundation_model_utils import (
    append_missing_optional_model_rows,
    get_model_label_maps,
    get_model_names,
    mark_optional_model_inclusion,
    n_genes_in_dataset_embedding_summary,
    to_filename,
)

logger = logging.getLogger(__name__)


def get_cache_path(
    dataset: str,
    category: str,
    cache_dir: Path,
    include_scgpt: bool,
    include_scfoundation: bool,
) -> Path:
    scgpt_flag = "with_scgpt" if include_scgpt else "without_scgpt"
    sf_flag = "with_scfoundation" if include_scfoundation else "without_scfoundation"
    return cache_dir / f"cross_model_{dataset}_{to_filename(category)}_{scgpt_flag}_{sf_flag}.pkl"


def get_all_categories(
    model_outputs_dir: Path,
    embedding_dataset: str,
    ignore_categories_with: List[str],
    include_scgpt: bool,
    include_scfoundation: bool,
    ignored_models: List[str],
) -> List[str]:
    model_names = get_model_names(
        include_scgpt=include_scgpt,
        include_scfoundation=include_scfoundation,
        ignored_models=ignored_models,
    )
    stores = [FoundationModelStore(Path(model_outputs_dir) / name) for name in model_names]
    category_sets = [set(store.list_categories(embedding_dataset)) for store in stores]
    common = category_sets[0].intersection(*category_sets[1:])
    return sorted([
        category
        for category in common
        if not any(substring in category for substring in ignore_categories_with)
    ])


def get_display_and_active_model_metadata(
    model_outputs_dir: Path,
    embedding_dataset: str,
    ignore_categories_with: List[str],
    include_scgpt: bool,
    include_scfoundation: bool,
    ignored_models: Optional[List[str]] = None,
) -> Tuple[pd.DataFrame, dict]:
    """Build the full metadata table and metadata for the active model subset.

    The display table always lists scGPT and scFoundation alongside the core models,
    adding placeholder rows when optional weights are absent and marking which
    optional families are included in the current analysis.
    """
    full_metadata = get_model_comparison_metadata(
        model_outputs_dir=model_outputs_dir,
        embedding_dataset=embedding_dataset,
        ignore_categories_with=ignore_categories_with,
        include_scgpt=True,
        include_scfoundation=True,
        ignored_models=ignored_models,
    )
    model_metadata_summary_full = mark_optional_model_inclusion(
        append_missing_optional_model_rows(full_metadata["model_metadata_summary"]),
        include_scgpt=include_scgpt,
        include_scfoundation=include_scfoundation,
    )

    if include_scgpt and include_scfoundation:
        active_metadata = full_metadata
    else:
        active_metadata = get_model_comparison_metadata(
            model_outputs_dir=model_outputs_dir,
            embedding_dataset=embedding_dataset,
            ignore_categories_with=ignore_categories_with,
            include_scgpt=include_scgpt,
            include_scfoundation=include_scfoundation,
            ignored_models=ignored_models,
        )

    return model_metadata_summary_full, active_metadata


def get_model_comparison_metadata(
    model_outputs_dir: Path,
    embedding_dataset: str,
    ignore_categories_with: List[str],
    include_scgpt: bool = True,
    include_scfoundation: bool = True,
    ignored_models: Optional[List[str]] = None,
) -> dict:
    model_comparison_metadata = dict()

    models = FoundationModels.load_multiple(
        model_outputs_dir,
        get_model_names(
            include_scgpt=include_scgpt,
            include_scfoundation=include_scfoundation,
            ignored_models=ignored_models,
        ),
        verbose=False,
    )
    model_comparison_metadata["model_order"] = models.model_names

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
    model_comparison_metadata["model_metadata_summary"] = model_metadata_summary

    a_category = models.models[0].store.list_categories(embedding_dataset)[0]
    attended_embeddings = AttentionPatternsInputs.from_expression(
        models, embedding_dataset, a_category, verbose=False
    )
    model_comparison_metadata["gene_ids"] = attended_embeddings.common_gene_ids

    return model_comparison_metadata


def get_model_layer_labels(cross_model_top_attentions: pd.DataFrame) -> pd.DataFrame:
    model_layers_df = cross_model_top_attentions[[FM_EDGELIST.MODEL, FM_EDGELIST.LAYER]].drop_duplicates()
    model_layers_df["label"] = model_layers_df.apply(
        lambda row: f"{row[FM_EDGELIST.MODEL]}-{row[FM_EDGELIST.LAYER]}", axis=1
    )
    return model_layers_df


def summarize_cross_model_attention_coherence(
    model_x_layer_rank_agreement: pd.DataFrame,
    model_sources: pd.DataFrame,
) -> pd.Series:
    return (
        model_x_layer_rank_agreement
        .merge(
            model_sources.rename(columns={"model": "query_model", "model_category": "query_model_category"}),
            on=["query_model"], how="left",
        )
        .merge(
            model_sources.rename(columns={"model": "eval_model", "model_category": "eval_model_category"}),
            on=["eval_model"], how="left",
        )
        .query("query_model_category != eval_model_category")
        .assign(log_quantile=lambda df: -np.log10(df["median_quantile"].clip(upper=0.5) * 2))
        .groupby(["query_model", "query_layer", "category"])["log_quantile"]
        .mean()
        .sort_values(ascending=False)
    )


def add_vertex_names_to_edgelist(
    edgelist: pd.DataFrame,
    gene_to_vertex_map: pd.DataFrame,
) -> pd.DataFrame:
    top_k_with_ids = (
        edgelist
        .merge(
            gene_to_vertex_map
            .rename(columns={ONTOLOGIES.ENSEMBL_GENE: FM_EDGELIST.FROM_GENE, NAPISTU_GRAPH_VERTICES.NAME: "from_vertex"})
            .drop(columns=[SBML_DFS.S_ID]),
            on=FM_EDGELIST.FROM_GENE, how="left",
        )
        .merge(
            gene_to_vertex_map
            .rename(columns={ONTOLOGIES.ENSEMBL_GENE: FM_EDGELIST.TO_GENE, NAPISTU_GRAPH_VERTICES.NAME: "to_vertex"})
            .drop(columns=[SBML_DFS.S_ID]),
            on=FM_EDGELIST.TO_GENE, how="left",
        )
    )
    invalid_edges = top_k_with_ids.isna().any(axis=1)
    if invalid_edges.sum() > 0:
        percent_invalid = (invalid_edges.sum() / len(top_k_with_ids)) * 100
        logger.warning(
            f"Dropping {invalid_edges.sum()} edges ({percent_invalid:.2f}%) which could not be mapped to vertices"
        )
    return top_k_with_ids.loc[~invalid_edges]


def calculate_background_edgelist_metrics(
    napistu_data,
    edge_prediction_task,
    vertex_names: List[str],
) -> Tuple[float, float]:
    gene_vertices = set(napistu_data.get_vertex_indices(vertex_names))
    src = napistu_data.edge_index[0]
    dst = napistu_data.edge_index[1]
    gene_vertices_tensor = torch.tensor(list(gene_vertices), dtype=torch.long)
    src_mask = torch.isin(src, gene_vertices_tensor)
    dst_mask = torch.isin(dst, gene_vertices_tensor)
    num_edges = (src_mask & dst_mask).sum().item()
    background_edge_rate = num_edges / (len(gene_vertices) - 1) ** 2

    all_pairs = pd.MultiIndex.from_product(
        [vertex_names, vertex_names], names=[NAPISTU_GRAPH_EDGES.FROM, NAPISTU_GRAPH_EDGES.TO]
    ).to_frame(index=False)
    all_pairs = all_pairs[all_pairs[NAPISTU_GRAPH_EDGES.FROM] != all_pairs[NAPISTU_GRAPH_EDGES.TO]]
    all_pairs = all_pairs.sample(min(10000, len(all_pairs)))

    background_edge_score = edge_prediction_task.predict_edge_scores(
        napistu_data,
        napistu_data.get_edge_indices(all_pairs, NAPISTU_GRAPH_EDGES.FROM, NAPISTU_GRAPH_EDGES.TO),
    ).mean().numpy()

    return background_edge_rate, background_edge_score


def compare_attention_and_napistu_graphs(
    top_k_attention: pd.DataFrame,
    top_k_attention_edgelist: pd.DataFrame,
    napistu_data,
    edge_prediction_task,
    model_order: List[str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    edge_indices = napistu_data.get_edge_indices(top_k_attention_edgelist, "from_vertex", "to_vertex")
    predictions = edge_prediction_task.predict_edge_scores(data=napistu_data, edge_index=edge_indices)

    top_k_attention_edgelist["prediction"] = predictions
    top_k_attention_edgelist["direct_edge_exists"] = napistu_data.has_edges(edge_indices)

    top_k_attention_w_metrics = (
        top_k_attention.reset_index()
        .merge(
            top_k_attention_edgelist[["from_gene", "to_gene", "prediction", "direct_edge_exists"]],
            on=["from_gene", "to_gene"], how="left",
        )
    )

    direct_edge_rate = (
        top_k_attention_w_metrics
        .value_counts(["model", "layer", "category", "direct_edge_exists"])
        .reset_index()
        .pivot(index=["model", "layer", "category"], columns="direct_edge_exists", values="count")
        .fillna(0)
        .assign(true_fraction=lambda df: df[True] / df.sum(axis=1))
        .reset_index()
        .assign(model_name=lambda df: pd.Categorical(df["model"], categories=model_order, ordered=True))
    )

    average_edge_prediction = (
        top_k_attention_w_metrics
        .groupby(["model", "layer", "category"])["prediction"]
        .median()
        .reset_index()
        .assign(model_name=lambda df: pd.Categorical(df["model"], categories=model_order, ordered=True))
    )

    return direct_edge_rate, average_edge_prediction


def plot_metric_by_layer(
    data,
    y_col,
    background_value,
    ylabel,
    model_order,
    model_metadata_summary,
    figsize=(14, 4),
    bar_width=0.6,
    group_gap=2,
):
    model_type_map, model_variant_map = get_model_label_maps(model_metadata_summary)
    fig, ax = plt.subplots(figsize=figsize)

    x_pos_map = {}
    current_x = 0
    model_label_positions = {}

    for model in model_order:
        model_data = data[data["model_name"] == model]
        layers = sorted(model_data["layer"].unique())
        model_start = current_x
        for layer in layers:
            x_pos_map[(model, layer)] = current_x
            current_x += 1
        model_label_positions[model] = (model_start + current_x - 1) / 2
        current_x += group_gap

    for model in model_order:
        model_data = data[data["model_name"] == model]
        for layer, group in model_data.groupby("layer"):
            x = x_pos_map[(model, layer)]
            mean_val = group[y_col].mean()
            min_val = group[y_col].min()
            max_val = group[y_col].max()
            ax.bar(x, mean_val, width=bar_width, color="#aaaaaa", alpha=0.8)
            ax.plot([x, x], [min_val, max_val], color="black", linewidth=1, zorder=5)
            ax.plot([x - 0.15, x + 0.15], [min_val, min_val], color="black", linewidth=1, zorder=5)
            ax.plot([x - 0.15, x + 0.15], [max_val, max_val], color="black", linewidth=1, zorder=5)

    ax.set_xticks([x_pos_map[(model, layer)] for model in model_order
                   for layer in sorted(data[data["model_name"] == model]["layer"].unique())])
    ax.set_xticklabels([str(layer) for model in model_order
                        for layer in sorted(data[data["model_name"] == model]["layer"].unique())],
                       fontsize=7)
    ax.axhline(background_value, color="gray", linestyle=":", linewidth=1.5, zorder=0)
    ax.set_xlabel("Layer")
    ax.set_ylabel(ylabel)
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
            ax.text(x_center, y_variant, variant, ha="center", va="bottom", fontsize=8, color="#555555")

    for model_type, positions in seen_types.items():
        ax.text(np.mean(positions), y_type, model_type, ha="center", va="bottom", fontsize=10, fontweight="bold")

    return fig, ax


def plot_stacked_histogram(ax, df, score_col, x_min, x_max, n_bins=40):
    bin_edges = np.linspace(x_min, x_max, n_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    width = bin_edges[1] - bin_edges[0]

    no_edge = df[~df["direct_edge_exists"]][score_col]
    has_edge = df[df["direct_edge_exists"]][score_col]
    no_edge_counts, _ = np.histogram(no_edge, bins=bin_edges)
    has_edge_counts, _ = np.histogram(has_edge, bins=bin_edges)

    ax.bar(bin_centers, has_edge_counts, width=width, label="Direct edge", alpha=0.8, color="#444444")
    ax.bar(bin_centers, no_edge_counts, width=width, bottom=has_edge_counts, label="No direct edge", alpha=0.8, color="#cccccc")

    medians = [(False, "No edge median", "--"), (True, "Direct edge median", "-.")]
    meds = [(df[df["direct_edge_exists"] == exists][score_col].median(), label, linestyle)
            for exists, label, linestyle in medians]
    meds_sorted = sorted(meds, key=lambda item: item[0])
    y_positions = [0.97, 0.83] if len(meds_sorted) == 2 else [0.97]

    for (med, label, linestyle), y_pos in zip(meds_sorted, y_positions):
        ax.axvline(med, linestyle=linestyle, linewidth=1.5, color="gray")
        ax.text(med + 0.01, y_pos, f"{label}: {med:.2f}",
                transform=ax.get_xaxis_transform(), fontsize=8, va="top")

    ax.set_xlabel(score_col.capitalize())
    ax.set_ylabel("Count")
    ax.legend(fontsize=8)
    sns.despine(ax=ax)


def summarize_edgelist_cosine_similarity(
    top_k_attention: pd.DataFrame,
    top_k_attention_edgelist: pd.DataFrame,
    napistu_data,
    napistu_gnn,
    model_order: List[str],
) -> Tuple[pd.DataFrame, float]:
    vertex_names = pd.concat([
        top_k_attention_edgelist["from_vertex"],
        top_k_attention_edgelist["to_vertex"],
    ]).unique().tolist()
    vertex_indices = napistu_data.get_vertex_indices(vertex_names)
    cosine_sim = 1 - compute_cosine_distances_torch(napistu_gnn.get_embeddings(napistu_data)[vertex_indices])
    vertex_to_idx = {name: i for i, name in enumerate(napistu_data.get_vertex_names()[vertex_indices])}

    from_idx = top_k_attention_edgelist["from_vertex"].map(vertex_to_idx)
    to_idx = top_k_attention_edgelist["to_vertex"].map(vertex_to_idx)

    working = top_k_attention_edgelist.copy()
    working["cosine_sim"] = cosine_sim[from_idx.values, to_idx.values]

    average_cosine_sim = (
        top_k_attention.reset_index()
        .merge(working, on=["from_gene", "to_gene"], how="inner")
        .groupby(["model", "layer", "category"])["cosine_sim"]
        .median()
        .reset_index()
        .assign(model_name=lambda df: pd.Categorical(df["model"], categories=model_order, ordered=True))
    )

    return average_cosine_sim, cosine_sim.mean()


def top_attention_to_napistu_edgelist(
    cross_model_top_attentions: pd.DataFrame,
    gene_to_vertex_map: pd.DataFrame,
    top_k: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    top_k_attention = (
        cross_model_top_attentions
        .query("attention_rank <= @top_k")
        .set_index(["model", "layer"])
        .sort_values("attention_rank")
        .sort_index()
    )
    top_k_attention_distinct_edges = top_k_attention[["from_gene", "to_gene"]].drop_duplicates()
    top_k_attention_edgelist = add_vertex_names_to_edgelist(top_k_attention_distinct_edges, gene_to_vertex_map)
    return top_k_attention, top_k_attention_edgelist
