"""Shared helpers for sc foundation model analysis notebooks."""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd

from napistu_torch.foundation_models.constants import FOUNDATION_MODEL_NAMES, MODEL_NICE_NAMES
from napistu_torch.foundation_models.foundation_models import FoundationModel, _get_disk_name

logger = logging.getLogger(__name__)


def to_filename(s: str) -> str:
    s = re.sub(r"[^\w\s-]", "", s)
    s = re.sub(r"\s+", "_", s)
    return s.strip("-_")


def get_model_names(
    include_scgpt: bool = True,
    include_scfoundation: bool = True,
    ignored_models: Optional[List[str]] = None,
) -> List[str]:
    ignored = ignored_models or []
    full_names = [
        _get_disk_name(model_name, model_variant)
        for model_name, model_variant in MODEL_NICE_NAMES
        if include_scgpt or model_name != FOUNDATION_MODEL_NAMES.SCGPT
        if include_scfoundation or model_name != FOUNDATION_MODEL_NAMES.SCFOUNDATION
    ]
    if ignored:
        invalid = [name for name in ignored if name not in full_names]
        if invalid:
            logger.warning(f"Ignored models not found in known model names: {invalid}")
        full_names = [name for name in full_names if name not in ignored]
    return full_names


def n_genes_in_dataset_embedding_summary(
    disk_name_or_model: str | FoundationModel,
    model_outputs_dir: Path,
    embedding_dataset: str,
    ignore_categories_with: List[str],
) -> int | str:
    """Distinct gene-token counts aggregated across analysed categories.

    Most models fix vocabulary size across cell-type clusters on a dataset.
    scFoundation can differ by category (expression-dependent token selection),
    so we report a numeric range ``min–max`` when counts are not identical.
    """
    if isinstance(disk_name_or_model, str):
        model = FoundationModel.load(Path(model_outputs_dir) / disk_name_or_model)
        disk_tag = disk_name_or_model
    else:
        model = disk_name_or_model
        disk_tag = model.disk_name

    categories = sorted([
        category
        for category in model.store.list_categories(embedding_dataset)
        if not any(substring in category for substring in ignore_categories_with)
    ])
    if not categories:
        raise ValueError(f"No categories for model {disk_tag!r} in dataset {embedding_dataset!r}")

    counts: List[int] = []
    for category in categories:
        residuals = model.load_category_residuals(embedding_dataset, category)
        counts.append(next(iter(residuals.values())).n_genes)

    uniq = sorted(set(counts))
    if len(uniq) == 1:
        return uniq[0]
    return f"{min(counts)}\u2013{max(counts)}"


def get_model_label_maps(model_metadata_summary: pd.DataFrame) -> Tuple[dict, dict]:
    model_type_map = model_metadata_summary.set_index("model")["model type"].to_dict()
    model_variant_map = model_metadata_summary.set_index("model")["model variant"].to_dict()
    return model_type_map, model_variant_map


OPTIONAL_MODEL_TYPES = (
    FOUNDATION_MODEL_NAMES.SCGPT,
    FOUNDATION_MODEL_NAMES.SCFOUNDATION,
)


def append_missing_optional_model_rows(
    model_metadata_summary: pd.DataFrame,
    optional_model_types: Tuple[str, ...] = OPTIONAL_MODEL_TYPES,
) -> pd.DataFrame:
    """Add placeholder rows for optional models whose weights were not loaded locally."""
    if "model type" not in model_metadata_summary.columns:
        return model_metadata_summary

    present_types = set(model_metadata_summary["model type"].astype(str))
    stub_rows = []
    for model_type in optional_model_types:
        if model_type in present_types:
            continue
        label = MODEL_NICE_NAMES.get((model_type, None), model_type)
        row = {column: pd.NA for column in model_metadata_summary.columns}
        row["model"] = label
        row["model type"] = model_type
        stub_rows.append(row)

    if not stub_rows:
        return model_metadata_summary
    return pd.concat([model_metadata_summary, pd.DataFrame(stub_rows)], ignore_index=True)


def mark_optional_model_inclusion(
    model_metadata_summary: pd.DataFrame,
    include_scgpt: bool,
    include_scfoundation: bool,
) -> pd.DataFrame:
    """Annotate whether each optional model family is included in the active analysis."""
    optional_flags = {
        FOUNDATION_MODEL_NAMES.SCGPT: include_scgpt,
        FOUNDATION_MODEL_NAMES.SCFOUNDATION: include_scfoundation,
    }

    def _inclusion_mark(model_type) -> str:
        model_type_str = "" if pd.isna(model_type) else str(model_type).strip()
        if model_type_str in optional_flags:
            return "✅" if optional_flags[model_type_str] else "❌"
        return "✅"

    result = model_metadata_summary.copy()
    result["included"] = result["model type"].apply(_inclusion_mark)
    return result
