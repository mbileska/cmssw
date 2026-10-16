#!/usr/bin/env python3
"""Stream a background NanoAOD production and build a broad QA plot campaign.

The campaign deliberately keeps three concepts separate:

* raw generated rows, which reveal how the input files were sampled;
* within-parent generator-reweighted shapes, which reveal the shape seen by a
  generator-weighted loss without inventing cross-process normalizations; and
* generator-reweighted subcategories, which test whether a parent background
  is internally multimodal after the sampling correction.

Cross sections are required before the four parents can be combined into a
physical rate prediction.  This script therefore normalizes every displayed
parent/category shape independently and never labels the mixture as a yield.
"""

from __future__ import annotations

import argparse
import csv
import html
import itertools
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

try:
    from .backgroundAuditCore import (
        AlignedReservoir,
        FeatureStats,
        PriorityReservoir,
        make_bins,
        normalized_histogram,
        safe_name,
        weighted_correlation,
        weighted_ks,
    )
except ImportError:  # Direct execution: python3 auditBackgroundNano.py
    from backgroundAuditCore import (
        AlignedReservoir,
        FeatureStats,
        PriorityReservoir,
        make_bins,
        normalized_histogram,
        safe_name,
        weighted_correlation,
        weighted_ks,
    )


PROCESS_COLORS = {
    0: "#4c78a8",
    1: "#f58518",
    2: "#54a24b",
    3: "#e45756",
}

CONTROL_BRANCHES = {
    "Background_process",
    "Background_sample",
    "Background_category",
    "Background_kinematicValue",
    "genWeight",
}

LABEL_BRANCHES = {
    "Background_process", "Background_sample", "Background_category",
}

SCORE_FRAGMENTS = ("AXO", "CICADA")

NURD_COLLECTIONS = (
    "ScoutingParticle",
    "ScoutingPFJet",
    "ScoutingPFJetRecluster",
    "ScoutingMuonVtx",
    "ScoutingMuonNoVtx",
    "ScoutingElectron",
    "ScoutingPhoton",
    "ScoutingMET",
)

PHYSICS_COLLECTIONS = NURD_COLLECTIONS + (
    "ScoutingTrack",
    "ScoutingFatPFJetRecluster",
    "ScoutingPrimaryVertex",
    "ScoutingRho",
    "Pileup",
    "L1Mu",
    "L1EG",
    "L1Jet",
    "L1Tau",
    "L1EtSum",
    "L1AXOInput",
    "L1CICADAInput",
)

IDENTIFIER_BRANCHES = {"run", "luminosityBlock", "event"}

OBJECT_SLOT_LIMITS = {
    "ScoutingPFJet": 10,
    "ScoutingPFJetRecluster": 10,
    "ScoutingMuonVtx": 4,
    "ScoutingMuonNoVtx": 4,
    "ScoutingElectron": 4,
    "ScoutingPhoton": 4,
}

MODEL_PDG_IDS = (211, 11, 13, 22, 130, 1, 2)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit all useful kinematics and scores in background NanoAOD files.")
    parser.add_argument(
        "inputs", nargs="+",
        help="ROOT files and/or @text-file lists (one file per line).")
    parser.add_argument(
        "--output-dir", default="production_reports/input_audit",
        help="New or existing directory for plots and machine-readable reports.")
    parser.add_argument(
        "--categories",
        default=str(Path(__file__).with_name("background_categories.json")),
        help="Process/category schema JSON.")
    parser.add_argument(
        "--expected-categories", default="0-23,27-36",
        help=("Required category IDs/ranges for this production, or 'all'. "
              "The default matches train3M_v1, which has no Z-to-nunu parent."))
    parser.add_argument(
        "--profile", choices=("nurd", "physics", "exhaustive"),
        default="exhaustive",
        help="nurd: model-source tables; physics/exhaustive: broader scouting/L1 QA.")
    parser.add_argument("--step-size", default="100 MB")
    parser.add_argument("--event-reservoir", type=int, default=50_000)
    parser.add_argument("--candidate-reservoir", type=int, default=30_000)
    parser.add_argument("--max-files", type=int, default=-1)
    parser.add_argument(
        "--max-plotted-features", type=int, default=-1,
        help="Plot only the first N features for a quick smoke test; still rank all.")
    parser.add_argument("--seed", type=int, default=20251016)
    parser.add_argument(
        "--min-category-entries", type=int, default=20,
        help="Minimum sampled entries required to draw/rank a subcategory.")
    parser.add_argument(
        "--skip-bad-files", action="store_true",
        help="Record unreadable/incompatible files and continue instead of failing.")
    parser.add_argument(
        "--pdf", action="store_true",
        help="Also save every plot as PDF (off by default; this can be large).")
    return parser.parse_args()


def expand_inputs(tokens: list[str]) -> list[str]:
    result: list[str] = []
    for token in tokens:
        if token.startswith("@"):
            path = Path(token[1:])
            for line in path.read_text().splitlines():
                value = line.strip()
                if value and not value.startswith("#"):
                    result.append(value)
        else:
            result.append(token)
    return list(dict.fromkeys(result))


def awkward_ndim(array) -> int:
    ndim = getattr(array, "ndim", None)
    if ndim is not None:
        return int(ndim)
    return int(array.layout.purelist_depth)


def to_numpy(array, awkward) -> np.ndarray:
    return np.asarray(awkward.to_numpy(array))


def flatten_per_event(array, awkward) -> tuple[np.ndarray, np.ndarray]:
    """Flatten any nested numeric branch while retaining scalar counts/event."""
    flattened = array
    while awkward_ndim(flattened) > 2:
        flattened = awkward.flatten(flattened, axis=2)
    counts = to_numpy(awkward.num(flattened, axis=1), awkward).astype(np.int64)
    values = to_numpy(
        awkward.flatten(flattened, axis=None), awkward).astype(np.float64)
    return values, counts


def is_numeric_typename(typename: str) -> bool:
    lowered = typename.lower()
    rejected = ("string", "char*", "tstring")
    return not any(fragment in lowered for fragment in rejected)


def discover_branches(tree, profile: str) -> list[str]:
    names = list(tree.keys())
    collections = NURD_COLLECTIONS if profile == "nurd" else PHYSICS_COLLECTIONS
    selected = []
    for name in names:
        branch = tree[name]
        if not is_numeric_typename(getattr(branch, "typename", "")):
            continue
        include = (
            name in CONTROL_BRANCHES
            or name in IDENTIFIER_BRANCHES
            or any(fragment.lower() in name.lower() for fragment in SCORE_FRAGMENTS)
            or any(name.startswith(f"{collection}_") for collection in collections)
            or any(name == f"n{collection}" for collection in collections)
        )
        if profile == "exhaustive":
            include = include or not name.startswith(("HLT_", "Flag_", "L1_"))
        if include:
            selected.append(name)
    missing = sorted(CONTROL_BRANCHES - set(selected))
    if missing:
        raise RuntimeError(f"required NanoAOD branches are missing: {missing}")
    return sorted(selected)


def load_schema(path: str) -> tuple[dict[int, str], dict[int, dict]]:
    payload = json.loads(Path(path).read_text())
    if "processes" in payload:
        processes = {
            int(key): str(value) for key, value in payload["processes"].items()
        }
        process_ids = None
    else:
        process_ids = {
            str(name): int(process_id)
            for name, process_id in payload["process_ids"].items()
        }
        display_names = {
            "qcd": "QCD", "wjets": "W+jets",
            "zjets": "Z+jets", "ttbar": "ttbar",
        }
        processes = {
            process_id: display_names.get(name, name)
            for name, process_id in process_ids.items()
        }
    categories = {}
    for item in payload["categories"]:
        entry = dict(item)
        if isinstance(entry["process"], str):
            if process_ids is None:
                raise ValueError("string category processes require process_ids")
            entry["process"] = process_ids[entry["process"]]
        else:
            entry["process"] = int(entry["process"])
        categories[int(entry["id"])] = entry
    return processes, categories


def parse_category_spec(specification: str, valid_ids: set[int]) -> set[int]:
    if specification.strip().lower() == "all":
        return set(valid_ids)
    result = set()
    for token in specification.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            lower, upper = map(int, token.split("-", 1))
            if upper < lower:
                raise ValueError(f"invalid descending category range: {token}")
            result.update(range(lower, upper + 1))
        else:
            result.add(int(token))
    unknown = result - valid_ids
    if unknown:
        raise ValueError(f"unknown expected category IDs: {sorted(unknown)}")
    return result


def collection_from_branch(name: str) -> str | None:
    for collection in sorted(PHYSICS_COLLECTIONS, key=len, reverse=True):
        if name.startswith(f"{collection}_"):
            return collection
    return None


def source_collection(feature: str) -> str:
    derived_prefix = feature.split("__", 1)[0]
    if derived_prefix in PHYSICS_COLLECTIONS:
        return derived_prefix
    return collection_from_branch(feature) or "event"


def derived_feature_names(
    jagged_branches: list[str],
) -> dict[str, set[str]]:
    available = set(jagged_branches)
    result: dict[str, set[str]] = {}
    for collection in PHYSICS_COLLECTIONS:
        pt_name = f"{collection}_pt"
        if pt_name not in available:
            continue
        fields = {"n", "sum_pt", "mean_pt", "leading_pt", "subleading_pt"}
        if f"{collection}_eta" in available:
            fields.add("max_abs_eta")
        if f"{collection}_energy" in available:
            fields.add("sum_energy")
        if f"{collection}_mass" in available:
            fields.add("sum_mass")
        result[collection] = fields
    return result


def derive_event_features(
    arrays: dict,
    specifications: dict[str, set[str]],
    awkward,
) -> dict[str, np.ndarray]:
    result: dict[str, np.ndarray] = {}
    for collection, fields in specifications.items():
        pt = arrays[f"{collection}_pt"]
        counts = to_numpy(awkward.num(pt, axis=1), awkward).astype(np.int64)
        result[f"{collection}__n"] = counts.astype(np.float64)
        result[f"{collection}__sum_pt"] = to_numpy(
            awkward.sum(pt, axis=1), awkward)
        result[f"{collection}__mean_pt"] = to_numpy(
            awkward.fill_none(awkward.mean(pt, axis=1), 0.0), awkward)
        sorted_pt = awkward.sort(pt, axis=1, ascending=False)
        order = awkward.argsort(pt, axis=1, ascending=False)
        result[f"{collection}__leading_pt"] = to_numpy(
            awkward.fill_none(awkward.firsts(sorted_pt), 0.0), awkward)
        result[f"{collection}__subleading_pt"] = to_numpy(
            awkward.fill_none(awkward.firsts(sorted_pt[:, 1:]), 0.0), awkward)
        if "max_abs_eta" in fields:
            eta = arrays[f"{collection}_eta"]
            result[f"{collection}__max_abs_eta"] = to_numpy(
                awkward.fill_none(awkward.max(abs(eta), axis=1), 0.0), awkward)
        if "sum_energy" in fields:
            result[f"{collection}__sum_energy"] = to_numpy(
                awkward.sum(arrays[f"{collection}_energy"], axis=1), awkward)
        if "sum_mass" in fields:
            result[f"{collection}__sum_mass"] = to_numpy(
                awkward.sum(arrays[f"{collection}_mass"], axis=1), awkward)
        slot_limit = OBJECT_SLOT_LIMITS.get(collection, 0)
        for field in ("pt", "eta", "phi"):
            branch = f"{collection}_{field}"
            if branch not in arrays or slot_limit <= 0:
                continue
            sorted_values = arrays[branch][order]
            padded = awkward.pad_none(
                sorted_values, slot_limit, axis=1, clip=True)
            for rank in range(slot_limit):
                result[f"{collection}__slot{rank}_{field}"] = to_numpy(
                    awkward.fill_none(padded[:, rank], 0.0), awkward)
    return result


def pf_model_candidate_features(arrays: dict, awkward) -> dict[str, np.ndarray]:
    """Reproduce deterministic PFPreProcessor inputs for direct PF rows."""
    required = {
        "ScoutingParticle_pt", "ScoutingParticle_eta",
        "ScoutingParticle_dxy", "ScoutingParticle_dxysig",
        "ScoutingParticle_pdgId",
    }
    if not required <= set(arrays):
        return {}
    pt = arrays["ScoutingParticle_pt"]
    eta = arrays["ScoutingParticle_eta"]
    dxy = arrays["ScoutingParticle_dxy"]
    pdg_id = arrays["ScoutingParticle_pdgId"]
    sum_pt = awkward.sum(pt, axis=1)
    energy = abs(pt * np.cosh(eta))
    sum_energy = awkward.sum(energy, axis=1)
    pt_fraction = awkward.where(sum_pt > 0, pt / sum_pt, np.nan)
    energy_fraction = awkward.where(
        sum_energy > 0, energy / sum_energy, np.nan)
    abs_pid = abs(pdg_id)
    positive_encoding = (
        (pdg_id == 11) | (pdg_id == 13) | (pdg_id == 211))
    negative_encoding = (
        (pdg_id == -11) | (pdg_id == -13) | (pdg_id == -211))
    charge_encoding = awkward.where(
        positive_encoding, 1.0, awkward.where(negative_encoding, -1.0, 0.0))
    result = {
        "ScoutingParticle_model_log_pt_fraction": to_numpy(
            awkward.flatten(np.log(pt_fraction), axis=None), awkward),
        "ScoutingParticle_model_tanh_dxy": to_numpy(
            awkward.flatten(np.tanh(dxy), axis=None), awkward),
        "ScoutingParticle_model_log_energy_fraction": to_numpy(
            awkward.flatten(np.log(energy_fraction), axis=None), awkward),
        "ScoutingParticle_model_charge_encoding": to_numpy(
            awkward.flatten(charge_encoding, axis=None), awkward),
        "ScoutingParticle_model_abs_pdgId": to_numpy(
            awkward.flatten(abs_pid, axis=None), awkward),
        "ScoutingParticle_model_is_pf_direct_rows": np.ones(
            int(awkward.sum(awkward.num(pt, axis=1))), dtype=np.float64),
    }
    for target in MODEL_PDG_IDS:
        result[f"ScoutingParticle_model_pid_is_{target}"] = to_numpy(
            awkward.flatten(abs_pid == target, axis=None), awkward).astype(np.float64)
    return result


def score_residuals(event_columns: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    result = {}
    lowered = {name.lower(): name for name in event_columns}
    for name in list(event_columns):
        low = name.lower()
        if "recomputed" in low or not low.endswith("_score"):
            continue
        prefix = low[:-len("_score")]
        candidates = (
            low + "recomputed",
            prefix + "_recomputedscore",
            prefix + "_recomputed_score",
        )
        match = next((lowered[item] for item in candidates if item in lowered), None)
        if match:
            result[f"{name}__minus__{match}"] = event_columns[name] - event_columns[match]
    return result


def update_group_stats(
    accumulator: dict[tuple[int, int], dict[str, float | int]],
    process: np.ndarray,
    category: np.ndarray,
    gen_weight: np.ndarray,
) -> None:
    for process_id, category_id in np.unique(
            np.column_stack((process, category)), axis=0):
        selected = (process == process_id) & (category == category_id)
        weights = gen_weight[selected]
        target = accumulator[(int(process_id), int(category_id))]
        target["events"] += int(selected.sum())
        target["sum_gen_weight"] += float(weights.sum())
        target["sum_abs_gen_weight"] += float(np.abs(weights).sum())
        target["sum_gen_weight2"] += float(np.square(weights).sum())
        target["negative_gen_weights"] += int(np.count_nonzero(weights < 0))
        target["zero_gen_weights"] += int(np.count_nonzero(weights == 0))


def finish_group_stats(
    accumulator: dict[tuple[int, int], dict[str, float | int]],
    processes: dict[int, str],
    categories: dict[int, dict],
) -> list[dict]:
    rows = []
    for (process_id, category_id), entry in sorted(accumulator.items()):
        row = dict(entry)
        row.update({
            "process_id": process_id,
            "process": processes.get(process_id, f"unknown_{process_id}"),
            "category_id": category_id,
            "category": categories.get(category_id, {}).get(
                "name", f"unknown_{category_id}"),
        })
        row["effective_events_signed"] = (
            float(entry["sum_gen_weight"]) ** 2
            / float(entry["sum_gen_weight2"])
            if float(entry["sum_gen_weight2"]) > 0 else 0.0
        )
        row["effective_fraction"] = (
            float(row["effective_events_signed"]) / int(row["events"])
            if int(row["events"]) else 0.0)
        rows.append(row)
    return rows


def feature_payload(
    feature: str,
    event_reservoir: AlignedReservoir,
    candidate_reservoirs: dict[str, PriorityReservoir],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, str]:
    if feature in event_reservoir.columns:
        return (
            event_reservoir.columns[feature],
            event_reservoir.process,
            event_reservoir.category,
            np.abs(event_reservoir.gen_weight),
            "event",
        )
    reservoir = candidate_reservoirs[feature]
    return (
        reservoir.values,
        reservoir.process,
        reservoir.category,
        np.abs(reservoir.event_balanced_weight),
        "candidate",
    )


def valid_group(
    selected: np.ndarray,
    values: np.ndarray,
    minimum: int,
) -> np.ndarray:
    selected = selected & np.isfinite(values)
    return selected if int(selected.sum()) >= minimum else np.zeros_like(selected)


def plot_hist_lines(
    axis,
    values: np.ndarray,
    bins: np.ndarray,
    groups: list[tuple[str, np.ndarray, str]],
    weights: np.ndarray | None,
    minimum: int,
) -> int:
    drawn = 0
    centers = 0.5 * (bins[:-1] + bins[1:])
    for label, selected, color in groups:
        selected = valid_group(selected, values, minimum)
        if not selected.any():
            continue
        selected_weights = None if weights is None else weights[selected]
        hist, _ = normalized_histogram(values[selected], bins, selected_weights)
        if np.any(hist > 0):
            axis.step(centers, hist, where="mid", label=label, color=color, lw=1.45)
            drawn += 1
    return drawn


def save_figure(figure, path: Path, save_pdf: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=150, bbox_inches="tight")
    if save_pdf:
        figure.savefig(path.with_suffix(".pdf"), bbox_inches="tight")


def process_groups(
    process: np.ndarray,
    processes: dict[int, str],
) -> list[tuple[str, np.ndarray, str]]:
    return [
        (processes.get(process_id, str(process_id)), process == process_id,
         PROCESS_COLORS.get(process_id, f"C{process_id % 10}"))
        for process_id in sorted(processes)
    ]


def category_groups(
    process_id: int,
    category: np.ndarray,
    categories: dict[int, dict],
    pyplot,
) -> list[tuple[str, np.ndarray, str]]:
    ids = [
        category_id for category_id, entry in sorted(categories.items())
        if int(entry["process"]) == process_id and np.any(category == category_id)
    ]
    colors = pyplot.cm.turbo(np.linspace(0.02, 0.98, max(len(ids), 1)))
    return [
        (categories[category_id]["name"], category == category_id, colors[index])
        for index, category_id in enumerate(ids)
    ]


def decorate_distribution_axis(axis, feature: str, log_y: bool = False) -> None:
    axis.set_xlabel(feature.replace("__", ": ").replace("_", " "))
    axis.set_ylabel("Fraction / bin")
    axis.grid(alpha=0.2)
    if log_y:
        axis.set_yscale("log")
        axis.set_ylim(bottom=1e-6)


def plot_parent_comparison(
    feature: str,
    values: np.ndarray,
    process: np.ndarray,
    weights: np.ndarray,
    scope: str,
    bins: np.ndarray,
    processes: dict[int, str],
    output: Path,
    pyplot,
    minimum: int,
    save_pdf: bool,
) -> None:
    groups = process_groups(process, processes)
    figure, axes = pyplot.subplots(2, 2, figsize=(13, 9), sharex=True)
    for column, (title, selected_weights) in enumerate((
        ("Raw generated rows", None),
        ("Within-parent |genWeight| shape", weights),
    )):
        for row, log_y in enumerate((False, True)):
            plot_hist_lines(
                axes[row, column], values, bins, groups,
                selected_weights, minimum)
            axes[row, column].set_title(title + (" (log y)" if log_y else ""))
            decorate_distribution_axis(axes[row, column], feature, log_y)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        figure.legend(
            handles, labels, loc="upper center", ncol=4,
            bbox_to_anchor=(0.5, 0.935))
    weight_note = (
        "Candidate weighted view gives each event total |genWeight|; "
        "raw view counts candidates."
        if scope == "candidate" else
        "Each parent curve is normalized independently; this is a shape comparison."
    )
    figure.suptitle(f"{feature}\n{weight_note}", y=0.995, fontsize=11)
    figure.tight_layout(rect=(0, 0, 1, 0.88))
    save_figure(figure, output, save_pdf)
    pyplot.close(figure)


def plot_reweighting_by_process(
    feature: str,
    values: np.ndarray,
    process: np.ndarray,
    weights: np.ndarray,
    scope: str,
    bins: np.ndarray,
    processes: dict[int, str],
    output: Path,
    pyplot,
    minimum: int,
    save_pdf: bool,
) -> None:
    figure, axes = pyplot.subplots(2, 2, figsize=(13, 9), sharex=True)
    for axis, process_id in zip(axes.flat, sorted(processes)):
        selected = valid_group(process == process_id, values, minimum)
        if not selected.any():
            axis.set_visible(False)
            continue
        groups = [("raw", selected, "#4c78a8")]
        plot_hist_lines(axis, values, bins, groups, None, minimum)
        groups = [("|genWeight|", selected, "#e45756")]
        plot_hist_lines(axis, values, bins, groups, weights, minimum)
        axis.set_title(processes[process_id])
        decorate_distribution_axis(axis, feature, log_y=True)
        axis.legend(fontsize=8)
    note = "event-balanced candidate weighting" if scope == "candidate" else "event weighting"
    figure.suptitle(f"Raw versus generator-reweighted: {feature} ({note})")
    figure.tight_layout()
    save_figure(figure, output, save_pdf)
    pyplot.close(figure)


def plot_subcategories(
    feature: str,
    values: np.ndarray,
    process: np.ndarray,
    category: np.ndarray,
    weights: np.ndarray,
    scope: str,
    bins: np.ndarray,
    processes: dict[int, str],
    categories: dict[int, dict],
    output_root: Path,
    pyplot,
    minimum: int,
    save_pdf: bool,
) -> list[Path]:
    outputs = []
    for process_id in sorted(processes):
        selected_process = process == process_id
        groups = category_groups(process_id, category, categories, pyplot)
        groups = [
            (name, mask & selected_process, color) for name, mask, color in groups
        ]
        groups = [
            item for item in groups
            if valid_group(item[1], values, minimum).any()
        ]
        if len(groups) < 2:
            continue
        figure, axes = pyplot.subplots(1, 2, figsize=(16, 6), sharex=True)
        plot_hist_lines(axes[0], values, bins, groups, None, minimum)
        plot_hist_lines(axes[1], values, bins, groups, weights, minimum)
        axes[0].set_title("Raw generated rows")
        axes[1].set_title("Within-parent |genWeight| shape")
        for axis in axes:
            decorate_distribution_axis(axis, feature, log_y=True)
            axis.legend(fontsize=7, ncol=1)
        note = (
            "raw counts candidates; weighted gives each event total |genWeight|"
            if scope == "candidate" else "every category curve has unit area")
        figure.suptitle(f"{processes[process_id]} subcategories: {feature}\n{note}")
        figure.tight_layout()
        path = (
            output_root / safe_name(processes[process_id])
            / f"{safe_name(feature)}.png")
        save_figure(figure, path, save_pdf)
        pyplot.close(figure)
        outputs.append(path)
    return outputs


def pairwise_max_ks(
    values: np.ndarray,
    labels: np.ndarray,
    weights: np.ndarray | None,
    minimum: int,
) -> tuple[float, str]:
    ids = [int(value) for value in np.unique(labels)]
    best = (-1.0, "")
    for first, second in itertools.combinations(ids, 2):
        first_mask = valid_group(labels == first, values, minimum)
        second_mask = valid_group(labels == second, values, minimum)
        if not first_mask.any() or not second_mask.any():
            continue
        distance = weighted_ks(
            values[first_mask], values[second_mask],
            None if weights is None else weights[first_mask],
            None if weights is None else weights[second_mask],
        )
        if math.isfinite(distance) and distance > best[0]:
            best = (float(distance), f"{first} vs {second}")
    return (float("nan"), "") if best[0] < 0 else best


def category_vs_rest_max_ks(
    values: np.ndarray,
    process: np.ndarray,
    category: np.ndarray,
    weights: np.ndarray | None,
    categories: dict[int, dict],
    minimum: int,
) -> tuple[float, str]:
    best = (-1.0, "")
    for category_id, entry in categories.items():
        parent = int(entry["process"])
        selected = valid_group(category == category_id, values, minimum)
        rest = valid_group((process == parent) & (category != category_id), values, minimum)
        if not selected.any() or not rest.any():
            continue
        distance = weighted_ks(
            values[selected], values[rest],
            None if weights is None else weights[selected],
            None if weights is None else weights[rest],
        )
        if math.isfinite(distance) and distance > best[0]:
            best = (float(distance), entry["name"])
    return (float("nan"), "") if best[0] < 0 else best


def feature_ranking_row(
    feature: str,
    values: np.ndarray,
    process: np.ndarray,
    category: np.ndarray,
    weights: np.ndarray,
    scope: str,
    categories: dict[int, dict],
    minimum: int,
) -> dict[str, object]:
    parent_raw, parent_raw_pair = pairwise_max_ks(
        values, process, None, minimum)
    parent_weighted, parent_weighted_pair = pairwise_max_ks(
        values, process, weights, minimum)
    sub_raw, sub_raw_category = category_vs_rest_max_ks(
        values, process, category, None, categories, minimum)
    sub_weighted, sub_weighted_category = category_vs_rest_max_ks(
        values, process, category, weights, categories, minimum)
    reweighting = []
    for process_id in np.unique(process):
        selected = valid_group(process == process_id, values, minimum)
        if selected.any():
            reweighting.append(weighted_ks(
                values[selected], values[selected],
                np.ones(int(selected.sum())), weights[selected]))
    return {
        "feature": feature,
        "scope": scope,
        "sampled_entries": int(values.size),
        "parent_raw_max_ks": parent_raw,
        "parent_raw_pair": parent_raw_pair,
        "parent_genweighted_max_ks": parent_weighted,
        "parent_genweighted_pair": parent_weighted_pair,
        "subcategory_raw_max_ks": sub_raw,
        "subcategory_raw_driver": sub_raw_category,
        "subcategory_genweighted_max_ks": sub_weighted,
        "subcategory_genweighted_driver": sub_weighted_category,
        "raw_vs_genweighted_max_ks": (
            float(np.nanmax(reweighting)) if reweighting else float("nan")),
    }


def subcategory_separation_rows(
    feature: str,
    values: np.ndarray,
    process: np.ndarray,
    category: np.ndarray,
    weights: np.ndarray,
    scope: str,
    categories: dict[int, dict],
    minimum: int,
) -> list[dict]:
    rows = []
    for category_id, entry in categories.items():
        parent = int(entry["process"])
        selected = valid_group(category == category_id, values, minimum)
        rest = valid_group(
            (process == parent) & (category != category_id), values, minimum)
        if not selected.any() or not rest.any():
            continue
        raw_distance = weighted_ks(values[selected], values[rest])
        weighted_distance = weighted_ks(
            values[selected], values[rest], weights[selected], weights[rest])
        rows.append({
            "process_id": parent,
            "category_id": category_id,
            "category": entry["name"],
            "feature": feature,
            "scope": scope,
            "model_source_feature": source_collection(feature) != "event",
            "sampled_category_entries": int(selected.sum()),
            "sampled_rest_entries": int(rest.sum()),
            "raw_ks_vs_parent_rest": raw_distance,
            "genweighted_ks_vs_parent_rest": weighted_distance,
            "weighted_minus_raw_ks": weighted_distance - raw_distance,
        })
    return rows


def summarize_subcategory_separation(
    separation_rows: list[dict],
    group_rows: list[dict],
    processes: dict[int, str],
) -> list[dict]:
    exact = {
        (int(row["process_id"]), int(row["category_id"])): row
        for row in group_rows
    }
    grouped: dict[tuple[int, int, str], list[dict]] = defaultdict(list)
    for row in separation_rows:
        if row["model_source_feature"]:
            grouped[(
                int(row["process_id"]), int(row["category_id"]),
                str(row["category"]))].append(row)
    summaries = []
    for (process_id, category_id, category_name), rows in sorted(grouped.items()):
        distances = np.array([
            row["genweighted_ks_vs_parent_rest"] for row in rows], dtype=float)
        raw_distances = np.array([
            row["raw_ks_vs_parent_rest"] for row in rows], dtype=float)
        order = np.argsort(distances)[::-1]
        strongest = rows[int(order[0])]
        stats = exact.get((process_id, category_id), {})
        events = int(stats.get("events", 0))
        effective = float(stats.get("effective_events_signed", 0.0))
        persistent = int(np.count_nonzero(distances >= 0.20))
        moderate = int(np.count_nonzero(distances >= 0.10))
        raw_only = int(np.count_nonzero(
            (raw_distances >= 0.20) & (distances < 0.10)))
        threshold = max(3, int(math.ceil(0.10 * distances.size)))
        if effective < 100 or events < 100:
            guidance = "statistics too small for a taxonomy decision"
        elif persistent >= threshold:
            guidance = (
                "investigate as a persistent internal domain; test held-out confusion "
                "before splitting or excluding")
        elif raw_only >= threshold:
            guidance = "mostly a raw-sampling mode; generator weighting reduces it"
        else:
            guidance = "no broad weighted separation across audited model-source features"
        summaries.append({
            "process_id": process_id,
            "process": processes.get(process_id, str(process_id)),
            "category_id": category_id,
            "category": category_name,
            "events": events,
            "effective_events_signed": effective,
            "audited_model_source_features": int(distances.size),
            "features_weighted_ks_ge_0p10": moderate,
            "features_weighted_ks_ge_0p20": persistent,
            "features_raw_ge_0p20_but_weighted_lt_0p10": raw_only,
            "median_raw_ks": float(np.median(raw_distances)),
            "median_genweighted_ks": float(np.median(distances)),
            "maximum_genweighted_ks": float(distances[order[0]]),
            "strongest_feature": strongest["feature"],
            "guidance": guidance,
        })
    return summaries


def plot_overview(
    rows: list[dict],
    processes: dict[int, str],
    output_dir: Path,
    pyplot,
    save_pdf: bool,
) -> list[Path]:
    outputs = []
    labels = [row["category"] for row in rows]
    counts = np.array([row["events"] for row in rows], dtype=float)
    mean_abs = np.array([
        row["sum_abs_gen_weight"] / max(row["events"], 1) for row in rows])
    ess_fraction = np.array([row["effective_fraction"] for row in rows])
    colors = [PROCESS_COLORS.get(int(row["process_id"]), "gray") for row in rows]
    x = np.arange(len(rows))

    figure, axes = pyplot.subplots(3, 1, figsize=(18, 14), sharex=True)
    axes[0].bar(x, counts, color=colors)
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Raw events")
    axes[1].bar(x, mean_abs, color=colors)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Mean |genWeight|")
    axes[2].bar(x, ess_fraction, color=colors)
    axes[2].set_ylabel("Signed ESS / N")
    axes[2].set_ylim(0, max(1.0, float(np.nanmax(ess_fraction)) * 1.05))
    axes[2].set_xticks(x, labels, rotation=75, ha="right", fontsize=8)
    for axis in axes:
        axis.grid(axis="y", alpha=0.2)
    figure.suptitle(
        "Background production composition and generator-weight diagnostics\n"
        "Mean weights are not comparable cross sections across parent samples")
    figure.tight_layout()
    path = output_dir / "category_counts_and_weights.png"
    save_figure(figure, path, save_pdf)
    pyplot.close(figure)
    outputs.append(path)

    process_rows = []
    for process_id in sorted(processes):
        selected = [row for row in rows if int(row["process_id"]) == process_id]
        process_rows.append((
            processes[process_id], sum(int(row["events"]) for row in selected)))
    figure, axis = pyplot.subplots(figsize=(8, 5))
    axis.bar(
        [item[0] for item in process_rows],
        [item[1] for item in process_rows],
        color=[PROCESS_COLORS.get(pid, "gray") for pid in sorted(processes)])
    axis.set_ylabel("Raw generated events")
    axis.set_title("Parent-background row counts (not physical yields)")
    axis.grid(axis="y", alpha=0.2)
    figure.tight_layout()
    path = output_dir / "parent_event_counts.png"
    save_figure(figure, path, save_pdf)
    pyplot.close(figure)
    outputs.append(path)
    return outputs


def parent_weight_summary(
    rows: list[dict], processes: dict[int, str],
) -> list[dict]:
    result = []
    for process_id in sorted(processes):
        selected = [row for row in rows if int(row["process_id"]) == process_id]
        events = sum(int(row["events"]) for row in selected)
        sum_weight = sum(float(row["sum_gen_weight"]) for row in selected)
        sum_abs = sum(float(row["sum_abs_gen_weight"]) for row in selected)
        sum_weight2 = sum(float(row["sum_gen_weight2"]) for row in selected)
        mean_abs = sum_abs / events if events else 0.0
        result.append({
            "process_id": process_id,
            "process": processes[process_id],
            "events": events,
            "sum_gen_weight": sum_weight,
            "sum_abs_gen_weight": sum_abs,
            "sum_gen_weight2": sum_weight2,
            "mean_abs_gen_weight": mean_abs,
            "recommended_parent_divisor": mean_abs if mean_abs > 0 else 1.0,
            "sum_abs_after_parent_normalization": (
                sum_abs / mean_abs if mean_abs > 0 else 0.0),
            "negative_gen_weights": sum(
                int(row["negative_gen_weights"]) for row in selected),
            "effective_events_signed": (
                sum_weight**2 / sum_weight2 if sum_weight2 > 0 else 0.0),
        })
    mean_events = np.mean([row["events"] for row in result]) if result else 0.0
    for row in result:
        row["illustrative_equal_parent_factor_after_parent_normalization"] = (
            float(mean_events) / int(row["events"]) if int(row["events"]) else 0.0)
    return result


def plot_weight_scenarios(
    rows: list[dict],
    output: Path,
    pyplot,
    save_pdf: bool,
) -> Path:
    labels = [row["process"] for row in rows]
    x = np.arange(len(rows))
    colors = [PROCESS_COLORS.get(int(row["process_id"]), "gray") for row in rows]
    raw_events = np.array([row["events"] for row in rows], dtype=float)
    raw_gen_mass = np.array([row["sum_abs_gen_weight"] for row in rows], dtype=float)
    normalized_gen_mass = np.array([
        row["sum_abs_after_parent_normalization"] for row in rows], dtype=float)
    equal_parent = np.ones(len(rows), dtype=float)
    panels = (
        (raw_events, "Raw row counts", False),
        (raw_gen_mass, "Raw |genWeight| mass (arbitrary cross-parent units)", True),
        (normalized_gen_mass, "Parent-normalized generator-only mass (AE input)", False),
        (equal_parent, "Equal-parent target (illustrative NURD class balance)", False),
    )
    figure, axes = pyplot.subplots(2, 2, figsize=(12, 9))
    for axis, (values, title, log_y) in zip(axes.flat, panels):
        axis.bar(x, values, color=colors)
        axis.set_xticks(x, labels, rotation=20)
        axis.set_title(title)
        axis.set_ylabel("Relative mass")
        if log_y and np.all(values > 0):
            axis.set_yscale("log")
        axis.grid(axis="y", alpha=0.2)
    figure.suptitle(
        "Weighting stages: generator-shape correction is separate from class balance")
    figure.tight_layout()
    save_figure(figure, output, save_pdf)
    pyplot.close(figure)
    return output


def plot_correlations(
    event_reservoir: AlignedReservoir,
    ranking: list[dict],
    processes: dict[int, str],
    output_dir: Path,
    pyplot,
    save_pdf: bool,
    maximum_features: int = 24,
) -> list[Path]:
    ranked_names = [
        row["feature"] for row in sorted(
            (row for row in ranking if row["scope"] == "event"),
            key=lambda row: float(row["parent_genweighted_max_ks"])
            if math.isfinite(float(row["parent_genweighted_max_ks"])) else -1,
            reverse=True)
    ]
    names = ranked_names[:maximum_features]
    if len(names) < 2:
        return []
    outputs = []
    for process_id in sorted(processes):
        selected = event_reservoir.process == process_id
        if int(selected.sum()) < 10:
            continue
        matrix = np.column_stack([
            event_reservoir.columns[name][selected] for name in names])
        finite_fraction = np.mean(np.isfinite(matrix), axis=0)
        keep = finite_fraction > 0.98
        local_names = [name for name, selected_name in zip(names, keep) if selected_name]
        matrix = matrix[:, keep]
        if matrix.shape[1] < 2:
            continue
        weights = np.abs(event_reservoir.gen_weight[selected])
        correlation = weighted_correlation(matrix, weights)
        figure, axis = pyplot.subplots(figsize=(14, 12))
        image = axis.imshow(correlation, vmin=-1, vmax=1, cmap="coolwarm")
        axis.set_xticks(
            np.arange(len(local_names)), local_names, rotation=90, fontsize=7)
        axis.set_yticks(np.arange(len(local_names)), local_names, fontsize=7)
        axis.set_title(
            f"{processes[process_id]} generator-weighted event-feature correlation")
        figure.colorbar(image, ax=axis, label="Weighted Pearson correlation")
        figure.tight_layout()
        path = output_dir / f"{safe_name(processes[process_id])}.png"
        save_figure(figure, path, save_pdf)
        pyplot.close(figure)
        outputs.append(path)
    return outputs


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def contract_payload(present: set[str]) -> dict:
    return {
        "status": "audit contract; converter implementation still required",
        "classifier_pf_tensor": {
            "shape": "[events, configurable padded candidates, 7]",
            "columns": ["pt", "eta", "phi", "dxy", "dxysig", "is_pf", "pdgId"],
            "direct_source": "ScoutingParticle",
            "direct_source_fields": ["pt", "eta", "phi", "dxy", "dxysig", "pdgId"],
            "legacy_observation": (
                "The old tensor also has is_pf=0 rows associated with high-level "
                "muon/electron/photon objects. Their exact merge and deduplication "
                "must be made explicit in the converter, not guessed by this audit."),
            "current_model_behavior": (
                "The model infers the padded candidate count from each converted tensor; "
                "valid rows require pt>0 and padding is pt=0."),
            "audited_transforms": [
                "log(pt / event sum pt)", "eta", "phi", "tanh(dxy)",
                "dxysig", "log(abs(pt*cosh(eta)) / event sum energy)",
                "model charge encoding", "direct-row is_pf", "absolute-pdgId one-hot"],
        },
        "autoencoder_object_tensor": {
            "shape": "[events, 23, 4] in the legacy sample",
            "columns": ["pt", "eta", "phi", "object_type"],
            "legacy_slots": {
                "0:10": "jets, type 0",
                "10:14": "muons, type 1",
                "14:18": "electrons, type 2",
                "18:22": "photons, type 3",
                "22:23": "MET, type 4",
            },
            "source_tables_to_audit": [
                "ScoutingPFJet", "ScoutingMuonVtx", "ScoutingMuonNoVtx",
                "ScoutingElectron", "ScoutingPhoton", "ScoutingMET"],
            "open_converter_decision": (
                "Choose and document the muon table, object selections, ordering, "
                "padding, truncation, and overlap treatment before conversion."),
            "audited_proxy": (
                "The audit plots per-collection candidate fields and pT-ranked slot0..N "
                "pt/eta/phi values using the legacy 10/4/4/4 slot capacities."),
        },
        "weight_tensor": {
            "source": "per-event NanoAOD genWeight for all four parents",
            "required_change": (
                "utils/hlt_weights.load_generator_weights currently applies an external "
                "weight vector only to one configured label and forces other labels to 1. "
                "It must support aligned all-event weights before this production is trained."),
            "parent_scale_policy": (
                "Normalize mean abs(genWeight) to one separately in each parent before "
                "combining parents. Raw generator scales are not cross-process yields."),
            "ae": (
                "parent-normalized generator-only event weighting; no class "
                "equalization in the loss"),
            "nurd": (
                "generator weights are subsequently combined with class and continuous-"
                "nuisance balancing for classifier/NURD losses"),
        },
        "present_source_tables": sorted({
            collection_from_branch(branch) for branch in present
            if collection_from_branch(branch) is not None
        }),
    }


def write_html_index(
    output_dir: Path,
    ranking: list[dict],
    subcategory_summary: list[dict],
    plot_paths: list[Path],
    process_count: int,
    event_count: int,
) -> None:
    relative_plots = sorted({str(path.relative_to(output_dir)) for path in plot_paths})
    top = sorted(
        ranking,
        key=lambda row: float(row["subcategory_genweighted_max_ks"])
        if math.isfinite(float(row["subcategory_genweighted_max_ks"])) else -1,
        reverse=True,
    )[:30]
    rows = "\n".join(
        "<tr>" + "".join(
            f"<td>{html.escape(str(row[key]))}</td>"
            for key in (
                "feature", "scope", "subcategory_genweighted_max_ks",
                "subcategory_genweighted_driver", "raw_vs_genweighted_max_ks"))
        + "</tr>" for row in top)
    category_rows = "\n".join(
        "<tr>" + "".join(
            f"<td>{html.escape(str(row[key]))}</td>"
            for key in (
                "process", "category", "events", "effective_events_signed",
                "features_weighted_ks_ge_0p20", "strongest_feature", "guidance"))
        + "</tr>" for row in subcategory_summary)
    cards = "\n".join(
        f'<a class="card" href="{html.escape(path)}">'
        f'<img loading="lazy" src="{html.escape(path)}"><span>{html.escape(path)}</span></a>'
        for path in relative_plots if path.endswith(".png"))
    document = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Background NanoAOD audit</title>
<style>
body{{font-family:sans-serif;margin:2rem;color:#222}} .note{{max-width:90rem}}
table{{border-collapse:collapse}} th,td{{border:1px solid #ccc;padding:.35rem}}
.grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(360px,1fr));gap:1rem}}
.card{{color:inherit;text-decoration:none;border:1px solid #ddd;padding:.5rem}}
.card img{{width:100%;height:260px;object-fit:contain}} .card span{{font-size:.78rem}}
</style></head><body>
<h1>Background NanoAOD audit</h1>
<p>{event_count:,} events across {process_count} parent backgrounds.</p>
<p class="note"><b>Interpretation:</b> raw curves show generated row sampling. Weighted
curves use within-parent |genWeight| shapes. Every curve has unit area; cross-process
physical yields require cross sections. Candidate weighted curves are event-balanced so
high-multiplicity events do not receive extra weight merely for containing more objects.</p>
<h2>Subcategory decision aid</h2>
<p class="note">This is a review queue, not an automatic split/exclusion decision. It
uses only model-source features and omits truth category variables.</p>
<table><tr><th>Parent</th><th>Category</th><th>events</th><th>signed ESS</th>
<th>features with weighted KS ≥ 0.20</th><th>strongest feature</th><th>guidance</th></tr>
{category_rows}</table>
<h2>Largest weighted within-parent subcategory differences</h2>
<table><tr><th>Feature</th><th>scope</th><th>KS</th><th>driver</th><th>raw/weighted KS</th></tr>
{rows}</table><h2>Plots</h2><div class="grid">{cards}</div></body></html>"""
    (output_dir / "index.html").write_text(document)


def main() -> int:
    args = parse_args()
    inputs = expand_inputs(args.inputs)
    if args.max_files > 0:
        inputs = inputs[:args.max_files]
    if not inputs:
        raise SystemExit("No input ROOT files were supplied.")

    try:
        import awkward
        import uproot
    except ImportError as error:
        raise SystemExit(
            "This audit needs uproot and awkward in the active environment.") from error
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as pyplot

    processes, categories = load_schema(args.categories)
    expected_categories = parse_category_spec(
        args.expected_categories, set(categories))
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    save_pdf = args.pdf

    first_tree = uproot.open(f"{inputs[0]}:Events")
    branches = discover_branches(first_tree, args.profile)
    probe = first_tree.arrays(branches, entry_stop=2, library="ak", how=dict)
    branch_ndim = {name: awkward_ndim(probe[name]) for name in branches}
    scalar_branches = [
        name for name in branches
        if branch_ndim[name] == 1
        and name not in LABEL_BRANCHES
        and name not in IDENTIFIER_BRANCHES
    ]
    jagged_branches = [name for name in branches if branch_ndim[name] > 1]
    derivations = derived_feature_names(jagged_branches)
    pf_derived_names = sorted(pf_model_candidate_features(probe, awkward))

    print(
        f"Discovered {len(branches)} selected branches: "
        f"{len(scalar_branches)} scalar QA, {len(jagged_branches)} jagged QA, "
        f"and {len(pf_derived_names)} model-derived PF features",
        flush=True)

    rng = np.random.default_rng(args.seed)
    event_reservoir = AlignedReservoir(args.event_reservoir, rng)
    candidate_reservoirs = {
        name: PriorityReservoir(args.candidate_reservoir, rng)
        for name in jagged_branches + pf_derived_names
    }
    feature_stats: dict[str, FeatureStats] = defaultdict(FeatureStats)
    group_stats = defaultdict(lambda: {
        "events": 0,
        "sum_gen_weight": 0.0,
        "sum_abs_gen_weight": 0.0,
        "sum_gen_weight2": 0.0,
        "negative_gen_weights": 0,
        "zero_gen_weights": 0,
    })
    processed_files: list[str] = []
    bad_files: list[dict[str, str]] = []
    event_count = 0
    maximum_score_residuals: dict[str, float] = defaultdict(float)

    for file_index, input_path in enumerate(inputs, start=1):
        try:
            with uproot.open(input_path) as source:
                tree = source["Events"]
                missing = sorted(set(branches) - set(tree.keys()))
                if missing:
                    raise RuntimeError(f"file is missing selected branches: {missing}")
                for arrays in tree.iterate(
                        branches, step_size=args.step_size, library="ak", how=dict):
                    process = to_numpy(arrays["Background_process"], awkward).astype(np.int8)
                    category = to_numpy(arrays["Background_category"], awkward).astype(np.int16)
                    gen_weight = to_numpy(arrays["genWeight"], awkward).astype(np.float64)
                    n_events = process.size
                    if category.size != n_events or gen_weight.size != n_events:
                        raise RuntimeError("event metadata branches do not align")
                    event_count += n_events
                    update_group_stats(group_stats, process, category, gen_weight)

                    event_columns = {
                        name: to_numpy(arrays[name], awkward).astype(np.float64)
                        for name in scalar_branches
                    }
                    event_columns["genWeight__log10_abs"] = np.log10(
                        np.where(np.abs(gen_weight) > 0, np.abs(gen_weight), np.nan))
                    event_columns["genWeight__sign"] = np.sign(gen_weight)
                    event_columns.update(
                        derive_event_features(arrays, derivations, awkward))
                    residuals = score_residuals(event_columns)
                    event_columns.update(residuals)
                    for name, values in residuals.items():
                        finite = np.abs(values[np.isfinite(values)])
                        if finite.size:
                            maximum_score_residuals[name] = max(
                                maximum_score_residuals[name], float(finite.max()))
                    for name, values in event_columns.items():
                        feature_stats[name].update(values)
                    event_reservoir.add(
                        event_columns, process, category, gen_weight)

                    for name in jagged_branches:
                        array = arrays[name]
                        flat, counts = flatten_per_event(array, awkward)
                        feature_stats[name].update(flat)
                        if not flat.size:
                            continue
                        repeated_process = np.repeat(process, counts)
                        repeated_category = np.repeat(category, counts)
                        repeated_weight = np.repeat(gen_weight, counts)
                        inverse_count = np.repeat(
                            np.divide(
                                1.0, counts,
                                out=np.zeros_like(counts, dtype=np.float64),
                                where=counts > 0),
                            counts)
                        candidate_reservoirs[name].add(
                            flat, repeated_process, repeated_category,
                            repeated_weight,
                            event_balanced_weight=repeated_weight * inverse_count)

                    derived_pf = pf_model_candidate_features(arrays, awkward)
                    if derived_pf:
                        counts = to_numpy(
                            awkward.num(
                                arrays["ScoutingParticle_pt"], axis=1), awkward
                        ).astype(np.int64)
                        repeated_process = np.repeat(process, counts)
                        repeated_category = np.repeat(category, counts)
                        repeated_weight = np.repeat(gen_weight, counts)
                        inverse_count = np.repeat(
                            np.divide(
                                1.0, counts,
                                out=np.zeros_like(counts, dtype=np.float64),
                                where=counts > 0),
                            counts)
                        for name, flat in derived_pf.items():
                            feature_stats[name].update(flat)
                            candidate_reservoirs[name].add(
                                flat, repeated_process, repeated_category,
                                repeated_weight,
                                event_balanced_weight=(
                                    repeated_weight * inverse_count))
            processed_files.append(input_path)
        except Exception as error:
            bad_files.append({"file": input_path, "error": str(error)})
            print(f"ERROR reading {input_path}: {error}", file=sys.stderr, flush=True)
            if not args.skip_bad_files:
                raise
        if file_index == 1 or file_index % 25 == 0 or file_index == len(inputs):
            print(
                f"[{file_index}/{len(inputs)}] files; {event_count:,} events; "
                f"bad={len(bad_files)}", flush=True)

    if not processed_files or not event_count:
        raise SystemExit("No events were successfully processed.")

    group_rows = finish_group_stats(group_stats, processes, categories)
    write_csv(output_dir / "category_weight_summary.csv", group_rows)
    parent_rows = parent_weight_summary(group_rows, processes)
    write_csv(output_dir / "parent_weight_summary.csv", parent_rows)
    plot_paths = plot_overview(
        group_rows, processes, output_dir / "overview", pyplot, save_pdf)
    plot_paths.append(plot_weight_scenarios(
        parent_rows, output_dir / "overview" / "weighting_scenarios.png",
        pyplot, save_pdf))

    all_features = sorted(
        list(event_reservoir.columns) + [
            name for name, reservoir in candidate_reservoirs.items()
            if reservoir.values.size
        ])
    ranking = []
    inventory = []
    subcategory_separation = []
    plotted_features = 0
    for feature_index, feature in enumerate(all_features, start=1):
        values, process, category, weights, scope = feature_payload(
            feature, event_reservoir, candidate_reservoirs)
        finite = np.isfinite(values)
        values, process, category, weights = (
            values[finite], process[finite], category[finite], weights[finite])
        if not values.size:
            continue
        bins, discrete = make_bins(values)
        row = feature_ranking_row(
            feature, values, process, category, weights, scope,
            categories, args.min_category_entries)
        ranking.append(row)
        subcategory_separation.extend(subcategory_separation_rows(
            feature, values, process, category, weights, scope,
            categories, args.min_category_entries))
        inventory.append({
            "feature": feature,
            "scope": scope,
            "source_collection": source_collection(feature),
            "discrete": discrete,
            "sampled_entries": int(values.size),
            **feature_stats[feature].as_dict(),
            "q001": float(np.quantile(values, 0.001)),
            "q01": float(np.quantile(values, 0.01)),
            "median": float(np.quantile(values, 0.5)),
            "q99": float(np.quantile(values, 0.99)),
            "q999": float(np.quantile(values, 0.999)),
        })

        if (args.max_plotted_features <= 0
                or plotted_features < args.max_plotted_features):
            parent_path = (
                output_dir / "features" / scope
                / f"{safe_name(feature)}__parents.png")
            plot_parent_comparison(
                feature, values, process, weights, scope, bins, processes,
                parent_path, pyplot, args.min_category_entries, save_pdf)
            plot_paths.append(parent_path)

            reweight_path = (
                output_dir / "reweighting" / scope
                / f"{safe_name(feature)}.png")
            plot_reweighting_by_process(
                feature, values, process, weights, scope, bins, processes,
                reweight_path, pyplot, args.min_category_entries, save_pdf)
            plot_paths.append(reweight_path)

            plot_paths.extend(plot_subcategories(
                feature, values, process, category, weights, scope, bins,
                processes, categories, output_dir / "subcategories" / scope,
                pyplot, args.min_category_entries, save_pdf))
            plotted_features += 1
        if feature_index == 1 or feature_index % 10 == 0 or feature_index == len(all_features):
            print(
                f"Audited {feature_index}/{len(all_features)} features; "
                f"plotted={plotted_features}",
                flush=True)

    ranking.sort(
        key=lambda row: float(row["subcategory_genweighted_max_ks"])
        if math.isfinite(float(row["subcategory_genweighted_max_ks"])) else -1,
        reverse=True)
    write_csv(output_dir / "feature_ranking.csv", ranking)
    write_csv(output_dir / "feature_inventory.csv", inventory)
    write_csv(
        output_dir / "subcategory_feature_separation.csv",
        subcategory_separation)
    subcategory_summary = summarize_subcategory_separation(
        subcategory_separation, group_rows, processes)
    write_csv(
        output_dir / "subcategory_decision_aid.csv", subcategory_summary)
    plot_paths.extend(plot_correlations(
        event_reservoir, ranking, processes, output_dir / "correlations",
        pyplot, save_pdf))

    observed_categories = {int(row["category_id"]) for row in group_rows}
    unknown_categories = sorted(observed_categories - set(categories))
    process_mismatches = [
        {
            "category_id": int(row["category_id"]),
            "observed_process": int(row["process_id"]),
            "expected_process": int(categories[int(row["category_id"])]["process"]),
        }
        for row in group_rows
        if int(row["category_id"]) in categories
        and int(row["process_id"]) != int(categories[int(row["category_id"])]["process"])
    ]
    score_families = {
        family for family in ("AXO", "CICADA")
        if any(family.lower() in name.lower() for name in maximum_score_residuals)
    }
    summary = {
        "ok": (
            not bad_files
            and not (expected_categories - observed_categories)
            and not unknown_categories
            and not process_mismatches
            and score_families == {"AXO", "CICADA"}
            and all(value <= 1e-6 for value in maximum_score_residuals.values())
        ),
        "input_files_requested": len(inputs),
        "input_files_processed": len(processed_files),
        "events": event_count,
        "selected_branches": branches,
        "scalar_features": len(event_reservoir.columns),
        "candidate_features": sum(
            reservoir.values.size > 0 for reservoir in candidate_reservoirs.values()),
        "plots_png": sum(path.suffix == ".png" for path in plot_paths),
        "features_plotted": plotted_features,
        "subcategory_decision_aid_rows": len(subcategory_summary),
        "profile": args.profile,
        "bad_files": bad_files,
        "missing_expected_categories": sorted(expected_categories - observed_categories),
        "known_but_not_expected_categories_absent": sorted(
            set(categories) - expected_categories - observed_categories),
        "unknown_categories": unknown_categories,
        "category_process_mismatches": process_mismatches,
        "maximum_score_residuals": dict(maximum_score_residuals),
        "score_closure_families_found": sorted(score_families),
        "weighting_interpretation": {
            "raw": "generated rows/candidates; no generator correction",
            "generator_weighted": (
                "absolute genWeight, normalized independently within every displayed "
                "curve; candidate views divide each event weight among its candidates"),
            "recommended_ae_input_weight": (
                "abs(genWeight) divided by mean abs(genWeight) separately inside each "
                "parent sample; this preserves within-parent ratios and raw parent row "
                "counts while removing arbitrary generator-unit differences"),
            "negative_weight_policy": (
                "the current training loader requires non-negative weights; if this report "
                "finds negative weights, an explicit signed-weight policy is required"),
            "cross_process_yields": (
                "not produced: parent-sample cross sections are required"),
            "class_balancing": (
                "not applied: this audit is intentionally upstream of NURD class and "
                "nuisance-stratum balancing"),
        },
        "decision_guidance": [
            "Treat raw-only modes that disappear after generator weighting as sampling artifacts.",
            "Consider a subcategory only when weighted differences recur across several model inputs and its effective statistics are adequate.",
            "Strong separation between parent labels is expected classifier information, not by itself a reason to remove a sample.",
            "Do not remove an overlapping or 'leaky' sample solely because it resembles another Standard Model parent; overlap can be physically valid.",
            "Use held-out classifier confusion and stability tests after conversion before changing the taxonomy.",
        ],
    }
    (output_dir / "audit_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n")
    (output_dir / "model_input_contract.json").write_text(
        json.dumps(contract_payload(set(branches)), indent=2, sort_keys=True) + "\n")
    write_html_index(
        output_dir, ranking, subcategory_summary, plot_paths,
        len(processes), event_count)

    print(
        f"Wrote {output_dir / 'index.html'} and {summary['plots_png']} PNG plots\n"
        f"Files: {len(processed_files)}/{len(inputs)}; events: {event_count}; "
        f"features: {len(ranking)}; ok={summary['ok']}",
        flush=True)
    return 0 if summary["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
