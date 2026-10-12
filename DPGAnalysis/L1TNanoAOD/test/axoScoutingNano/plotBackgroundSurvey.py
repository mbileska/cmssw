#!/usr/bin/env python3
"""Make unweighted QA plots and a JSON summary from background survey NanoAOD files."""

import argparse
import json
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot


BRANCHES = (
    "genWeight",
    "L1AXO_score",
    "L1AXO_scoreRecomputed",
    "L1CICADA_score",
    "L1CICADA_scoreRecomputed",
    "Background_process",
    "Background_sample",
    "Background_category",
    "Background_kinematicValue",
)

SAMPLE_FIELDS = (
    "L1AXO_score",
    "L1AXO_scoreRecomputed",
    "L1CICADA_score",
    "L1CICADA_scoreRecomputed",
    "Background_process",
    "Background_sample",
    "Background_category",
    "Background_kinematicValue",
    "genWeight",
)


class PrioritySample:
    """Uniform, mergeable reservoir based on random priorities."""

    def __init__(self, capacity, rng):
        self.capacity = capacity
        self.rng = rng
        self.priorities = np.empty(0, dtype=np.float64)
        self.data = {field: np.empty(0) for field in SAMPLE_FIELDS}

    def add(self, arrays, selection=None):
        if selection is None:
            selection = slice(None)
        incoming = {field: np.asarray(arrays[field])[selection] for field in SAMPLE_FIELDS}
        size = len(incoming[SAMPLE_FIELDS[0]])
        if size == 0:
            return
        priorities = self.rng.random(size)
        if len(self.priorities):
            priorities = np.concatenate((self.priorities, priorities))
            incoming = {
                field: np.concatenate((self.data[field], incoming[field]))
                for field in SAMPLE_FIELDS
            }
        if len(priorities) > self.capacity:
            keep = np.argpartition(priorities, -self.capacity)[-self.capacity :]
            priorities = priorities[keep]
            incoming = {field: values[keep] for field, values in incoming.items()}
        self.priorities = priorities
        self.data = incoming


def parse_args():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__, fromfile_prefix_chars="@")
    parser.add_argument(
        "inputs",
        nargs="+",
        help="NanoAOD ROOT files or xrootd URLs; use @filelist.txt for one path per line",
    )
    parser.add_argument(
        "--categories",
        default=str(here / "background_categories.json"),
        help="exclusive-category schema JSON",
    )
    parser.add_argument("--output-dir", default="survey_plots", help="plot/report directory")
    parser.add_argument("--step-size", default="100 MB", help="uproot iteration step size")
    parser.add_argument(
        "--max-global-points",
        type=int,
        default=200000,
        help="event-representative reservoir size for process plots",
    )
    parser.add_argument(
        "--max-points-per-category",
        type=int,
        default=20000,
        help="reservoir size for each category plot",
    )
    parser.add_argument("--score-tolerance", type=float, default=1.0e-6)
    parser.add_argument("--seed", type=int, default=12345)
    return parser.parse_args()


def new_stats():
    return {
        "events": 0,
        "sum_gen_weight": 0.0,
        "sum_gen_weight2": 0.0,
        "sum_abs_gen_weight": 0.0,
        "negative_gen_weights": 0,
        "zero_gen_weights": 0,
        "sum_axo_score": 0.0,
        "sum_cicada_score": 0.0,
        "max_abs_axo_difference": 0.0,
        "max_abs_cicada_difference": 0.0,
    }


def update_stats(stats, arrays, selection):
    weights = np.asarray(arrays["genWeight"])[selection].astype(np.float64, copy=False)
    axo = np.asarray(arrays["L1AXO_score"])[selection].astype(np.float64, copy=False)
    axo_recomputed = np.asarray(arrays["L1AXO_scoreRecomputed"])[selection].astype(
        np.float64, copy=False
    )
    cicada = np.asarray(arrays["L1CICADA_score"])[selection].astype(np.float64, copy=False)
    cicada_recomputed = np.asarray(arrays["L1CICADA_scoreRecomputed"])[selection].astype(
        np.float64, copy=False
    )
    stats["events"] += len(weights)
    stats["sum_gen_weight"] += float(np.sum(weights, dtype=np.float64))
    stats["sum_gen_weight2"] += float(np.sum(weights * weights, dtype=np.float64))
    stats["sum_abs_gen_weight"] += float(np.sum(np.abs(weights), dtype=np.float64))
    stats["negative_gen_weights"] += int(np.count_nonzero(weights < 0))
    stats["zero_gen_weights"] += int(np.count_nonzero(weights == 0))
    stats["sum_axo_score"] += float(np.sum(axo, dtype=np.float64))
    stats["sum_cicada_score"] += float(np.sum(cicada, dtype=np.float64))
    stats["max_abs_axo_difference"] = max(
        stats["max_abs_axo_difference"],
        float(np.max(np.abs(axo_recomputed - axo), initial=0.0)),
    )
    stats["max_abs_cicada_difference"] = max(
        stats["max_abs_cicada_difference"],
        float(np.max(np.abs(cicada_recomputed - cicada), initial=0.0)),
    )


def finite(values):
    values = np.asarray(values, dtype=np.float64)
    return values[np.isfinite(values)]


def robust_range(values, lower=0.002, upper=0.998):
    values = finite(values)
    if not len(values):
        return (0.0, 1.0)
    low, high = np.quantile(values, [lower, upper])
    if low == high:
        width = max(abs(float(low)) * 0.05, 1.0e-6)
        low -= width
        high += width
    return float(low), float(high)


def safe_name(name):
    return "".join(character if character.isalnum() else "_" for character in name).strip("_")


def add_sample_quantiles(entry, sample):
    for output_name, branch in (
        ("axo_score", "L1AXO_score"),
        ("cicada_score", "L1CICADA_score"),
        ("kinematic_value", "Background_kinematicValue"),
    ):
        values = finite(sample.data[branch])
        entry[f"sampled_{output_name}_quantiles"] = (
            {
                "q00": float(np.min(values)),
                "q25": float(np.quantile(values, 0.25)),
                "q50": float(np.quantile(values, 0.50)),
                "q75": float(np.quantile(values, 0.75)),
                "q99": float(np.quantile(values, 0.99)),
                "q100": float(np.max(values)),
            }
            if len(values)
            else None
        )


def plot_category_counts(output_dir, observed_ids, category_names, category_stats):
    counts = [category_stats[category_id]["events"] for category_id in observed_ids]
    labels = [f"{category_id}: {category_names[category_id]}" for category_id in observed_ids]
    width = max(10, 0.38 * len(observed_ids))
    figure, axis = plt.subplots(figsize=(width, 6))
    axis.bar(np.arange(len(labels)), counts, color="#4472C4")
    axis.set_yscale("log")
    axis.set_ylabel("Events")
    axis.set_title("Exclusive background-category populations")
    axis.set_xticks(np.arange(len(labels)), labels, rotation=65, ha="right", fontsize=8)
    axis.grid(axis="y", which="both", alpha=0.25)
    figure.tight_layout()
    figure.savefig(output_dir / "category_counts.png", dpi=180)
    plt.close(figure)


def plot_process_scores(output_dir, global_sample, process_names):
    process_ids = np.asarray(global_sample.data["Background_process"], dtype=np.int16)
    axo_all = finite(global_sample.data["L1AXO_score"])
    cicada_all = finite(global_sample.data["L1CICADA_score"])
    ranges = {
        "L1AXO_score": robust_range(axo_all),
        "L1CICADA_score": robust_range(cicada_all),
    }
    figure, axes = plt.subplots(2, 1, figsize=(9, 9))
    for axis, branch, title in (
        (axes[0], "L1AXO_score", "AXOL1TL score"),
        (axes[1], "L1CICADA_score", "CICADA Model v2.2.0 score"),
    ):
        values = np.asarray(global_sample.data[branch], dtype=np.float64)
        for process_id, process_name in sorted(process_names.items()):
            selected = (process_ids == process_id) & np.isfinite(values)
            if np.any(selected):
                axis.hist(
                    values[selected],
                    bins=70,
                    range=ranges[branch],
                    density=True,
                    histtype="step",
                    linewidth=1.6,
                    label=process_name,
                )
        axis.set_xlabel(title)
        axis.set_ylabel("Normalized events")
        axis.set_yscale("log")
        axis.grid(alpha=0.25)
        axis.legend()
    figure.suptitle("Unweighted anomaly-score distributions (event-representative sample)")
    figure.tight_layout()
    figure.savefig(output_dir / "score_distributions_by_process.png", dpi=180)
    plt.close(figure)


def plot_score_correlation(output_dir, global_sample):
    axo = np.asarray(global_sample.data["L1AXO_score"], dtype=np.float64)
    cicada = np.asarray(global_sample.data["L1CICADA_score"], dtype=np.float64)
    selected = np.isfinite(axo) & np.isfinite(cicada)
    figure, axis = plt.subplots(figsize=(8, 7))
    if np.any(selected):
        plot = axis.hexbin(
            axo[selected],
            cicada[selected],
            gridsize=65,
            bins="log",
            mincnt=1,
            extent=(*robust_range(axo[selected]), *robust_range(cicada[selected])),
            cmap="viridis",
        )
        figure.colorbar(plot, ax=axis, label="log sampled event density")
    axis.set_xlabel("AXOL1TL score")
    axis.set_ylabel("CICADA Model v2.2.0 score")
    axis.set_title("AXOL1TL versus CICADA (event-representative sample)")
    axis.grid(alpha=0.2)
    figure.tight_layout()
    figure.savefig(output_dir / "axo_vs_cicada.png", dpi=180)
    plt.close(figure)


def plot_closure(output_dir, process_stats, process_names, tolerance):
    process_ids = [
        process_id
        for process_id in sorted(process_names)
        if process_stats[process_id]["events"]
    ]
    labels = [process_names[process_id] for process_id in process_ids]
    positions = np.arange(len(labels))
    width = 0.36
    axo = [process_stats[process_id]["max_abs_axo_difference"] for process_id in process_ids]
    cicada = [process_stats[process_id]["max_abs_cicada_difference"] for process_id in process_ids]
    figure, axis = plt.subplots(figsize=(9, 5))
    axis.bar(positions - width / 2, axo, width, label="AXOL1TL")
    axis.bar(positions + width / 2, cicada, width, label="CICADA")
    axis.axhline(tolerance, color="black", linestyle="--", label=f"tolerance {tolerance:g}")
    axis.set_xticks(positions, labels)
    axis.set_ylabel("Maximum |recomputed - official|")
    axis.set_yscale("symlog", linthresh=max(tolerance / 10, 1.0e-15))
    axis.set_title("Two-way anomaly-score closure")
    axis.grid(axis="y", alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(output_dir / "score_closure.png", dpi=180)
    plt.close(figure)


def plot_category_shapes(output_dir, category_samples, categories, process_names):
    by_process = defaultdict(list)
    for entry in categories:
        if entry["id"] in category_samples and len(category_samples[entry["id"]].priorities):
            by_process[entry["process"]].append(entry)

    for process, entries in by_process.items():
        all_axo = np.concatenate(
            [finite(category_samples[entry["id"]].data["L1AXO_score"]) for entry in entries]
        )
        all_cicada = np.concatenate(
            [finite(category_samples[entry["id"]].data["L1CICADA_score"]) for entry in entries]
        )
        figure, axes = plt.subplots(2, 1, figsize=(10, 9))
        for entry in entries:
            sample = category_samples[entry["id"]]
            for axis, branch, value_range in (
                (axes[0], "L1AXO_score", robust_range(all_axo)),
                (axes[1], "L1CICADA_score", robust_range(all_cicada)),
            ):
                values = finite(sample.data[branch])
                if len(values):
                    axis.hist(
                        values,
                        bins=60,
                        range=value_range,
                        density=True,
                        histtype="step",
                        linewidth=1.25,
                        label=f"{entry['id']}: {entry['name']}",
                    )
        axes[0].set_xlabel("AXOL1TL score")
        axes[1].set_xlabel("CICADA Model v2.2.0 score")
        for axis in axes:
            axis.set_ylabel("Normalized events")
            axis.set_yscale("log")
            axis.grid(alpha=0.25)
            axis.legend(fontsize=7, ncol=2)
        process_id = entries[0]["process_id"]
        display_name = process_names.get(process_id, process)
        figure.suptitle(f"{display_name}: unweighted score shapes by exclusive category")
        figure.tight_layout()
        figure.savefig(output_dir / f"category_scores_{safe_name(process)}.png", dpi=180)
        plt.close(figure)

        if process == "ttbar":
            continue
        all_kinematic = np.concatenate(
            [
                finite(category_samples[entry["id"]].data["Background_kinematicValue"])
                for entry in entries
            ]
        )
        figure, axis = plt.subplots(figsize=(10, 6))
        for entry in entries:
            values = finite(category_samples[entry["id"]].data["Background_kinematicValue"])
            if len(values):
                axis.hist(
                    values,
                    bins=60,
                    range=robust_range(all_kinematic, lower=0.0, upper=0.995),
                    density=True,
                    histtype="step",
                    linewidth=1.25,
                    label=f"{entry['id']}: {entry['name']}",
                )
        axis.set_xlabel("Generator hard pT (QCD) or LHE HT (W/Z+jets) [GeV]")
        axis.set_ylabel("Normalized events")
        axis.set_yscale("log")
        axis.set_title(f"{display_name}: category-defining kinematic variable")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=7, ncol=2)
        figure.tight_layout()
        figure.savefig(output_dir / f"kinematic_shapes_{safe_name(process)}.png", dpi=180)
        plt.close(figure)


def plot_weight_shapes(output_dir, global_sample, process_names):
    process_ids = np.asarray(global_sample.data["Background_process"], dtype=np.int16)
    weights = np.asarray(global_sample.data["genWeight"], dtype=np.float64)
    log_weights = np.log10(
        np.abs(weights), where=np.abs(weights) > 0, out=np.full_like(weights, np.nan)
    )
    value_range = robust_range(log_weights)
    figure, axis = plt.subplots(figsize=(9, 6))
    for process_id, process_name in sorted(process_names.items()):
        selected = (process_ids == process_id) & np.isfinite(log_weights)
        if np.any(selected):
            axis.hist(
                log_weights[selected],
                bins=70,
                range=value_range,
                density=True,
                histtype="step",
                linewidth=1.5,
                label=process_name,
            )
    axis.set_xlabel(r"$\log_{10}(|\mathrm{genWeight}|)$")
    axis.set_ylabel("Normalized events")
    axis.set_yscale("log")
    axis.set_title("Generator-weight magnitude (event-representative sample)")
    axis.grid(alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(output_dir / "genweight_magnitudes.png", dpi=180)
    plt.close(figure)


def main():
    args = parse_args()
    if args.max_global_points <= 0 or args.max_points_per_category <= 0:
        raise SystemExit("ERROR: reservoir sizes must be positive")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    definition = json.loads(Path(args.categories).read_text())
    categories = definition["categories"]
    display_names = {"qcd": "QCD", "wjets": "W+jets", "zjets": "Z+jets", "ttbar": "TTbar"}
    process_names_by_id = {
        process_id: display_names.get(process, process)
        for process, process_id in definition["process_ids"].items()
    }
    for entry in categories:
        entry["process_id"] = definition["process_ids"][entry["process"]]
    category_names = {entry["id"]: entry["name"] for entry in categories}
    valid_category_ids = set(category_names)

    rng = np.random.default_rng(args.seed)
    global_sample = PrioritySample(args.max_global_points, rng)
    category_samples = {
        category_id: PrioritySample(args.max_points_per_category, rng)
        for category_id in valid_category_ids
    }
    overall_stats = new_stats()
    category_stats = defaultdict(new_stats)
    process_stats = defaultdict(new_stats)
    observed_samples = set()
    unknown_categories = set()

    for file_index, input_name in enumerate(args.inputs, start=1):
        print(f"[{file_index}/{len(args.inputs)}] {input_name}", flush=True)
        source = f"{input_name}:Events"
        for arrays in uproot.iterate(
            source,
            expressions=BRANCHES,
            step_size=args.step_size,
            library="np",
            how=dict,
        ):
            size = len(arrays["Background_category"])
            if not size:
                continue
            all_events = np.ones(size, dtype=bool)
            update_stats(overall_stats, arrays, all_events)
            global_sample.add(arrays)
            observed_samples.update(np.unique(arrays["Background_sample"]).astype(int).tolist())

            for process_id in np.unique(arrays["Background_process"]).astype(int):
                selected = np.asarray(arrays["Background_process"]) == process_id
                update_stats(process_stats[process_id], arrays, selected)
            for category_id in np.unique(arrays["Background_category"]).astype(int):
                selected = np.asarray(arrays["Background_category"]) == category_id
                if category_id not in valid_category_ids:
                    unknown_categories.add(category_id)
                    continue
                update_stats(category_stats[category_id], arrays, selected)
                category_samples[category_id].add(arrays, selected)

    observed_ids = sorted(
        category_id for category_id in valid_category_ids if category_stats[category_id]["events"]
    )
    if not overall_stats["events"]:
        raise SystemExit("ERROR: no events found")

    report_categories = {}
    for category_id in sorted(valid_category_ids):
        entry = dict(category_stats[category_id])
        entry["name"] = category_names[category_id]
        entry["process"] = next(
            category["process"] for category in categories if category["id"] == category_id
        )
        if entry["events"]:
            entry["mean_axo_score"] = entry["sum_axo_score"] / entry["events"]
            entry["mean_cicada_score"] = entry["sum_cicada_score"] / entry["events"]
            entry["negative_gen_weight_fraction"] = entry["negative_gen_weights"] / entry["events"]
            add_sample_quantiles(entry, category_samples[category_id])
        report_categories[str(category_id)] = entry

    report_processes = {}
    for process_id, process_name in sorted(process_names_by_id.items()):
        entry = dict(process_stats[process_id])
        entry["name"] = process_name
        if entry["events"]:
            entry["mean_axo_score"] = entry["sum_axo_score"] / entry["events"]
            entry["mean_cicada_score"] = entry["sum_cicada_score"] / entry["events"]
            entry["negative_gen_weight_fraction"] = entry["negative_gen_weights"] / entry["events"]
        report_processes[str(process_id)] = entry

    closure_ok = (
        overall_stats["max_abs_axo_difference"] <= args.score_tolerance
        and overall_stats["max_abs_cicada_difference"] <= args.score_tolerance
    )
    report = {
        "inputs": args.inputs,
        "files": len(args.inputs),
        "events": overall_stats["events"],
        "observed_parent_sample_ids": sorted(observed_samples),
        "observed_category_ids": observed_ids,
        "unobserved_category_ids": sorted(valid_category_ids - set(observed_ids)),
        "unknown_category_ids": sorted(unknown_categories),
        "score_tolerance": args.score_tolerance,
        "score_closure_ok": closure_ok,
        "normalization_note": (
            "Plots are unweighted shape/count QA. Do not compare physical process rates until each "
            "parent sample has a cross section and is normalized by its own sum of signed "
            "genWeight."
        ),
        "reservoir_note": (
            "Counts, weight sums, means, and maximum closure differences are exact. Distribution "
            "plots and reported quantiles use deterministic uniform reservoir samples."
        ),
        "overall": overall_stats,
        "processes": report_processes,
        "categories": report_categories,
    }
    report_path = output_dir / "survey_summary.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    plot_category_counts(output_dir, observed_ids, category_names, category_stats)
    plot_process_scores(output_dir, global_sample, process_names_by_id)
    plot_score_correlation(output_dir, global_sample)
    plot_closure(output_dir, process_stats, process_names_by_id, args.score_tolerance)
    plot_category_shapes(output_dir, category_samples, categories, process_names_by_id)
    plot_weight_shapes(output_dir, global_sample, process_names_by_id)

    print(f"Wrote {output_dir / 'survey_summary.json'} and PNG plots")
    print(f"Files: {len(args.inputs)}; events: {overall_stats['events']}")
    print(
        "Closure: "
        f"AXOL1TL={overall_stats['max_abs_axo_difference']:.6g}, "
        f"CICADA={overall_stats['max_abs_cicada_difference']:.6g}, "
        f"ok={closure_ok}"
    )
    print("Category counts:")
    for category_id in observed_ids:
        name = category_names[category_id]
        events = category_stats[category_id]["events"]
        print(f"  {category_id:2d} {name:30s} {events:10d}")
    if unknown_categories:
        raise SystemExit(f"ERROR: unknown category IDs found: {sorted(unknown_categories)}")
    if not closure_ok:
        raise SystemExit("ERROR: score closure exceeds tolerance")


if __name__ == "__main__":
    main()
