#!/usr/bin/env python3
"""Validate score closure, exclusive categories, event identity, and generator weights."""

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import awkward as ak
import numpy as np
import uproot


REQUIRED_BRANCHES = {
    "run",
    "luminosityBlock",
    "event",
    "genWeight",
    "L1AXO_score",
    "L1AXO_scoreRecomputed",
    "L1AXO_scoreDifference",
    "nL1AXOInput",
    "L1AXOInput_value2",
    "L1CICADA_score",
    "L1CICADA_scoreRecomputed",
    "L1CICADA_scoreDifference",
    "nL1CICADAInput",
    "L1CICADAInput_et",
    "L1CICADAInput_gctEta",
    "L1CICADAInput_gctPhi",
    "L1CICADAInput_present",
    "Background_process",
    "Background_sample",
    "Background_decay",
    "Background_kinematicBin",
    "Background_category",
    "Background_kinematicValue",
}


def parse_args():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", help="NanoAOD ROOT files or xrootd URLs")
    parser.add_argument(
        "--categories",
        default=str(here / "background_categories.json"),
        help="category schema JSON",
    )
    parser.add_argument(
        "--samples",
        default=str(here / "background_samples_2025.json"),
        help="parent-sample and cross-section manifest JSON",
    )
    parser.add_argument("--output", default="trial_validation.json", help="JSON report")
    parser.add_argument("--score-tolerance", type=float, default=1.0e-6)
    parser.add_argument(
        "--expected-categories",
        default="",
        help="required category IDs: 'all' or comma-separated IDs/ranges such as 0-23,27-36",
    )
    parser.add_argument(
        "--min-events-per-populated-category",
        type=int,
        default=0,
        help="fail when a populated (or explicitly expected) category has fewer events",
    )
    parser.add_argument(
        "--require-cross-sections",
        action="store_true",
        help="fail unless every parent sample present has a cross section in the sample manifest",
    )
    parser.add_argument("--step-size", default="100 MB", help="uproot iteration step size")
    return parser.parse_args()


def parse_expected_categories(specification, valid_ids):
    if not specification:
        return set()
    if specification.lower() == "all":
        return set(valid_ids)
    result = set()
    for item in specification.split(","):
        item = item.strip()
        if not item:
            continue
        if "-" in item:
            lower_text, upper_text = item.split("-", 1)
            lower, upper = int(lower_text), int(upper_text)
            if lower > upper:
                raise ValueError(f"invalid descending category range {item}")
            result.update(range(lower, upper + 1))
        else:
            result.add(int(item))
    unknown = result - valid_ids
    if unknown:
        raise ValueError(f"unknown expected category IDs: {sorted(unknown)}")
    return result


def category_is_consistent(arrays, definition):
    process_ids = definition["process_ids"]
    categories = {entry["id"]: entry for entry in definition["categories"]}
    result = np.ones(len(arrays), dtype=bool)
    ids = ak.to_numpy(arrays["Background_category"])
    processes = ak.to_numpy(arrays["Background_process"])
    decays = ak.to_numpy(arrays["Background_decay"]).astype(np.int16)
    bins = ak.to_numpy(arrays["Background_kinematicBin"]).astype(np.int16)
    values = ak.to_numpy(arrays["Background_kinematicValue"])
    expected_ids = np.full(len(ids), -1, dtype=np.int16)
    expected_ids[processes == 0] = bins[processes == 0]
    expected_ids[processes == 1] = 6 + (decays[processes == 1] - 1) * 3 + bins[processes == 1]
    expected_ids[processes == 2] = 15 + (decays[processes == 2] - 1) * 3 + bins[processes == 2]
    expected_ids[processes == 3] = 27 + decays[processes == 3]
    result &= ids == expected_ids
    for category_id in np.unique(ids):
        selected = ids == category_id
        entry = categories.get(int(category_id))
        if entry is None:
            result[selected] = False
            continue
        result[selected] &= processes[selected] == process_ids[entry["process"]]
        bounds = entry["kinematic_bin"]
        if bounds is not None:
            result[selected] &= values[selected] >= bounds[0]
            if bounds[1] is not None:
                result[selected] &= values[selected] < bounds[1]
    return result


def main():
    args = parse_args()
    definition = json.loads(Path(args.categories).read_text())
    valid_ids = {entry["id"] for entry in definition["categories"]}
    try:
        expected_ids = parse_expected_categories(args.expected_categories, valid_ids)
    except ValueError as error:
        raise SystemExit(f"ERROR: {error}") from error
    sample_definition = json.loads(Path(args.samples).read_text())
    sample_entries = {entry["id"]: entry for entry in sample_definition["samples"]}
    if len(sample_entries) != len(sample_definition["samples"]):
        raise SystemExit("ERROR: duplicate parent sample ID in sample manifest")
    process_ids = definition["process_ids"]
    failures = []
    seen_events = set()
    counts = defaultdict(int)
    sumw = defaultdict(float)
    sumw2 = defaultdict(float)
    sample_counts = defaultdict(int)
    sample_sumw = defaultdict(float)
    sample_sumw2 = defaultdict(float)
    sample_files = defaultdict(int)
    events_total = 0
    max_axo_difference = 0.0
    max_cicada_difference = 0.0
    runs_sumw_total = 0.0

    for input_name in args.inputs:
        runs_sumw = None
        with uproot.open(input_name) as root_file:
            if "Events" not in root_file:
                failures.append(f"{input_name}: missing Events tree")
                continue
            available = set(root_file["Events"].keys())
            missing = sorted(REQUIRED_BRANCHES - available)
            if missing:
                failures.append(f"{input_name}: missing branches: {', '.join(missing)}")
                continue
            if "Runs" not in root_file or "genEventSumw" not in root_file["Runs"]:
                failures.append(f"{input_name}: missing Runs/genEventSumw")
            else:
                runs_sumw = float(np.sum(root_file["Runs"]["genEventSumw"].array(library="np"), dtype=np.float64))
                runs_sumw_total += runs_sumw

        source = f"{input_name}:Events"
        file_event_sumw = 0.0
        file_event_sum_absw = 0.0
        file_sample_ids = set()
        for arrays in uproot.iterate(source, expressions=sorted(REQUIRED_BRANCHES), step_size=args.step_size):
            chunk_size = len(arrays)
            events_total += chunk_size
            axo_score = ak.to_numpy(arrays["L1AXO_score"])
            axo_recomputed = ak.to_numpy(arrays["L1AXO_scoreRecomputed"])
            axo_stored_difference = ak.to_numpy(arrays["L1AXO_scoreDifference"])
            cicada_score = ak.to_numpy(arrays["L1CICADA_score"])
            cicada_recomputed = ak.to_numpy(arrays["L1CICADA_scoreRecomputed"])
            cicada_stored_difference = ak.to_numpy(arrays["L1CICADA_scoreDifference"])
            score_fields = (
                axo_score,
                axo_recomputed,
                axo_stored_difference,
                cicada_score,
                cicada_recomputed,
                cicada_stored_difference,
            )
            if not all(np.all(np.isfinite(field)) for field in score_fields):
                failures.append(f"{input_name}: non-finite anomaly score or score difference")
            axo_actual_difference = axo_recomputed - axo_score
            cicada_actual_difference = cicada_recomputed - cicada_score
            axo_difference = np.abs(axo_actual_difference)
            cicada_difference = np.abs(cicada_actual_difference)
            finite_axo_difference = axo_difference[np.isfinite(axo_difference)]
            finite_cicada_difference = cicada_difference[np.isfinite(cicada_difference)]
            max_axo_difference = max(
                max_axo_difference, float(np.max(finite_axo_difference, initial=0.0))
            )
            max_cicada_difference = max(
                max_cicada_difference, float(np.max(finite_cicada_difference, initial=0.0))
            )
            if np.any(axo_difference > args.score_tolerance):
                failures.append(f"{input_name}: AXOL1TL closure exceeds tolerance")
            if np.any(cicada_difference > args.score_tolerance):
                failures.append(f"{input_name}: CICADA closure exceeds tolerance")
            if np.any(np.abs(axo_stored_difference - axo_actual_difference) > args.score_tolerance):
                failures.append(f"{input_name}: stored AXOL1TL difference is inconsistent with the two scores")
            if np.any(np.abs(cicada_stored_difference - cicada_actual_difference) > args.score_tolerance):
                failures.append(f"{input_name}: stored CICADA difference is inconsistent with the two scores")

            if np.any(ak.to_numpy(arrays["nL1AXOInput"]) != 57):
                failures.append(f"{input_name}: AXOL1TL tensor does not have 57 values for every event")
            if np.any(ak.to_numpy(arrays["nL1CICADAInput"]) != 252):
                failures.append(f"{input_name}: CICADA tensor does not have 252 values for every event")
            present = ak.to_numpy(ak.sum(arrays["L1CICADAInput_present"], axis=1))
            if np.any(present != 252):
                failures.append(f"{input_name}: CICADA region grid is incomplete or duplicated")
            expected_eta = np.tile(np.arange(4, 18, dtype=np.uint8), 18)
            expected_phi = np.repeat(np.arange(18, dtype=np.uint8), 14)
            eta_grid = ak.to_numpy(arrays["L1CICADAInput_gctEta"])
            phi_grid = ak.to_numpy(arrays["L1CICADAInput_gctPhi"])
            if np.any(eta_grid != expected_eta) or np.any(phi_grid != expected_phi):
                failures.append(f"{input_name}: CICADA tensor coordinates are not in canonical phi-major order")
            cicada_et = ak.to_numpy(arrays["L1CICADAInput_et"])
            if np.any(cicada_et > 1023):
                failures.append(f"{input_name}: CICADA tensor contains an ET outside the 10-bit input range")

            weights = ak.to_numpy(arrays["genWeight"])
            if not np.all(np.isfinite(weights)):
                failures.append(f"{input_name}: non-finite genWeight")
            file_event_sumw += float(np.sum(weights, dtype=np.float64))
            file_event_sum_absw += float(np.sum(np.abs(weights), dtype=np.float64))
            category_ids = ak.to_numpy(arrays["Background_category"])
            if not set(np.unique(category_ids)).issubset(valid_ids):
                failures.append(f"{input_name}: category outside declared schema")
            if not np.all(category_is_consistent(arrays, definition)):
                failures.append(f"{input_name}: category/process/kinematic fields are inconsistent")

            samples = ak.to_numpy(arrays["Background_sample"])
            file_sample_ids.update(map(int, np.unique(samples)))
            processes = ak.to_numpy(arrays["Background_process"])
            for sample_id in np.unique(samples):
                sample_id = int(sample_id)
                sample_entry = sample_entries.get(sample_id)
                if sample_entry is None:
                    failures.append(f"{input_name}: parent sample {sample_id} is absent from {args.samples}")
                    continue
                expected_process = process_ids[sample_entry["process"]]
                selected_sample = samples == sample_id
                if np.any(processes[selected_sample] != expected_process):
                    failures.append(
                        f"{input_name}: parent sample {sample_id} does not match its configured process "
                        f"{sample_entry['process']}"
                    )
                allowed_categories = set(sample_entry["allowed_category_ids"])
                observed_categories = set(map(int, np.unique(category_ids[selected_sample])))
                if not observed_categories.issubset(allowed_categories):
                    failures.append(
                        f"{input_name}: parent sample {sample_id} leaked outside its allowed categories: "
                        f"{sorted(observed_categories - allowed_categories)}"
                    )
            runs = ak.to_numpy(arrays["run"])
            lumis = ak.to_numpy(arrays["luminosityBlock"])
            event_ids = ak.to_numpy(arrays["event"])
            for event_key in zip(samples.tolist(), runs.tolist(), lumis.tolist(), event_ids.tolist()):
                if event_key in seen_events:
                    failures.append(f"duplicate event {event_key}")
                seen_events.add(event_key)

            for category_id in np.unique(category_ids):
                selected_weights = weights[category_ids == category_id]
                category_id = int(category_id)
                counts[category_id] += len(selected_weights)
                sumw[category_id] += float(np.sum(selected_weights, dtype=np.float64))
                sumw2[category_id] += float(np.sum(np.square(selected_weights), dtype=np.float64))
            for sample_id in np.unique(samples):
                selected_weights = weights[samples == sample_id]
                sample_id = int(sample_id)
                sample_counts[sample_id] += len(selected_weights)
                sample_sumw[sample_id] += float(np.sum(selected_weights, dtype=np.float64))
                sample_sumw2[sample_id] += float(np.sum(np.square(selected_weights), dtype=np.float64))

        if runs_sumw is not None:
            tolerance = 1.0e-5 * max(file_event_sum_absw, 1.0)
            if abs(file_event_sumw - runs_sumw) > tolerance:
                failures.append(
                    f"{input_name}: Events/genWeight sum {file_event_sumw} disagrees with "
                    f"Runs/genEventSumw {runs_sumw}"
                )
        if len(file_sample_ids) != 1:
            failures.append(
                f"{input_name}: expected exactly one parent sample per output file, found {sorted(file_sample_ids)}"
            )
        else:
            sample_files[next(iter(file_sample_ids))] += 1

    if events_total == 0:
        failures.append("no events were read")

    count_check_ids = expected_ids if expected_ids else set(counts)
    for category_id in sorted(count_check_ids):
        count = counts[category_id]
        if count < args.min_events_per_populated_category:
            failures.append(
                f"category {category_id} has {count} events, below {args.min_events_per_populated_category}"
            )
    missing_expected = expected_ids - set(counts)
    if missing_expected:
        failures.append(f"expected categories are absent: {sorted(missing_expected)}")

    missing_cross_sections = sorted(
        sample_id
        for sample_id in sample_counts
        if sample_id not in sample_entries or sample_entries[sample_id].get("cross_section_pb") is None
    )
    if args.require_cross_sections and missing_cross_sections:
        failures.append(f"parent samples lack cross sections: {missing_cross_sections}")

    names = {entry["id"]: entry["name"] for entry in definition["categories"]}

    def parent_report(sample_id):
        entry = sample_entries.get(sample_id, {})
        allowed = entry.get("allowed_category_ids", [])
        least_category = min(allowed, key=lambda category_id: counts[category_id]) if allowed else None
        least_count = counts[least_category] if least_category is not None else 0
        if args.min_events_per_populated_category > 0 and least_count > 0:
            scale_factor = max(1.0, args.min_events_per_populated_category / least_count)
            estimated_files = math.ceil(sample_files[sample_id] * scale_factor)
        else:
            scale_factor = None
            estimated_files = None
        return {
            "name": entry.get("name"),
            "output_files": sample_files[sample_id],
            "events": sample_counts[sample_id],
            "sum_gen_weight": sample_sumw[sample_id],
            "sum_gen_weight2": sample_sumw2[sample_id],
            "cross_section_pb": entry.get("cross_section_pb"),
            "allowed_category_ids": allowed,
            "least_populated_category": least_category,
            "least_populated_category_events": least_count,
            "scale_factor_to_target_at_current_mix": scale_factor,
            "estimated_total_output_files_at_current_mix": estimated_files,
        }

    report = {
        "ok": not failures,
        "inputs": args.inputs,
        "category_schema": str(Path(args.categories).resolve()),
        "sample_manifest": str(Path(args.samples).resolve()),
        "events": events_total,
        "unique_events": len(seen_events),
        "expected_categories": sorted(expected_ids),
        "max_abs_score_difference": {
            "AXOL1TL": max_axo_difference,
            "CICADA": max_cicada_difference,
        },
        "categories": {
            str(category_id): {
                "name": names[category_id],
                "events": counts[category_id],
                "sum_gen_weight": sumw[category_id],
                "sum_gen_weight2": sumw2[category_id],
            }
            for category_id in sorted(counts)
        },
        "parent_samples": {str(sample_id): parent_report(sample_id) for sample_id in sorted(sample_counts)},
        "normalization": {
            "formula": sample_definition["normalization"],
            "events_sum_gen_weight": float(sum(sample_sumw.values())),
            "runs_sum_gen_weight": runs_sumw_total,
            "parents_missing_cross_section": missing_cross_sections,
            "ready": not missing_cross_sections,
        },
        "failures": sorted(set(failures)),
    }
    Path(args.output).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
