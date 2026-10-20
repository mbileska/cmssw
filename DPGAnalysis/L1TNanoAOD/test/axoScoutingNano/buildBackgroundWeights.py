#!/usr/bin/env python3
"""Build aligned all-background physics weights from a validated production.

The output NPZ stores one row per NanoAOD event in input-file order.  Its
``event_id`` is the four-column key
``(Background_sample, run, luminosityBlock, event)``.  A consumer must compare
that key with the converted tensor sample before using ``weights``.
"""

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


BRANCHES = (
    "run",
    "luminosityBlock",
    "event",
    "genWeight",
    "Background_sample",
    "Background_category",
)


def parse_args():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__, fromfile_prefix_chars="@")
    parser.add_argument(
        "inputs",
        nargs="+",
        help="validated NanoAOD files; use @filelist.txt for one path per line",
    )
    parser.add_argument(
        "--validation-report",
        required=True,
        help="successful validateTrialNano.py JSON for exactly these inputs",
    )
    parser.add_argument(
        "--reference-validation",
        help=(
            "optional earlier validation JSON (for example survey100) used to "
            "check that missing production files did not change category fractions"
        ),
    )
    parser.add_argument(
        "--samples",
        default=str(here / "background_samples_2025.json"),
        help="parent cross-section manifest",
    )
    parser.add_argument(
        "--output",
        default="production_reports/background_physics_weights.npz",
        help="aligned output NPZ",
    )
    parser.add_argument("--step-size", default="100 MB")
    parser.add_argument("--expected-files", type=int, default=0)
    parser.add_argument("--expected-events", type=int, default=0)
    parser.add_argument(
        "--max-reference-count-fraction-shift",
        type=float,
        default=0.02,
        help="largest allowed absolute category raw-fraction shift per parent",
    )
    parser.add_argument(
        "--max-reference-sumw-fraction-shift",
        type=float,
        default=0.02,
        help="largest allowed absolute category sum(genWeight)-fraction shift per parent",
    )
    return parser.parse_args()


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path):
    with Path(path).open() as handle:
        return json.load(handle)


def sample_entries(manifest):
    entries = {int(entry["id"]): entry for entry in manifest["samples"]}
    if len(entries) != len(manifest["samples"]):
        raise ValueError("duplicate parent sample ID in manifest")
    for sample_id, entry in entries.items():
        cross_section = entry.get("cross_section_pb")
        if cross_section is None or not math.isfinite(float(cross_section)) or cross_section <= 0:
            raise ValueError(f"parent sample {sample_id} lacks a positive cross section")
    return entries


def validate_report(report, entries, inputs):
    if not report.get("ok") or report.get("failures"):
        raise ValueError("validation report did not pass")
    if list(report.get("inputs", [])) != list(inputs):
        raise ValueError("input order does not exactly match the validation report")
    if report.get("events") != report.get("unique_events"):
        raise ValueError("validation report contains duplicate events")

    category_owner = {}
    for sample_id, entry in entries.items():
        for category_id in entry["allowed_category_ids"]:
            if category_id in category_owner:
                raise ValueError(f"category {category_id} belongs to multiple parents")
            category_owner[int(category_id)] = sample_id

    categories = report["categories"]
    for sample_id_text, parent in report["parent_samples"].items():
        sample_id = int(sample_id_text)
        if sample_id not in entries:
            raise ValueError(f"report parent {sample_id} is absent from manifest")
        allowed = entries[sample_id]["allowed_category_ids"]
        category_events = sum(
            int(categories.get(str(category_id), {}).get("events", 0))
            for category_id in allowed
        )
        category_sumw = sum(
            float(categories.get(str(category_id), {}).get("sum_gen_weight", 0.0))
            for category_id in allowed
        )
        if category_events != int(parent["events"]):
            raise ValueError(
                f"parent {sample_id} event closure failed: categories={category_events}, "
                f"parent={parent['events']}"
            )
        tolerance = 1.0e-9 * max(abs(float(parent["sum_gen_weight"])), 1.0)
        if abs(category_sumw - float(parent["sum_gen_weight"])) > tolerance:
            raise ValueError(
                f"parent {sample_id} sum(genWeight) closure failed: "
                f"categories={category_sumw}, parent={parent['sum_gen_weight']}"
            )

    for category_id_text in categories:
        category_id = int(category_id_text)
        if category_id not in category_owner:
            raise ValueError(f"reported category {category_id} has no parent")
    return category_owner


def compare_category_fractions(report, reference, entries, count_limit, sumw_limit):
    comparison = {
        "reference_events": int(reference["events"]),
        "production_events": int(report["events"]),
        "thresholds": {
            "max_abs_count_fraction_shift": float(count_limit),
            "max_abs_sumw_fraction_shift": float(sumw_limit),
        },
        "parents": {},
    }
    overall = True
    for sample_id, entry in sorted(entries.items()):
        parent_key = str(sample_id)
        if parent_key not in report["parent_samples"] or parent_key not in reference["parent_samples"]:
            raise ValueError(f"parent {sample_id} is absent from one comparison report")
        current_parent = report["parent_samples"][parent_key]
        reference_parent = reference["parent_samples"][parent_key]
        current_events = float(current_parent["events"])
        reference_events = float(reference_parent["events"])
        current_sumw = float(current_parent["sum_gen_weight"])
        reference_sumw = float(reference_parent["sum_gen_weight"])
        rows = []
        for category_id in entry["allowed_category_ids"]:
            current = report["categories"].get(str(category_id))
            earlier = reference["categories"].get(str(category_id))
            if current is None or earlier is None:
                raise ValueError(
                    f"category {category_id} is absent from production or reference report"
                )
            count_fraction = float(current["events"]) / current_events
            reference_count_fraction = float(earlier["events"]) / reference_events
            sumw_fraction = float(current["sum_gen_weight"]) / current_sumw
            reference_sumw_fraction = float(earlier["sum_gen_weight"]) / reference_sumw
            rows.append({
                "category_id": int(category_id),
                "name": current["name"],
                "count_fraction": count_fraction,
                "reference_count_fraction": reference_count_fraction,
                "count_fraction_shift": count_fraction - reference_count_fraction,
                "sumw_fraction": sumw_fraction,
                "reference_sumw_fraction": reference_sumw_fraction,
                "sumw_fraction_shift": sumw_fraction - reference_sumw_fraction,
            })
        max_count = max(abs(row["count_fraction_shift"]) for row in rows)
        max_sumw = max(abs(row["sumw_fraction_shift"]) for row in rows)
        passed = max_count <= count_limit and max_sumw <= sumw_limit
        overall &= passed
        comparison["parents"][parent_key] = {
            "name": entry["name"],
            "passed": bool(passed),
            "max_abs_count_fraction_shift": max_count,
            "max_abs_sumw_fraction_shift": max_sumw,
            "categories": rows,
        }
    comparison["passed"] = bool(overall)
    return comparison


def normalization_table(report, entries):
    factors = {}
    parent_rows = {}
    category_rows = []
    for sample_id, entry in sorted(entries.items()):
        parent = report["parent_samples"][str(sample_id)]
        sumw = float(parent["sum_gen_weight"])
        if not math.isfinite(sumw) or sumw <= 0.0:
            raise ValueError(f"parent {sample_id} has non-positive sum(genWeight)")
        cross_section = float(entry["cross_section_pb"])
        factor = cross_section / sumw
        factors[sample_id] = factor
        category_xsec_sum = 0.0
        for category_id in entry["allowed_category_ids"]:
            category = report["categories"][str(category_id)]
            category_xsec = factor * float(category["sum_gen_weight"])
            category_xsec_sum += category_xsec
            category_rows.append({
                "parent_sample_id": sample_id,
                "parent": entry["name"],
                "category_id": category_id,
                "category": category["name"],
                "events": int(category["events"]),
                "sum_gen_weight": float(category["sum_gen_weight"]),
                "effective_cross_section_pb": category_xsec,
                "fraction_of_parent_cross_section": category_xsec / cross_section,
            })
        if not math.isclose(category_xsec_sum, cross_section, rel_tol=1.0e-10, abs_tol=1.0e-8):
            raise ValueError(f"parent {sample_id} category cross-section closure failed")
        parent_rows[str(sample_id)] = {
            "name": entry["name"],
            "events": int(parent["events"]),
            "sum_gen_weight": sumw,
            "sum_gen_weight2": float(parent["sum_gen_weight2"]),
            "cross_section_pb": cross_section,
            "pb_per_unit_gen_weight": factor,
            "cross_section_kind": entry.get("cross_section_kind"),
            "cross_section_source": entry.get("cross_section_source"),
            "category_cross_section_sum_pb": category_xsec_sum,
        }
    return factors, parent_rows, category_rows


def write_category_csv(path, rows):
    columns = (
        "parent_sample_id",
        "parent",
        "category_id",
        "category",
        "events",
        "sum_gen_weight",
        "effective_cross_section_pb",
        "fraction_of_parent_cross_section",
    )
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main():
    import uproot

    args = parse_args()
    report = load_json(args.validation_report)
    manifest = load_json(args.samples)
    entries = sample_entries(manifest)
    category_owner = validate_report(report, entries, args.inputs)
    if args.expected_files and len(args.inputs) != args.expected_files:
        raise SystemExit(
            f"ERROR: expected {args.expected_files} input files, found {len(args.inputs)}"
        )
    if args.expected_events and int(report["events"]) != args.expected_events:
        raise SystemExit(
            f"ERROR: expected {args.expected_events} events, report has {report['events']}"
        )

    reference_comparison = None
    if args.reference_validation:
        reference = load_json(args.reference_validation)
        if not reference.get("ok") or reference.get("failures"):
            raise SystemExit("ERROR: reference validation report did not pass")
        reference_comparison = compare_category_fractions(
            report,
            reference,
            entries,
            args.max_reference_count_fraction_shift,
            args.max_reference_sumw_fraction_shift,
        )
        if not reference_comparison["passed"]:
            raise SystemExit(
                "ERROR: production/reference category fractions exceed the configured "
                "representativeness thresholds"
            )

    factors, parent_rows, category_rows = normalization_table(report, entries)
    event_id_parts = []
    weight_parts = []
    gen_weight_parts = []
    sample_parts = []
    category_parts = []
    observed_events = defaultdict(int)
    observed_sumw = defaultdict(float)
    observed_weight_mass = defaultdict(float)

    for file_index, input_name in enumerate(args.inputs, start=1):
        for arrays in uproot.iterate(
            f"{input_name}:Events",
            expressions=BRANCHES,
            step_size=args.step_size,
            library="np",
        ):
            samples = np.asarray(arrays["Background_sample"], dtype=np.int16)
            categories = np.asarray(arrays["Background_category"], dtype=np.int16)
            gen_weights = np.asarray(arrays["genWeight"], dtype=np.float64)
            if not np.all(np.isfinite(gen_weights)):
                raise RuntimeError(f"{input_name}: non-finite genWeight")
            if np.any(gen_weights < 0.0):
                raise RuntimeError(
                    f"{input_name}: negative genWeight cannot be used as a non-negative "
                    "AE/NURD loss weight without an explicit signed-weight policy"
                )
            physics_weights = np.empty(gen_weights.shape, dtype=np.float64)
            for sample_id in np.unique(samples):
                sample_id = int(sample_id)
                if sample_id not in factors:
                    raise RuntimeError(f"{input_name}: unknown parent sample {sample_id}")
                selected = samples == sample_id
                if any(category_owner.get(int(value)) != sample_id for value in np.unique(categories[selected])):
                    raise RuntimeError(f"{input_name}: category leaked across parent samples")
                physics_weights[selected] = gen_weights[selected] * factors[sample_id]
                observed_events[sample_id] += int(np.count_nonzero(selected))
                observed_sumw[sample_id] += float(np.sum(gen_weights[selected], dtype=np.float64))
                observed_weight_mass[sample_id] += float(
                    np.sum(physics_weights[selected], dtype=np.float64)
                )

            event_id = np.column_stack((
                samples.astype(np.int64),
                np.asarray(arrays["run"], dtype=np.int64),
                np.asarray(arrays["luminosityBlock"], dtype=np.int64),
                np.asarray(arrays["event"], dtype=np.int64),
            ))
            event_id_parts.append(event_id)
            weight_parts.append(physics_weights)
            gen_weight_parts.append(gen_weights)
            sample_parts.append(samples)
            category_parts.append(categories)
        if file_index % 100 == 0 or file_index == len(args.inputs):
            print(f"Read {file_index}/{len(args.inputs)} files", flush=True)

    event_id = np.concatenate(event_id_parts)
    weights = np.concatenate(weight_parts)
    gen_weight = np.concatenate(gen_weight_parts)
    background_sample = np.concatenate(sample_parts)
    background_category = np.concatenate(category_parts)
    if weights.size != int(report["events"]):
        raise RuntimeError(
            f"read {weights.size} events but validation report contains {report['events']}"
        )
    for sample_id, parent in parent_rows.items():
        numeric_id = int(sample_id)
        tolerance = 1.0e-9 * max(abs(parent["sum_gen_weight"]), 1.0)
        if observed_events[numeric_id] != parent["events"]:
            raise RuntimeError(f"parent {numeric_id} event count changed while exporting")
        if abs(observed_sumw[numeric_id] - parent["sum_gen_weight"]) > tolerance:
            raise RuntimeError(f"parent {numeric_id} sum(genWeight) changed while exporting")
        if not math.isclose(
            observed_weight_mass[numeric_id],
            parent["cross_section_pb"],
            rel_tol=1.0e-9,
            abs_tol=1.0e-6,
        ):
            raise RuntimeError(f"parent {numeric_id} physics-weight closure failed")

    output = Path(args.output)
    if output.suffix != ".npz":
        raise SystemExit("ERROR: --output must end in .npz")
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        weights=weights,
        event_id=event_id,
        gen_weight=gen_weight,
        background_sample=background_sample,
        background_category=background_category,
        weight_scope=np.asarray("all_events"),
        weight_units=np.asarray("pb"),
        schema_version=np.asarray(1, dtype=np.int16),
    )
    category_csv = output.with_name(f"{output.stem}_categories.csv")
    write_category_csv(category_csv, category_rows)
    output_report = output.with_suffix(".json")
    payload = {
        "ok": True,
        "schema_version": 1,
        "formula": "cross_section_pb * genWeight / parent_sum_gen_weight",
        "weight_scope": "all_events",
        "weight_units": "pb (equivalent to yields at 1/pb)",
        "events": int(weights.size),
        "files": len(args.inputs),
        "event_id_columns": ["Background_sample", "run", "luminosityBlock", "event"],
        "event_id_unique": int(np.unique(event_id, axis=0).shape[0]),
        "negative_weights": int(np.count_nonzero(weights < 0.0)),
        "zero_weights": int(np.count_nonzero(weights == 0.0)),
        "parents": parent_rows,
        "categories_csv": str(category_csv),
        "validation_report": str(Path(args.validation_report).resolve()),
        "validation_report_sha256": sha256_file(args.validation_report),
        "sample_manifest": str(Path(args.samples).resolve()),
        "sample_manifest_sha256": sha256_file(args.samples),
        "reference_validation": (
            str(Path(args.reference_validation).resolve())
            if args.reference_validation else None
        ),
        "reference_comparison": reference_comparison,
        "npz": str(output),
    }
    if payload["event_id_unique"] != payload["events"]:
        raise RuntimeError("event IDs became duplicated while exporting")
    output_report.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {output} ({weights.size} aligned all-background weights)")
    print(f"Wrote {output_report} and {category_csv}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (KeyError, ValueError, RuntimeError) as error:
        raise SystemExit(f"ERROR: {error}") from error
