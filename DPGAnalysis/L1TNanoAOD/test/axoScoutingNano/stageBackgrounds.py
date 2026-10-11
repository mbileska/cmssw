#!/usr/bin/env python3
"""Stage NanoAOD into exclusive fine-category Parquet shards without changing labels."""

import argparse
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import awkward as ak
import numpy as np
import uproot


def parse_args():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", help="validated NanoAOD ROOT files or xrootd URLs")
    parser.add_argument("--output", required=True, help="output directory")
    parser.add_argument(
        "--validation-report",
        required=True,
        help="successful JSON report from validateTrialNano.py for exactly these inputs",
    )
    parser.add_argument(
        "--categories",
        default=str(here / "background_categories.json"),
        help="exclusive-category and training-group map",
    )
    parser.add_argument(
        "--samples",
        default=str(here / "background_samples_2025.json"),
        help="parent-sample and cross-section manifest",
    )
    parser.add_argument("--step-size", default="100 MB", help="uproot iteration step size")
    return parser.parse_args()


def main():
    args = parse_args()
    output = Path(args.output)
    if output.exists() and any(output.iterdir()):
        raise SystemExit(f"ERROR: output directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    definition = json.loads(Path(args.categories).read_text())
    sample_definition = json.loads(Path(args.samples).read_text())
    validation = json.loads(Path(args.validation_report).read_text())
    if not validation.get("ok"):
        raise SystemExit(f"ERROR: validation report did not pass: {args.validation_report}")
    if validation.get("inputs") != args.inputs:
        raise SystemExit("ERROR: staging inputs do not exactly match the validation report")
    if validation.get("category_schema") != str(Path(args.categories).resolve()):
        raise SystemExit("ERROR: category schema does not match the validation report")
    if validation.get("sample_manifest") != str(Path(args.samples).resolve()):
        raise SystemExit("ERROR: sample manifest does not match the validation report")
    sample_entries = {entry["id"]: entry for entry in sample_definition["samples"]}
    category_names = {entry["id"]: entry["name"] for entry in definition["categories"]}
    counts = defaultdict(int)
    sumw = defaultdict(float)
    sumw2 = defaultdict(float)
    sample_counts = defaultdict(int)
    sample_sumw = defaultdict(float)
    sample_sumw2 = defaultdict(float)
    outputs = []

    for input_name in args.inputs:
        source_hash = hashlib.sha256(input_name.encode()).hexdigest()[:12]
        for chunk_index, arrays in enumerate(
            uproot.iterate(f"{input_name}:Events", step_size=args.step_size, library="ak")
        ):
            category_ids = ak.to_numpy(arrays["Background_category"])
            sample_ids = ak.to_numpy(arrays["Background_sample"])
            weights = ak.to_numpy(arrays["genWeight"])
            unknown_samples = set(map(int, np.unique(sample_ids))) - set(sample_entries)
            if unknown_samples:
                raise RuntimeError(f"parent samples are absent from {args.samples}: {sorted(unknown_samples)}")
            for category_id in np.unique(category_ids):
                category_id = int(category_id)
                if category_id not in category_names:
                    raise RuntimeError(f"category {category_id} is not declared in {args.categories}")
                selected = category_ids == category_id
                shard = arrays[selected]
                category_dir = output / f"{category_id:02d}_{category_names[category_id]}"
                category_dir.mkdir(exist_ok=True)
                destination = category_dir / f"part-{source_hash}-{chunk_index:06d}.parquet"
                ak.to_parquet(shard, destination)
                selected_weights = weights[selected]
                counts[category_id] += len(shard)
                sumw[category_id] += float(np.sum(selected_weights, dtype=np.float64))
                sumw2[category_id] += float(np.sum(np.square(selected_weights), dtype=np.float64))
                outputs.append(str(destination.relative_to(output)))
            for sample_id in np.unique(sample_ids):
                selected_weights = weights[sample_ids == sample_id]
                sample_id = int(sample_id)
                sample_counts[sample_id] += len(selected_weights)
                sample_sumw[sample_id] += float(np.sum(selected_weights, dtype=np.float64))
                sample_sumw2[sample_id] += float(np.sum(np.square(selected_weights), dtype=np.float64))

    manifest = {
        "schema_version": 1,
        "category_schema": str(Path(args.categories).resolve()),
        "sample_manifest": str(Path(args.samples).resolve()),
        "validation_report": str(Path(args.validation_report).resolve()),
        "inputs": args.inputs,
        "training_groups": definition["training_groups"],
        "normalization": sample_definition["normalization"],
        "categories": {
            str(category_id): {
                "name": category_names[category_id],
                "events": counts[category_id],
                "sum_gen_weight": sumw[category_id],
                "sum_gen_weight2": sumw2[category_id],
            }
            for category_id in sorted(counts)
        },
        "parent_samples": {
            str(sample_id): {
                "name": sample_entries.get(sample_id, {}).get("name"),
                "events": sample_counts[sample_id],
                "sum_gen_weight": sample_sumw[sample_id],
                "sum_gen_weight2": sample_sumw2[sample_id],
                "cross_section_pb": sample_entries.get(sample_id, {}).get("cross_section_pb"),
                "allowed_category_ids": sample_entries.get(sample_id, {}).get("allowed_category_ids"),
                "normalization_per_pb_inverse": (
                    sample_entries[sample_id]["cross_section_pb"] / sample_sumw[sample_id]
                    if sample_id in sample_entries
                    and sample_entries[sample_id].get("cross_section_pb") is not None
                    and sample_sumw[sample_id] != 0.0
                    else None
                ),
            }
            for sample_id in sorted(sample_counts)
        },
        "files": sorted(outputs),
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {len(outputs)} exclusive shards and {output / 'manifest.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
