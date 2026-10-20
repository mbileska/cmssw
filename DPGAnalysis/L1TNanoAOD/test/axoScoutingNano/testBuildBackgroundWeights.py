import pytest

from buildBackgroundWeights import (
    compare_category_fractions,
    normalization_table,
    sample_entries,
    validate_report,
)


def manifest():
    return {
        "samples": [
            {
                "id": 0,
                "name": "first",
                "cross_section_pb": 100.0,
                "allowed_category_ids": [0, 1],
            },
            {
                "id": 1,
                "name": "second",
                "cross_section_pb": 20.0,
                "allowed_category_ids": [2],
            },
        ]
    }


def report(first_counts=(30, 10), first_sumw=(6.0, 2.0)):
    return {
        "ok": True,
        "failures": [],
        "inputs": ["one.root", "two.root"],
        "events": sum(first_counts) + 5,
        "unique_events": sum(first_counts) + 5,
        "categories": {
            "0": {
                "name": "low",
                "events": first_counts[0],
                "sum_gen_weight": first_sumw[0],
                "sum_gen_weight2": 1.0,
            },
            "1": {
                "name": "high",
                "events": first_counts[1],
                "sum_gen_weight": first_sumw[1],
                "sum_gen_weight2": 1.0,
            },
            "2": {
                "name": "only",
                "events": 5,
                "sum_gen_weight": 5.0,
                "sum_gen_weight2": 5.0,
            },
        },
        "parent_samples": {
            "0": {
                "events": sum(first_counts),
                "sum_gen_weight": sum(first_sumw),
                "sum_gen_weight2": 2.0,
            },
            "1": {
                "events": 5,
                "sum_gen_weight": 5.0,
                "sum_gen_weight2": 5.0,
            },
        },
    }


def test_category_cross_sections_close_to_each_parent():
    entries = sample_entries(manifest())
    current = report()
    owners = validate_report(current, entries, current["inputs"])
    assert owners == {0: 0, 1: 0, 2: 1}

    factors, parents, categories = normalization_table(current, entries)
    assert factors == pytest.approx({0: 12.5, 1: 4.0})
    assert parents["0"]["category_cross_section_sum_pb"] == pytest.approx(100.0)
    assert parents["1"]["category_cross_section_sum_pb"] == pytest.approx(20.0)
    category_xsecs = {
        row["category_id"]: row["effective_cross_section_pb"] for row in categories
    }
    assert category_xsecs == pytest.approx({0: 75.0, 1: 25.0, 2: 20.0})


def test_validation_rejects_category_parent_nonclosure():
    entries = sample_entries(manifest())
    broken = report()
    broken["categories"]["1"]["sum_gen_weight"] = 1.0
    with pytest.raises(ValueError, match=r"sum\(genWeight\) closure"):
        validate_report(broken, entries, broken["inputs"])


def test_reference_comparison_checks_raw_and_weighted_category_fractions():
    entries = sample_entries(manifest())
    current = report()
    compatible = report(first_counts=(301, 99), first_sumw=(60.1, 19.9))
    comparison = compare_category_fractions(
        current, compatible, entries, count_limit=0.01, sumw_limit=0.01
    )
    assert comparison["passed"]

    shifted = report(first_counts=(20, 20), first_sumw=(4.0, 4.0))
    comparison = compare_category_fractions(
        current, shifted, entries, count_limit=0.01, sumw_limit=0.01
    )
    assert not comparison["passed"]
    assert not comparison["parents"]["0"]["passed"]
