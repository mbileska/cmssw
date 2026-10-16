import numpy as np

from auditBackgroundNano import (
    load_schema,
    parse_category_spec,
    score_residuals,
    subcategory_separation_rows,
)
from backgroundAuditCore import (
    AlignedReservoir,
    PriorityReservoir,
    effective_events,
    make_bins,
    normalized_histogram,
    parent_normalized_weights,
    weighted_correlation,
    weighted_ks,
)


def test_weighted_ks_identical_and_separated_samples():
    values = np.linspace(0.0, 1.0, 101)
    assert weighted_ks(values, values) == 0.0
    assert weighted_ks(values, values + 10.0) == 1.0


def test_parent_normalization_preserves_ratios_and_unit_mean_absolute_weight():
    weights = np.array([1.0, 3.0, 100.0, 300.0])
    processes = np.array([0, 0, 1, 1])
    normalized, divisors = parent_normalized_weights(weights, processes)
    assert divisors == {0: 2.0, 1: 200.0}
    assert np.allclose(normalized, [0.5, 1.5, 0.5, 1.5])
    for process in (0, 1):
        assert np.mean(np.abs(normalized[processes == process])) == 1.0


def test_effective_events_penalizes_concentrated_weights():
    assert effective_events(np.ones(10)) == 10.0
    assert effective_events(np.array([10.0, 0.0, 0.0])) == 1.0


def test_make_bins_recognizes_small_integer_support():
    bins, discrete = make_bins(np.array([0, 1, 1, 2, 2, 2]))
    assert discrete
    assert np.allclose(bins, [-0.5, 0.5, 1.5, 2.5])


def test_normalized_histogram_folds_underflow_and_overflow_into_edges():
    histogram, _ = normalized_histogram(
        np.array([-100.0, 0.25, 0.75, 100.0]),
        np.array([0.0, 0.5, 1.0]),
    )
    assert np.allclose(histogram, [0.5, 0.5])
    assert np.isclose(histogram.sum(), 1.0)


def test_priority_reservoir_keeps_payloads_aligned_and_bounded():
    reservoir = PriorityReservoir(20, np.random.default_rng(7))
    values = np.arange(100, dtype=float)
    reservoir.add(
        values,
        process=(values % 4).astype(int),
        category=(values % 10).astype(int),
        gen_weight=values + 1,
        event_balanced_weight=(values + 1) / 2,
    )
    assert reservoir.values.size == 20
    assert np.allclose(reservoir.gen_weight, reservoir.values + 1)
    assert np.allclose(reservoir.event_balanced_weight, (reservoir.values + 1) / 2)
    assert np.array_equal(reservoir.process, (reservoir.values % 4).astype(int))


def test_aligned_reservoir_keeps_columns_and_labels_aligned():
    reservoir = AlignedReservoir(25, np.random.default_rng(11))
    values = np.arange(100, dtype=float)
    reservoir.add(
        {"x": values, "twice_x": 2 * values},
        process=(values % 4).astype(int),
        category=(values % 8).astype(int),
        gen_weight=values + 0.5,
    )
    assert reservoir.process.size == 25
    assert np.allclose(reservoir.columns["twice_x"], 2 * reservoir.columns["x"])
    assert np.allclose(reservoir.gen_weight, reservoir.columns["x"] + 0.5)
    assert np.array_equal(
        reservoir.process, (reservoir.columns["x"] % 4).astype(int))


def test_weighted_correlation_recovers_linear_relationship():
    x = np.arange(20, dtype=float)
    matrix = np.column_stack((x, 2 * x, -x))
    correlation = weighted_correlation(matrix, np.linspace(1.0, 2.0, 20))
    assert np.allclose(np.diag(correlation), 1.0)
    assert np.isclose(correlation[0, 1], 1.0)
    assert np.isclose(correlation[0, 2], -1.0)


def test_official_score_branch_names_form_closure_residuals():
    columns = {
        "L1AXO_score": np.array([1.0, 2.0]),
        "L1AXO_scoreRecomputed": np.array([1.0, 2.0]),
        "L1CICADA_score": np.array([3.0, 4.0]),
        "L1CICADA_scoreRecomputed": np.array([3.0, 4.0]),
    }
    residuals = score_residuals(columns)
    assert set(residuals) == {
        "L1AXO_score__minus__L1AXO_scoreRecomputed",
        "L1CICADA_score__minus__L1CICADA_scoreRecomputed",
    }
    assert all(np.all(values == 0.0) for values in residuals.values())


def test_subcategory_ranking_retains_raw_and_weighted_views():
    values = np.concatenate((np.zeros(20), np.ones(20)))
    process = np.zeros(40, dtype=int)
    category = np.concatenate((np.zeros(20), np.ones(20))).astype(int)
    rows = subcategory_separation_rows(
        "ScoutingParticle_pt", values, process, category, np.ones(40),
        "candidate",
        {
            0: {"name": "low", "process": 0},
            1: {"name": "high", "process": 0},
        },
        minimum=5,
    )
    assert len(rows) == 2
    assert all(row["model_source_feature"] for row in rows)
    assert all(row["raw_ks_vs_parent_rest"] == 1.0 for row in rows)
    assert all(row["genweighted_ks_vs_parent_rest"] == 1.0 for row in rows)


def test_cmssw_category_schema_and_current_production_expectation():
    from pathlib import Path

    schema = Path(__file__).with_name("background_categories.json")
    processes, categories = load_schema(str(schema))
    assert processes == {0: "QCD", 1: "W+jets", 2: "Z+jets", 3: "ttbar"}
    assert categories[0]["process"] == 0
    assert categories[36]["process"] == 3
    expected = parse_category_spec("0-23,27-36", set(categories))
    assert expected == (set(range(24)) | set(range(27, 37)))
