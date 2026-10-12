
"""Regression and integrity tests for AeroFold benchmarks."""

import csv
import math
from pathlib import Path

import pytest


BENCHMARK = (
    Path(__file__).resolve().parents[1]
    / "results"
    / "multifamily"
    / "benchmark_summary.csv"
)

EXPECTED_SPLITS = {
    "random",
    "mutation_extrapolation",
}


@pytest.fixture(scope="module")
def results():
    assert BENCHMARK.is_file(), "Benchmark CSV is missing"

    with BENCHMARK.open(newline="") as f:
        rows = list(csv.DictReader(f))

    assert rows, "Benchmark CSV is empty"
    return rows


def test_expected_columns(results):
    required = {
        "task",
        "family",
        "antigen",
        "split",
        "model",
        "spearman_mean",
        "top_recall_mean",
        "top_enrichment_mean",
    }

    assert required <= set(results[0])


def test_nine_antigen_tasks(results):
    tasks = {
        row["task"]
        for row in results
        if row["model"] == "AeroFold"
    }

    assert len(tasks) == 9


def test_both_evaluation_splits(results):
    by_task = {}

    for row in results:
        if row["model"] == "AeroFold":
            by_task.setdefault(row["task"], set()).add(
                row["split"]
            )

    assert len(by_task) == 9
    assert all(
        splits == EXPECTED_SPLITS
        for splits in by_task.values()
    )


def test_no_duplicate_benchmarks(results):
    keys = [
        (r["task"], r["split"], r["model"])
        for r in results
    ]

    assert len(keys) == len(set(keys))


@pytest.mark.parametrize(
    "column,lower,upper",
    [
        ("spearman_mean", -1.0, 1.0),
        ("top_recall_mean", 0.0, 1.0),
        ("top_enrichment_mean", 0.0, None),
    ],
)
def test_metric_ranges(results, column, lower, upper):
    for row in results:
        raw = row[column].strip()

        # Spearman is undefined for constant labels.
        if raw.lower() in {"", "nan"}:
            assert column == "spearman_mean"
            continue

        value = float(raw)

        assert math.isfinite(value), (
            f"Invalid {column} for {row['task']}"
        )
        assert value >= lower

        if upper is not None:
            assert value <= upper


@pytest.mark.parametrize(
    "metric",
    ["spearman", "top_recall", "top_enrichment"],
)
def test_summary_ordering(results, metric):
    for row in results:
        raw_values = [
            row[f"{metric}_{suffix}"]
            for suffix in ("min", "mean", "max")
        ]

        values = [float(v) for v in raw_values]

        # An undefined Spearman correlation is allowed.
        if metric == "spearman" and all(
            math.isnan(v) for v in values
        ):
            continue

        assert all(math.isfinite(v) for v in values)
        low, mean, high = values

        assert low - 1e-9 <= mean <= high + 1e-9


def test_model_scores_not_all_missing(results):
    scores = [
        float(r["spearman_mean"])
        for r in results
        if r["model"] == "AeroFold"
    ]

    assert scores
    assert any(math.isfinite(v) for v in scores)
