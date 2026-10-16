"""Small, dependency-light helpers for the background NanoAOD audit."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Iterable

import numpy as np


def effective_events(weights: np.ndarray) -> float:
    weights = np.asarray(weights, dtype=np.float64)
    weights = weights[np.isfinite(weights)]
    denominator = float(np.dot(weights, weights))
    numerator = float(weights.sum()) ** 2
    return numerator / denominator if denominator > 0 else 0.0


def weighted_ks(
    first: np.ndarray,
    second: np.ndarray,
    first_weights: np.ndarray | None = None,
    second_weights: np.ndarray | None = None,
) -> float:
    """Two-sample KS distance with non-negative sample weights."""
    first = np.asarray(first, dtype=np.float64).reshape(-1)
    second = np.asarray(second, dtype=np.float64).reshape(-1)
    if first_weights is None:
        first_weights = np.ones(first.size, dtype=np.float64)
    if second_weights is None:
        second_weights = np.ones(second.size, dtype=np.float64)
    first_weights = np.asarray(first_weights, dtype=np.float64).reshape(-1)
    second_weights = np.asarray(second_weights, dtype=np.float64).reshape(-1)
    first_ok = np.isfinite(first) & np.isfinite(first_weights) & (first_weights > 0)
    second_ok = np.isfinite(second) & np.isfinite(second_weights) & (second_weights > 0)
    first, first_weights = first[first_ok], first_weights[first_ok]
    second, second_weights = second[second_ok], second_weights[second_ok]
    if not first.size or not second.size:
        return float("nan")
    first_order = np.argsort(first, kind="stable")
    second_order = np.argsort(second, kind="stable")
    first, first_weights = first[first_order], first_weights[first_order]
    second, second_weights = second[second_order], second_weights[second_order]
    first_cdf = np.cumsum(first_weights) / first_weights.sum()
    second_cdf = np.cumsum(second_weights) / second_weights.sum()
    grid = np.unique(np.concatenate((first, second)))
    first_index = np.searchsorted(first, grid, side="right") - 1
    second_index = np.searchsorted(second, grid, side="right") - 1
    first_values = np.where(
        first_index >= 0, first_cdf[np.maximum(first_index, 0)], 0.0)
    second_values = np.where(
        second_index >= 0, second_cdf[np.maximum(second_index, 0)], 0.0)
    return float(np.max(np.abs(first_values - second_values)))


def parent_normalized_weights(
    weights: np.ndarray,
    processes: np.ndarray,
) -> tuple[np.ndarray, dict[int, float]]:
    """Normalize mean absolute generator weight to one inside each parent.

    The operation preserves all within-parent generator-weight ratios while
    removing arbitrary normalization differences between generator samples.
    It does not equalize the number of rows or the total mass of the parents.
    """
    weights = np.asarray(weights, dtype=np.float64).reshape(-1)
    processes = np.asarray(processes, dtype=np.int64).reshape(-1)
    if weights.size != processes.size:
        raise ValueError("weights and processes must have equal length")
    result = np.zeros_like(weights)
    divisors: dict[int, float] = {}
    for process in np.unique(processes):
        selected = processes == process
        finite = selected & np.isfinite(weights)
        divisor = float(np.mean(np.abs(weights[finite]))) if finite.any() else 0.0
        if not math.isfinite(divisor) or divisor <= 0:
            divisor = 1.0
        divisors[int(process)] = divisor
        result[selected] = weights[selected] / divisor
    return result, divisors


def safe_name(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_.")
    return value or "feature"


def normalized_histogram(
    values: np.ndarray,
    bins: np.ndarray,
    weights: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=np.float64).reshape(-1)
    if weights is None:
        weights = np.ones(values.size, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64).reshape(-1)
    selected = np.isfinite(values) & np.isfinite(weights)
    selected_values = values[selected]
    if selected_values.size:
        # Robust plotting bins deliberately use central quantiles. Fold the
        # true tails into the edge bins instead of silently dropping them.
        selected_values = np.clip(
            selected_values,
            np.nextafter(bins[0], bins[-1]),
            np.nextafter(bins[-1], bins[0]),
        )
    hist, _ = np.histogram(
        selected_values, bins=bins, weights=weights[selected])
    sumw2, _ = np.histogram(
        selected_values, bins=bins, weights=np.square(weights[selected]))
    normalization = float(np.sum(hist))
    if abs(normalization) <= 1e-30:
        normalization = float(np.sum(np.abs(hist)))
    if abs(normalization) <= 1e-30:
        return hist.astype(float), sumw2.astype(float)
    return hist.astype(float) / normalization, sumw2.astype(float) / normalization**2


def make_bins(values: np.ndarray, maximum_bins: int = 70) -> tuple[np.ndarray, bool]:
    """Return robust bins and whether the feature should be treated as discrete."""
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if not values.size:
        return np.linspace(0.0, 1.0, 2), False
    unique = np.unique(values)
    integer_like = np.all(np.abs(unique - np.rint(unique)) < 1e-7)
    if integer_like and unique.size <= 40:
        rounded = np.rint(unique).astype(np.int64)
        low, high = int(rounded.min()), int(rounded.max())
        if high - low <= 80:
            return np.arange(low - 0.5, high + 1.5), True
    low, high = np.quantile(values, [0.001, 0.999])
    if not high > low:
        low, high = float(values.min()), float(values.max())
    if not high > low:
        width = max(abs(float(low)) * 0.05, 1.0)
        low, high = float(low - width), float(high + width)
    margin = 0.025 * (high - low)
    return np.linspace(low - margin, high + margin, maximum_bins + 1), False


def weighted_correlation(matrix: np.ndarray, weights: np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64).reshape(-1)
    selected = np.isfinite(matrix).all(axis=1) & np.isfinite(weights) & (weights > 0)
    matrix, weights = matrix[selected], weights[selected]
    if matrix.shape[0] < 2:
        return np.full((matrix.shape[1], matrix.shape[1]), np.nan)
    weights = weights / weights.sum()
    mean = np.sum(matrix * weights[:, None], axis=0)
    centered = matrix - mean
    covariance = (centered * weights[:, None]).T @ centered
    scale = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    denominator = np.outer(scale, scale)
    return np.divide(
        covariance,
        denominator,
        out=np.zeros_like(covariance),
        where=denominator > 0,
    )


@dataclass
class FeatureStats:
    entries: int = 0
    finite: int = 0
    nonfinite: int = 0
    zero: int = 0
    minimum: float = math.inf
    maximum: float = -math.inf

    def update(self, values: np.ndarray) -> None:
        values = np.asarray(values, dtype=np.float64).reshape(-1)
        finite = np.isfinite(values)
        good = values[finite]
        self.entries += int(values.size)
        self.finite += int(good.size)
        self.nonfinite += int(values.size - good.size)
        if good.size:
            self.zero += int(np.count_nonzero(good == 0))
            self.minimum = min(self.minimum, float(good.min()))
            self.maximum = max(self.maximum, float(good.max()))

    def as_dict(self) -> dict[str, int | float | None]:
        return {
            "entries": self.entries,
            "finite": self.finite,
            "nonfinite": self.nonfinite,
            "zero": self.zero,
            "minimum": self.minimum if math.isfinite(self.minimum) else None,
            "maximum": self.maximum if math.isfinite(self.maximum) else None,
        }


class PriorityReservoir:
    """Uniform mergeable reservoir carrying labels and two weight views."""

    def __init__(self, capacity: int, rng: np.random.Generator):
        self.capacity = int(capacity)
        if self.capacity <= 0:
            raise ValueError("reservoir capacity must be positive")
        self.rng = rng
        self.priority = np.empty(0, dtype=np.float64)
        self.values = np.empty(0, dtype=np.float64)
        self.process = np.empty(0, dtype=np.int8)
        self.category = np.empty(0, dtype=np.int16)
        self.gen_weight = np.empty(0, dtype=np.float64)
        self.event_balanced_weight = np.empty(0, dtype=np.float64)

    def add(
        self,
        values: np.ndarray,
        process: np.ndarray,
        category: np.ndarray,
        gen_weight: np.ndarray,
        event_balanced_weight: np.ndarray | None = None,
    ) -> None:
        values = np.asarray(values, dtype=np.float64).reshape(-1)
        process = np.asarray(process, dtype=np.int8).reshape(-1)
        category = np.asarray(category, dtype=np.int16).reshape(-1)
        gen_weight = np.asarray(gen_weight, dtype=np.float64).reshape(-1)
        if event_balanced_weight is None:
            event_balanced_weight = gen_weight
        event_balanced_weight = np.asarray(
            event_balanced_weight, dtype=np.float64).reshape(-1)
        lengths = {
            values.size, process.size, category.size,
            gen_weight.size, event_balanced_weight.size,
        }
        if len(lengths) != 1:
            raise ValueError("reservoir arrays must have equal length")
        selected = np.isfinite(values) & np.isfinite(gen_weight)
        if not selected.any():
            return
        values, process, category = values[selected], process[selected], category[selected]
        gen_weight = gen_weight[selected]
        event_balanced_weight = event_balanced_weight[selected]
        priority = self.rng.random(values.size)
        if self.priority.size >= self.capacity and self.capacity > 0:
            competitive = priority > float(self.priority.min())
            if not competitive.any():
                return
            priority, values = priority[competitive], values[competitive]
            process, category = process[competitive], category[competitive]
            gen_weight = gen_weight[competitive]
            event_balanced_weight = event_balanced_weight[competitive]
        if values.size > self.capacity:
            keep = np.argpartition(priority, -self.capacity)[-self.capacity:]
            priority, values = priority[keep], values[keep]
            process, category = process[keep], category[keep]
            gen_weight = gen_weight[keep]
            event_balanced_weight = event_balanced_weight[keep]
        priority = np.concatenate((self.priority, priority))
        values = np.concatenate((self.values, values))
        process = np.concatenate((self.process, process))
        category = np.concatenate((self.category, category))
        gen_weight = np.concatenate((self.gen_weight, gen_weight))
        event_balanced_weight = np.concatenate(
            (self.event_balanced_weight, event_balanced_weight))
        if priority.size > self.capacity:
            keep = np.argpartition(priority, -self.capacity)[-self.capacity:]
            priority, values = priority[keep], values[keep]
            process, category = process[keep], category[keep]
            gen_weight = gen_weight[keep]
            event_balanced_weight = event_balanced_weight[keep]
        self.priority, self.values = priority, values
        self.process, self.category = process, category
        self.gen_weight = gen_weight
        self.event_balanced_weight = event_balanced_weight

    def quantiles(self, quantiles: Iterable[float]) -> list[float]:
        if not self.values.size:
            return [float("nan") for _ in quantiles]
        return np.quantile(self.values, list(quantiles)).tolist()


class AlignedReservoir:
    """Uniform reservoir for aligned event-level columns.

    Unlike one reservoir per feature, this keeps the same events for every
    column, which makes correlation matrices and multivariate follow-up plots
    meaningful.
    """

    def __init__(self, capacity: int, rng: np.random.Generator):
        self.capacity = int(capacity)
        if self.capacity <= 0:
            raise ValueError("reservoir capacity must be positive")
        self.rng = rng
        self.priority = np.empty(0, dtype=np.float64)
        self.process = np.empty(0, dtype=np.int8)
        self.category = np.empty(0, dtype=np.int16)
        self.gen_weight = np.empty(0, dtype=np.float64)
        self.columns: dict[str, np.ndarray] = {}

    def add(
        self,
        columns: dict[str, np.ndarray],
        process: np.ndarray,
        category: np.ndarray,
        gen_weight: np.ndarray,
    ) -> None:
        if not columns:
            return
        process = np.asarray(process, dtype=np.int8).reshape(-1)
        category = np.asarray(category, dtype=np.int16).reshape(-1)
        gen_weight = np.asarray(gen_weight, dtype=np.float64).reshape(-1)
        converted = {
            name: np.asarray(values, dtype=np.float64).reshape(-1)
            for name, values in columns.items()
        }
        lengths = {process.size, category.size, gen_weight.size}
        lengths.update(values.size for values in converted.values())
        if len(lengths) != 1:
            raise ValueError("aligned reservoir arrays must have equal length")
        if self.columns and set(converted) != set(self.columns):
            missing = sorted(set(self.columns) - set(converted))
            added = sorted(set(converted) - set(self.columns))
            raise ValueError(
                f"aligned columns changed between chunks; missing={missing}, added={added}")

        priority = self.rng.random(process.size)
        if self.priority.size >= self.capacity and self.capacity > 0:
            competitive = priority > float(self.priority.min())
            if not competitive.any():
                return
            priority = priority[competitive]
            process, category, gen_weight = (
                process[competitive], category[competitive], gen_weight[competitive])
            converted = {
                name: values[competitive] for name, values in converted.items()}
        if process.size > self.capacity:
            keep = np.argpartition(priority, -self.capacity)[-self.capacity:]
            priority = priority[keep]
            process, category, gen_weight = (
                process[keep], category[keep], gen_weight[keep])
            converted = {name: values[keep] for name, values in converted.items()}

        self.priority = np.concatenate((self.priority, priority))
        self.process = np.concatenate((self.process, process))
        self.category = np.concatenate((self.category, category))
        self.gen_weight = np.concatenate((self.gen_weight, gen_weight))
        if not self.columns:
            self.columns = converted
        else:
            self.columns = {
                name: np.concatenate((self.columns[name], converted[name]))
                for name in self.columns
            }

        if self.priority.size > self.capacity:
            keep = np.argpartition(
                self.priority, -self.capacity)[-self.capacity:]
            self.priority = self.priority[keep]
            self.process = self.process[keep]
            self.category = self.category[keep]
            self.gen_weight = self.gen_weight[keep]
            self.columns = {
                name: values[keep] for name, values in self.columns.items()}
