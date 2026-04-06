"""
Базовые unit-тесты для src/consensus.py
"""

import numpy as np
import pandas as pd
import pytest

from src.consensus import (
    CONSENSUS_COLORS,
    cluster_boundaries,
    compute_consensus,
    extract_boundaries,
    save_consensus_bed,
)


# ──────────────────────────────────────────────────────────────────────────────
# Фикстуры
# ──────────────────────────────────────────────────────────────────────────────

@pytest.fixture
def simple_domains():
    """Простой набор доменов."""
    return pd.DataFrame({
        "chrom": ["chr1"] * 4,
        "start": [0, 100_000, 300_000, 600_000],
        "end":   [100_000, 300_000, 600_000, 1_000_000],
    })


@pytest.fixture
def resolution():
    return 25_000


# ──────────────────────────────────────────────────────────────────────────────
# extract_boundaries
# ──────────────────────────────────────────────────────────────────────────────

class TestExtractBoundaries:
    def test_basic(self, simple_domains, resolution):
        bnd = extract_boundaries(simple_domains, resolution)
        assert isinstance(bnd, np.ndarray)
        assert 0 in bnd
        assert 100_000 in bnd
        assert 1_000_000 in bnd

    def test_sorted(self, simple_domains, resolution):
        bnd = extract_boundaries(simple_domains, resolution)
        assert list(bnd) == sorted(bnd)

    def test_unique(self, resolution):
        df = pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "start": [0, 100_000],
            "end":   [100_000, 200_000],
        })
        bnd = extract_boundaries(df, resolution)
        # Граница 100_000 должна встречаться один раз
        assert np.sum(bnd == 100_000) == 1

    def test_empty(self, resolution):
        df = pd.DataFrame(columns=["chrom", "start", "end"])
        bnd = extract_boundaries(df, resolution)
        assert len(bnd) == 0


# ──────────────────────────────────────────────────────────────────────────────
# cluster_boundaries
# ──────────────────────────────────────────────────────────────────────────────

class TestClusterBoundaries:
    def test_no_overlap(self):
        positions = np.array([0, 100_000, 200_000])
        tolerance = 25_000
        clusters = cluster_boundaries(positions, tolerance)
        assert len(clusters) == 3

    def test_all_merge(self):
        positions = np.array([0, 10_000, 20_000])
        tolerance = 25_000
        clusters = cluster_boundaries(positions, tolerance)
        assert len(clusters) == 1
        assert 0 in clusters[0]
        assert 10_000 in clusters[0]

    def test_partial_merge(self):
        positions = np.array([0, 20_000, 300_000])
        tolerance = 25_000
        clusters = cluster_boundaries(positions, tolerance)
        assert len(clusters) == 2

    def test_single(self):
        positions = np.array([50_000])
        clusters  = cluster_boundaries(positions, 25_000)
        assert len(clusters) == 1
        assert clusters[0] == [50_000]

    def test_empty(self):
        clusters = cluster_boundaries(np.array([]), 25_000)
        assert clusters == []


# ──────────────────────────────────────────────────────────────────────────────
# compute_consensus
# ──────────────────────────────────────────────────────────────────────────────

class TestComputeConsensus:
    @pytest.fixture
    def four_algo_results(self):
        """4 алгоритма с одинаковыми границами → support=4."""
        domains = pd.DataFrame({
            "chrom": ["chr17"] * 3,
            "start": [0,        500_000, 1_000_000],
            "end":   [500_000, 1_000_000, 2_000_000],
        })
        return {
            "armatus": domains,
            "topdom":  domains,
            "scktld":  domains,
            "coitad":  domains,
        }

    def test_consensus_columns(self, four_algo_results):
        df = compute_consensus(four_algo_results, "chr17", 25_000)
        assert set(df.columns) >= {"chrom", "position", "support", "color"}

    def test_consensus_support_4(self, four_algo_results):
        df = compute_consensus(four_algo_results, "chr17", 25_000)
        assert (df["support"] == 4).any()

    def test_consensus_colors(self, four_algo_results):
        df = compute_consensus(four_algo_results, "chr17", 25_000)
        for _, row in df.iterrows():
            if row["support"] in CONSENSUS_COLORS:
                assert row["color"] == CONSENSUS_COLORS[row["support"]]

    def test_min_support_filter(self, four_algo_results):
        """support < min_support не должны попасть в результат."""
        df = compute_consensus(four_algo_results, "chr17", 25_000, min_support=3)
        assert (df["support"] >= 3).all()

    def test_two_algo_weak_consensus(self):
        d1 = pd.DataFrame({"chrom": ["chr1"], "start": [0],    "end": [500_000]})
        d2 = pd.DataFrame({"chrom": ["chr1"], "start": [0],    "end": [500_000]})
        d3 = pd.DataFrame({"chrom": ["chr1"], "start": [0],    "end": [1_000_000]})
        results = {"algo_a": d1, "algo_b": d2, "algo_c": d3}
        df = compute_consensus(results, "chr1", 50_000, min_support=2)
        assert not df.empty
        assert (df["support"] >= 2).all()

    def test_empty_results(self):
        empty = pd.DataFrame(columns=["chrom", "start", "end"])
        results = {"a": empty, "b": empty}
        df = compute_consensus(results, "chr1", 25_000)
        assert df.empty

    def test_no_consensus_when_far_apart(self):
        """Если границы далеко друг от друга, консенсуса нет."""
        d1 = pd.DataFrame({"chrom": ["chr1"], "start": [0],         "end": [500_000]})
        d2 = pd.DataFrame({"chrom": ["chr1"], "start": [10_000_000], "end": [15_000_000]})
        results = {"a": d1, "b": d2}
        df = compute_consensus(results, "chr1", 25_000, min_support=2)
        # Никакого алгоритма не поддерживает обе границы одновременно
        assert df[df["support"] >= 2].empty


# ──────────────────────────────────────────────────────────────────────────────
# Цветовая схема
# ──────────────────────────────────────────────────────────────────────────────

class TestConsensusColors:
    def test_color_scheme(self):
        assert CONSENSUS_COLORS[2] == "#FFD700"
        assert CONSENSUS_COLORS[3] == "#FF8C00"
        assert CONSENSUS_COLORS[4] == "#00C800"

    def test_no_color_for_1(self):
        assert 1 not in CONSENSUS_COLORS


# ──────────────────────────────────────────────────────────────────────────────
# Интеграционный тест
# ──────────────────────────────────────────────────────────────────────────────

class TestConsensusIntegration:
    def test_tolerance_bin_matching(self):
        """Границы в пределах ±1 bin должны кластеризоваться."""
        resolution = 50_000
        # Algo A: граница в 500_000
        # Algo B: граница в 530_000 (немного смещена, в пределах 1 bin = 50_000)
        d_a = pd.DataFrame({"chrom": ["chr1", "chr1"],
                             "start": [0,       500_000],
                             "end":   [500_000, 1_000_000]})
        d_b = pd.DataFrame({"chrom": ["chr1", "chr1"],
                             "start": [0,       530_000],
                             "end":   [530_000, 1_000_000]})
        results = {"algo_a": d_a, "algo_b": d_b}
        df = compute_consensus(results, "chr1", resolution,
                               tolerance_bins=1, min_support=2)
        # Обе границы должны кластеризоваться → support=2
        # (примерно на позиции 500_000-530_000)
        assert not df.empty
        assert (df["support"] == 2).any()