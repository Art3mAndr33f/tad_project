"""
Базовые unit-тесты для src/statistics.py
"""

import numpy as np
import pandas as pd
import pytest

from src.statistics import (
    HG19_CHROM_SIZES,
    boundary_overlap_rate,
    compare_with_reference,
    compute_basic_stats,
    compute_pairwise_matrix,
    jaccard_domains,
)


# ──────────────────────────────────────────────────────────────────────────────
# Фикстуры
# ──────────────────────────────────────────────────────────────────────────────

@pytest.fixture
def domains_a():
    return pd.DataFrame({
        "chrom": ["chr1"] * 4,
        "start": [0,      500_000, 1_000_000, 2_000_000],
        "end":   [500_000, 1_000_000, 2_000_000, 3_000_000],
    })


@pytest.fixture
def domains_b():
    return pd.DataFrame({
        "chrom": ["chr1"] * 3,
        "start": [0,      500_000, 1_500_000],
        "end":   [500_000, 1_500_000, 3_000_000],
    })


@pytest.fixture
def domains_identical(domains_a):
    return domains_a.copy()


@pytest.fixture
def resolution():
    return 50_000


# ──────────────────────────────────────────────────────────────────────────────
# compute_basic_stats
# ──────────────────────────────────────────────────────────────────────────────

class TestComputeBasicStats:
    def test_basic_n_tads(self, domains_a, resolution):
        stats = compute_basic_stats(domains_a, "chr1", resolution)
        assert stats["n_tads"] == 4

    def test_basic_median_size(self, domains_a, resolution):
        stats = compute_basic_stats(domains_a, "chr1", resolution)
        assert stats["median_size_kb"] == pytest.approx(500.0)

    def test_basic_coverage(self, domains_a, resolution):
        stats = compute_basic_stats(domains_a, "chr1", resolution)
        # 3_000_000 bp из ~249_250_621 ≈ 1.2%
        expected = 100.0 * 3_000_000 / HG19_CHROM_SIZES["chr1"]
        assert stats["coverage_pct"] == pytest.approx(expected, rel=1e-3)

    def test_empty_returns_zeros(self, resolution):
        df = pd.DataFrame(columns=["chrom", "start", "end"])
        stats = compute_basic_stats(df, "chr1", resolution)
        assert stats["n_tads"] == 0
        assert stats["coverage_pct"] == 0.0
        assert np.isnan(stats["median_size_kb"])

    def test_std_positive(self, resolution):
        df = pd.DataFrame({
            "chrom": ["chr1"] * 3,
            "start": [0,       100_000, 500_000],
            "end":   [100_000, 500_000, 1_000_000],
        })
        stats = compute_basic_stats(df, "chr1", resolution)
        assert stats["std_size_kb"] > 0


# ──────────────────────────────────────────────────────────────────────────────
# jaccard_domains
# ──────────────────────────────────────────────────────────────────────────────

class TestJaccardDomains:
    def test_identical(self, domains_a, domains_identical):
        j = jaccard_domains(domains_a, domains_identical)
        assert j == pytest.approx(1.0)

    def test_no_overlap(self, domains_a):
        df_b = pd.DataFrame({
            "chrom": ["chr2"] * 2,
            "start": [5_000_000, 6_000_000],
            "end":   [6_000_000, 7_000_000],
        })
        j = jaccard_domains(domains_a, df_b)
        assert j == pytest.approx(0.0)

    def test_partial(self, domains_a, domains_b):
        j = jaccard_domains(domains_a, domains_b)
        assert 0.0 < j < 1.0

    def test_both_empty(self):
        df = pd.DataFrame(columns=["chrom", "start", "end"])
        j = jaccard_domains(df, df)
        assert j == pytest.approx(1.0)

    def test_one_empty(self, domains_a):
        df_empty = pd.DataFrame(columns=["chrom", "start", "end"])
        j = jaccard_domains(domains_a, df_empty)
        assert j == pytest.approx(0.0)

    def test_symmetry(self, domains_a, domains_b):
        assert jaccard_domains(domains_a, domains_b) == pytest.approx(
            jaccard_domains(domains_b, domains_a)
        )


# ──────────────────────────────────────────────────────────────────────────────
# boundary_overlap_rate
# ──────────────────────────────────────────────────────────────────────────────

class TestBoundaryOverlapRate:
    def test_identical_100pct(self, domains_a, resolution):
        rate = boundary_overlap_rate(domains_a, domains_a, resolution)
        assert rate == pytest.approx(100.0)

    def test_no_overlap(self, domains_a, resolution):
        df_b = pd.DataFrame({
            "chrom": ["chr1"] * 2,
            "start": [50_000_000, 100_000_000],
            "end":   [100_000_000, 150_000_000],
        })
        rate = boundary_overlap_rate(domains_a, df_b, resolution)
        assert rate == pytest.approx(0.0)

    def test_tolerance(self, resolution):
        """Граница в пределах ±1 bin должна засчитываться."""
        d_a = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [500_000]})
        d_b = pd.DataFrame({
            "chrom": ["chr1"],
            "start": [40_000],   # в пределах 50_000 от 0
            "end":   [540_000],  # в пределах 50_000 от 500_000
        })
        rate = boundary_overlap_rate(d_a, d_b, resolution, tolerance_bins=1)
        assert rate > 0.0

    def test_empty_a_returns_nan(self, domains_b, resolution):
        df_empty = pd.DataFrame(columns=["chrom", "start", "end"])
        rate = boundary_overlap_rate(df_empty, domains_b, resolution)
        assert np.isnan(rate)

    def test_empty_b_returns_zero(self, domains_a, resolution):
        df_empty = pd.DataFrame(columns=["chrom", "start", "end"])
        rate = boundary_overlap_rate(domains_a, df_empty, resolution)
        assert rate == pytest.approx(0.0)


# ──────────────────────────────────────────────────────────────────────────────
# compute_pairwise_matrix
# ──────────────────────────────────────────────────────────────────────────────

class TestComputePairwiseMatrix:
    def test_shape(self, domains_a, domains_b, resolution):
        algo_results = {"algo_a": domains_a, "algo_b": domains_b}
        jac_df, ovl_df = compute_pairwise_matrix(algo_results, resolution)
        assert jac_df.shape == (2, 2)
        assert ovl_df.shape == (2, 2)

    def test_diagonal_jaccard(self, domains_a, resolution):
        algo_results = {"a": domains_a, "b": domains_a.copy()}
        jac_df, _ = compute_pairwise_matrix(algo_results, resolution)
        assert jac_df.loc["a", "a"] == pytest.approx(1.0)
        assert jac_df.loc["b", "b"] == pytest.approx(1.0)

    def test_diagonal_overlap(self, domains_a, resolution):
        algo_results = {"a": domains_a}
        _, ovl_df = compute_pairwise_matrix(algo_results, resolution)
        assert ovl_df.loc["a", "a"] == pytest.approx(100.0)

    def test_symmetric_jaccard(self, domains_a, domains_b, resolution):
        algo_results = {"a": domains_a, "b": domains_b}
        jac_df, _ = compute_pairwise_matrix(algo_results, resolution)
        assert jac_df.loc["a", "b"] == pytest.approx(jac_df.loc["b", "a"])


# ──────────────────────────────────────────────────────────────────────────────
# compare_with_reference
# ──────────────────────────────────────────────────────────────────────────────

class TestCompareWithReference:
    def test_perfect_recall(self, domains_a, resolution):
        """Если pred == ref, recall = 1."""
        res = compare_with_reference(domains_a, domains_a, resolution, [1])
        assert res["tol1bin"]["recall"] == pytest.approx(1.0)

    def test_perfect_precision(self, domains_a, resolution):
        res = compare_with_reference(domains_a, domains_a, resolution, [1])
        assert res["tol1bin"]["precision"] == pytest.approx(1.0)

    def test_f1_perfect(self, domains_a, resolution):
        res = compare_with_reference(domains_a, domains_a, resolution, [1])
        assert res["tol1bin"]["f1"] == pytest.approx(1.0, rel=1e-3)

    def test_zero_recall_no_overlap(self, domains_a, resolution):
        df_ref = pd.DataFrame({
            "chrom": ["chr1"],
            "start": [100_000_000],
            "end":   [200_000_000],
        })
        res = compare_with_reference(domains_a, df_ref, resolution, [1])
        assert res["tol1bin"]["recall"] == pytest.approx(0.0)

    def test_empty_pred(self, domains_a, resolution):
        df_empty = pd.DataFrame(columns=["chrom", "start", "end"])
        res = compare_with_reference(df_empty, domains_a, resolution, [1])
        assert res["tol1bin"]["recall"] == pytest.approx(0.0)

    def test_multiple_tolerances(self, domains_a, resolution):
        res = compare_with_reference(domains_a, domains_a, resolution, [1, 2])
        assert "tol1bin" in res
        assert "tol2bin" in res

    def test_tolerance_2_more_permissive(self, resolution):
        """При большем допуске recall не меньше чем при меньшем."""
        d_pred = pd.DataFrame({
            "chrom": ["chr1"],
            "start": [80_000],   # смещено на 80kb от 0
            "end":   [580_000],
        })
        d_ref = pd.DataFrame({
            "chrom": ["chr1"],
            "start": [0],
            "end":   [500_000],
        })
        res = compare_with_reference(d_pred, d_ref, resolution, [1, 2])
        # tol2bin ≥ tol1bin
        assert res["tol2bin"]["recall"] >= res["tol1bin"]["recall"]


# ──────────────────────────────────────────────────────────────────────────────
# hg19 chrom sizes
# ──────────────────────────────────────────────────────────────────────────────

class TestHg19ChromSizes:
    def test_chr1_size(self):
        assert HG19_CHROM_SIZES["chr1"] == 249_250_621

    def test_all_autosomes_present(self):
        for i in range(1, 23):
            assert f"chr{i}" in HG19_CHROM_SIZES

    def test_chrX_present(self):
        assert "chrX" in HG19_CHROM_SIZES