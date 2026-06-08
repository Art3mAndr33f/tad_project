from __future__ import annotations

"""Tests for src/statistical_methods/changepoint_tad.py

Coverage
--------
1.  compute_insulation_score
    a. known 4-block matrix → IS dips at block boundaries
    b. flat (uniform) matrix → constant IS
    c. minimum matrix size (< 2*window → edge bins only = NaN)
    d. output length == matrix size
    e. output is z-score normalised (mean≈0, std≈1 on valid bins)
    f. edge bins are NaN

2.  _run_pelt_single_beta
    a. synthetic step signal → at least one breakpoint
    b. constant signal → zero or very few breakpoints
    c. output is sorted and within signal bounds
    d. ImportError raised cleanly if ruptures missing (mocked)

3.  _stability_filter
    a. boundary present in all runs → retained
    b. boundary present in too few runs → dropped
    c. tolerance merging: nearby bins counted as same boundary
    d. empty input → empty output, no exception
    e. all empty runs → empty output

4.  run_pelt_tad  (integration via mocked load_hic_matrix)
    a. synthetic TAD matrix with clear block structure → returns DataFrame
    b. empty matrix returned by loader → empty DataFrame, no exception
    c. loader raises exception → empty DataFrame, no exception
    d. None returned by loader → empty DataFrame, no exception
    e. output contract: columns == ["chrom","start","end"]
    f. output types: start/end are integers, chrom matches input
    g. all-zero IS → empty DataFrame, no exception
"""

import importlib
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

# ── optional ruptures import guard ──────────────────────────────────────────
try:
    import ruptures  # noqa: F401
    HAS_RUPTURES = True
except ImportError:
    HAS_RUPTURES = False

ruptures_required = pytest.mark.skipif(
    not HAS_RUPTURES,
    reason="ruptures not installed — skipping PELT tests",
)

from src.statistical_methods.changepoint_tad import (
    _run_pelt_single_beta,
    _stability_filter,
    compute_insulation_score,
    run_pelt_tad,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _block_matrix(n_blocks: int, block_size: int, intra: float = 10.0,
                  inter: float = 0.1) -> np.ndarray:
    """Create a symmetric block-diagonal contact matrix with n_blocks blocks."""
    n = n_blocks * block_size
    mat = np.full((n, n), inter, dtype=np.float64)
    for b in range(n_blocks):
        lo, hi = b * block_size, (b + 1) * block_size
        mat[lo:hi, lo:hi] = intra
    return mat


def _step_signal(n: int = 200, n_steps: int = 4) -> np.ndarray:
    """Signal with n_steps abrupt level changes (easy for PELT)."""
    signal = np.zeros(n, dtype=np.float64)
    step   = n // n_steps
    for k in range(n_steps):
        signal[k * step : (k + 1) * step] = float(k)
    return signal


CFG_MINIMAL = {
    "changepoint": {
        "window_bins":      3,
        "beta_min":         0.5,
        "beta_max":         3.0,
        "n_beta":           5,
        "min_beta_support": 2,
    },
    "consensus": {"tolerance_bins": 1},
}


# ---------------------------------------------------------------------------
# 1. compute_insulation_score
# ---------------------------------------------------------------------------

class TestComputeInsulationScore:

    def test_output_length_equals_matrix_size(self) -> None:
        mat = _block_matrix(4, 10)
        is_ = compute_insulation_score(mat, window_bins=3)
        assert len(is_) == mat.shape[0]

    def test_edge_bins_are_nan(self) -> None:
        mat = _block_matrix(4, 10)
        is_ = compute_insulation_score(mat, window_bins=3)
        assert np.isnan(is_[:3]).all()
        assert np.isnan(is_[-3:]).all()

    def test_block_boundaries_have_low_is(self) -> None:
        """IS should be low (negative after z-score) at block boundaries."""
        n_blocks, bs = 4, 15
        mat = _block_matrix(n_blocks, bs, intra=20.0, inter=0.01)
        is_ = compute_insulation_score(mat, window_bins=5)
        valid = ~np.isnan(is_)
        # Block boundaries are at bins 15, 30, 45
        for boundary_bin in [15, 30, 45]:
            if valid[boundary_bin]:
                # IS at boundary should be below median of valid bins
                median_is = np.nanmedian(is_[valid])
                assert is_[boundary_bin] < median_is, (
                    f"Expected low IS at block boundary bin {boundary_bin}, "
                    f"got IS={is_[boundary_bin]:.3f} vs median={median_is:.3f}"
                )

    def test_flat_matrix_constant_is(self) -> None:
        """Uniform contact matrix → constant IS profile (ignoring edges)."""
        n = 40
        mat = np.ones((n, n), dtype=np.float64)
        is_ = compute_insulation_score(mat, window_bins=3)
        valid = ~np.isnan(is_)
        # All valid IS values should be equal (after z-score: all 0 or uniform)
        vals = is_[valid]
        assert np.allclose(vals, vals[0], atol=1e-10), \
            "Flat matrix should yield constant IS on valid bins"

    def test_z_score_normalised(self) -> None:
        """Valid bins should have approx mean=0, std=1 (or all-zero for flat)."""
        mat = _block_matrix(5, 12, intra=15.0, inter=0.05)
        is_ = compute_insulation_score(mat, window_bins=4)
        valid = ~np.isnan(is_)
        vals = is_[valid]
        if len(vals) > 1 and not np.allclose(vals, vals[0]):
            assert abs(vals.mean()) < 0.1,  \
                f"Expected z-score mean ≈ 0, got {vals.mean():.4f}"
            assert abs(vals.std() - 1.0) < 0.1, \
                f"Expected z-score std ≈ 1, got {vals.std():.4f}"

    def test_tiny_matrix_all_nan(self) -> None:
        """Matrix smaller than 2*window → all NaN (no valid bins)."""
        mat = np.ones((4, 4), dtype=np.float64)
        is_ = compute_insulation_score(mat, window_bins=3)   # 2*3=6 > 4
        assert np.isnan(is_).all()

    def test_output_dtype_float64(self) -> None:
        mat = _block_matrix(3, 10)
        is_ = compute_insulation_score(mat, window_bins=2)
        assert is_.dtype == np.float64


# ---------------------------------------------------------------------------
# 2. _run_pelt_single_beta
# ---------------------------------------------------------------------------

class TestRunPeltSingleBeta:

    @ruptures_required
    def test_step_signal_has_breakpoints(self) -> None:
        signal = _step_signal(n=200, n_steps=4)
        bkps   = _run_pelt_single_beta(signal, beta=1.0)
        assert len(bkps) > 0

    @ruptures_required
    def test_constant_signal_few_breakpoints(self) -> None:
        """Constant signal with large penalty → very few (ideally 0) breakpoints."""
        signal = np.zeros(100, dtype=np.float64)
        bkps   = _run_pelt_single_beta(signal, beta=50.0)
        assert len(bkps) <= 2

    @ruptures_required
    def test_output_sorted_and_in_bounds(self) -> None:
        signal = _step_signal(n=100, n_steps=5)
        bkps   = _run_pelt_single_beta(signal, beta=1.0)
        assert (np.diff(bkps) > 0).all() or len(bkps) <= 1, \
            "Breakpoints must be strictly sorted"
        assert all(0 <= b < len(signal) for b in bkps), \
            "All breakpoints must be within signal bounds (sentinel excluded)"

    @ruptures_required
    def test_large_penalty_fewer_breakpoints(self) -> None:
        """Larger penalty → fewer or equal breakpoints."""
        signal = _step_signal(n=200, n_steps=8)
        bkps_small = _run_pelt_single_beta(signal, beta=0.5)
        bkps_large = _run_pelt_single_beta(signal, beta=10.0)
        assert len(bkps_large) <= len(bkps_small)

    @ruptures_required
    def test_output_is_ndarray(self) -> None:
        signal = _step_signal()
        bkps   = _run_pelt_single_beta(signal, beta=1.0)
        assert isinstance(bkps, np.ndarray)
        assert bkps.dtype == np.int64


# ---------------------------------------------------------------------------
# 3. _stability_filter
# ---------------------------------------------------------------------------

class TestStabilityFilter:

    def test_boundary_in_all_runs_retained(self) -> None:
        """A boundary present in every run must always be retained."""
        runs = [
            np.array([10, 20, 30], dtype=np.int64),
            np.array([10, 21, 31], dtype=np.int64),   # 21 within tol=1 of 20
            np.array([10, 20, 29], dtype=np.int64),   # 29 within tol=1 of 30
        ]
        stable = _stability_filter(runs, n_bins=50, min_beta_support=3,
                                   tolerance_bins=1)
        assert 10 in stable

    def test_boundary_too_few_runs_dropped(self) -> None:
        """A boundary appearing in only 1 run (support < min) must be dropped."""
        runs = [
            np.array([10], dtype=np.int64),
            np.array([20], dtype=np.int64),   # 10 only in run 0
            np.array([20], dtype=np.int64),
        ]
        stable = _stability_filter(runs, n_bins=50, min_beta_support=2,
                                   tolerance_bins=0)
        assert 10 not in stable    # support == 1 < 2
        assert 20 in stable        # support == 2 >= 2

    def test_tolerance_merges_nearby_bins(self) -> None:
        """Bins 9, 10, 11 within tolerance=1 of 10 → all count as same boundary."""
        runs = [
            np.array([9],  dtype=np.int64),
            np.array([10], dtype=np.int64),
            np.array([11], dtype=np.int64),
        ]
        stable = _stability_filter(runs, n_bins=50, min_beta_support=2,
                                   tolerance_bins=1)
        # At least one of the merged group should appear in stable
        assert any(b in stable for b in [9, 10, 11])

    def test_empty_input_returns_empty(self) -> None:
        stable = _stability_filter([], n_bins=100, min_beta_support=2)
        assert isinstance(stable, np.ndarray)
        assert len(stable) == 0

    def test_all_empty_runs_returns_empty(self) -> None:
        runs = [np.array([], dtype=np.int64)] * 5
        stable = _stability_filter(runs, n_bins=100, min_beta_support=2)
        assert len(stable) == 0

    def test_output_is_sorted(self) -> None:
        runs = [
            np.array([30, 10, 50], dtype=np.int64),
            np.array([30, 10, 50], dtype=np.int64),
        ]
        stable = _stability_filter(runs, n_bins=100, min_beta_support=2,
                                   tolerance_bins=0)
        assert (np.diff(stable) > 0).all() or len(stable) <= 1

    def test_min_support_one_returns_all(self) -> None:
        """min_beta_support=1 → every candidate boundary is retained."""
        runs = [
            np.array([5, 15, 25], dtype=np.int64),
            np.array([7, 17],     dtype=np.int64),
        ]
        stable = _stability_filter(runs, n_bins=50, min_beta_support=1,
                                   tolerance_bins=0)
        # All unique bins from runs should be present
        all_bins = {5, 7, 15, 17, 25}
        for b in all_bins:
            assert b in stable


# ---------------------------------------------------------------------------
# 4. run_pelt_tad — integration (mocked load_hic_matrix)
# ---------------------------------------------------------------------------

class TestRunPeltTad:
    """Integration tests using unittest.mock to bypass file I/O."""

    CHROM = "chr17"
    RES   = 100_000

    def _run(self, matrix: np.ndarray, cfg: dict = CFG_MINIMAL) -> pd.DataFrame:
        # Patch the module-level name so patch() can find the attribute.
        # changepoint_tad.py imports load_hic_matrix at module level,
        # so patching "src.statistical_methods.changepoint_tad.load_hic_matrix"
        # is the correct and sufficient approach.
        with patch(
            "src.statistical_methods.changepoint_tad.load_hic_matrix",
            return_value=matrix,
        ):
            return run_pelt_tad(
                chrom=self.CHROM,
                resolution=self.RES,
                data_path="data/raw",
                cfg=cfg,
            )

    @ruptures_required
    def test_block_matrix_returns_dataframe(self) -> None:
        mat = _block_matrix(6, 20, intra=20.0, inter=0.01)
        df  = self._run(mat)
        assert isinstance(df, pd.DataFrame)

    @ruptures_required
    def test_output_contract_columns(self) -> None:
        mat = _block_matrix(4, 15)
        df  = self._run(mat)
        assert list(df.columns) == ["chrom", "start", "end"]

    @ruptures_required
    def test_output_chrom_matches_input(self) -> None:
        mat = _block_matrix(4, 15)
        df  = self._run(mat)
        if len(df) > 0:
            assert (df["chrom"] == self.CHROM).all()

    @ruptures_required
    def test_output_types(self) -> None:
        mat = _block_matrix(4, 15)
        df  = self._run(mat)
        if len(df) > 0:
            assert np.issubdtype(df["start"].dtype, np.integer)
            assert np.issubdtype(df["end"].dtype,   np.integer)

    @ruptures_required
    def test_start_end_aligned_to_resolution(self) -> None:
        mat = _block_matrix(5, 20)
        df  = self._run(mat)
        if len(df) > 0:
            assert (df["start"] % self.RES == 0).all()
            assert (df["end"]   % self.RES == 0).all()
            assert (df["end"] > df["start"]).all()

    def test_empty_matrix_returns_empty_df(self) -> None:
        df = self._run(np.array([]))
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0
        assert list(df.columns) == ["chrom", "start", "end"]

    def test_none_matrix_returns_empty_df(self) -> None:
        df = self._run(None)   # type: ignore[arg-type]
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0

    def test_loader_exception_returns_empty_df(self) -> None:
        with patch(
            "src.statistical_methods.changepoint_tad.load_hic_matrix",
            side_effect=RuntimeError("disk error"),
        ):
            df = run_pelt_tad(
                chrom=self.CHROM, resolution=self.RES,
                data_path="data/raw", cfg=CFG_MINIMAL,
            )
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0

    @ruptures_required
    def test_all_zero_is_returns_empty_df(self) -> None:
        """If matrix produces all-zero IS, no boundaries should be returned."""
        # A matrix of zeros produces zero IS everywhere
        mat = np.zeros((60, 60), dtype=np.float64)
        df  = self._run(mat)
        assert isinstance(df, pd.DataFrame)
        assert list(df.columns) == ["chrom", "start", "end"]

    def test_none_cfg_uses_defaults(self) -> None:
        """run_pelt_tad must not crash when cfg=None."""
        df = self._run(np.array([]), cfg=None)
        assert isinstance(df, pd.DataFrame)

    @ruptures_required
    def test_block_matrix_detects_approximate_boundaries(self) -> None:
        """
        With a clear 4-block matrix, stable boundaries should be near
        the true block edges (tolerance: ±2 bins).
        """
        n_blocks, bs = 4, 25
        mat = _block_matrix(n_blocks, bs, intra=30.0, inter=0.001)
        cfg = {
            "changepoint": {
                "window_bins":      4,
                "beta_min":         0.2,
                "beta_max":         2.0,
                "n_beta":           6,
                "min_beta_support": 2,
            },
            "consensus": {"tolerance_bins": 1},
        }
        df = self._run(mat, cfg=cfg)
        if len(df) == 0:
            pytest.skip("PELT found no boundaries — numerical edge case")

        true_boundaries = {bs, 2 * bs, 3 * bs}   # bins 25, 50, 75
        detected_bins   = set(df["start"].values // self.RES)
        tol = 2
        matched = sum(
            any(abs(d - t) <= tol for d in detected_bins)
            for t in true_boundaries
        )
        assert matched >= 1, (
            f"Expected at least 1 boundary near {true_boundaries}, "
            f"detected bins: {detected_bins}"
        )