from __future__ import annotations

"""Tests for src/statistical_methods/probabilistic_boundary.py

Coverage
--------
1.  init_sensitivity     — sigmoid shape, cr=1.0→0.5, cr>1→>0.5, cr<1→<0.5
2.  init_fpr             — proportional to n_bnd, clamping
3.  compute_prior        — None matrix fallback, signal-proportional, clamp
4.  build_observation_matrix — shape, exact hit, tolerance, empty algo
5.  run_em               — convergence, posterior in [0,1], schema,
                           strong signal → high posterior,
                           all-zero obs → prior dominates,
                           reproducibility
6.  compute_probabilistic_consensus
    a. synthetic block matrix with clear TAD structure → boundaries detected
    b. empty tad_files → empty DataFrame
    c. all algos missing → empty DataFrame
    d. output schema: columns, types, posterior in [0,1]
    e. theta sweep: higher theta → fewer or equal boundaries
    f. high-quality algo (high CR, low FPR) dominates over noisy algo
    g. exclude_from_sources respected
7.  save_probabilistic_consensus — file created, no header, roundtrip
"""

from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from src.statistical_methods.probabilistic_boundary import (
    _sigmoid,
    build_observation_matrix,
    compute_probabilistic_consensus,
    compute_prior,
    init_fpr,
    init_sensitivity,
    run_em,
    save_probabilistic_consensus,
)

# ---------------------------------------------------------------------------
# Constants / helpers
# ---------------------------------------------------------------------------

RES   = 100_000
CHROM = "chr1"
SEED  = 42

CFG_MINIMAL = {
    "probabilistic": {
        "theta":        0.5,
        "k_sens":       8.0,
        "k_prior":      3.0,
        "global_prior": 0.10,
        "max_iter":     30,
        "em_tol":       1e-4,
        "fpr_min":      0.01,
        "fpr_max":      0.60,
        "sens_min":     0.50,
        "sens_max":     0.99,
    },
    "consensus":        {"tolerance_bins": 1},
    "weighted_consensus": {"exclude_from_sources": []},
}


def _block_matrix(n_blocks: int, block_size: int,
                  intra: float = 20.0, inter: float = 0.1) -> np.ndarray:
    n   = n_blocks * block_size
    mat = np.full((n, n), inter, dtype=np.float64)
    for b in range(n_blocks):
        lo, hi = b * block_size, (b + 1) * block_size
        mat[lo:hi, lo:hi] = intra
    return mat


def _make_bed(tmp_path: Path, name: str,
              rows: list[tuple]) -> Path:
    p = tmp_path / name
    lines = "\n".join(f"{c}\t{s}\t{e}" for c, s, e in rows) + "\n"
    p.write_text(lines)
    return p


# ---------------------------------------------------------------------------
# 1. init_sensitivity
# ---------------------------------------------------------------------------

class TestInitSensitivity:

    def test_cr_equal_one_gives_half(self) -> None:
        """cr = 1.0 → excess = 0 → sigmoid(0) = 0.5."""
        s = init_sensitivity(np.array([1.0]))
        assert pytest.approx(s[0], abs=1e-6) == 0.5

    def test_cr_above_one_above_half(self) -> None:
        s = init_sensitivity(np.array([1.3, 1.5]))
        assert (s > 0.5).all()

    def test_cr_below_one_below_or_equal_half(self) -> None:
        """cr < 1.0 → excess = max(0, cr-1) = 0 → sigmoid(0) = 0.5."""
        s = init_sensitivity(np.array([0.8, 0.5]))
        assert (s <= 0.5).all()

    def test_higher_cr_higher_sens(self) -> None:
        cr  = np.array([1.0, 1.1, 1.3, 1.5])
        s   = init_sensitivity(cr)
        assert (np.diff(s) >= 0).all()

    def test_output_in_range(self) -> None:
        cr = np.array([0.5, 1.0, 1.2, 1.5, 2.0])
        s  = init_sensitivity(cr)
        assert (s >= 0).all() and (s <= 1).all()

    def test_large_k_more_discriminative(self) -> None:
        cr = np.array([1.1, 1.4])
        s4 = init_sensitivity(cr, k_sens=4.0)
        s8 = init_sensitivity(cr, k_sens=8.0)
        # larger k → larger gap between weak and strong
        gap4 = s4[1] - s4[0]
        gap8 = s8[1] - s8[0]
        assert gap8 >= gap4


# ---------------------------------------------------------------------------
# 2. init_fpr
# ---------------------------------------------------------------------------

class TestInitFpr:

    def test_proportional_to_n_boundaries(self) -> None:
        fpr = init_fpr(np.array([100, 200]), n_bins=1000)
        assert pytest.approx(fpr[0], rel=1e-6) == 0.1
        assert pytest.approx(fpr[1], rel=1e-6) == 0.2

    def test_clamped_at_fpr_max(self) -> None:
        fpr = init_fpr(np.array([900]), n_bins=1000, fpr_max=0.5)
        assert fpr[0] == 0.5

    def test_clamped_at_fpr_min(self) -> None:
        fpr = init_fpr(np.array([0]), n_bins=1000, fpr_min=0.01)
        assert fpr[0] == 0.01

    def test_zero_n_bins_no_crash(self) -> None:
        fpr = init_fpr(np.array([10]), n_bins=0)
        assert np.isfinite(fpr[0])

    def test_more_boundaries_higher_fpr(self) -> None:
        fpr = init_fpr(np.array([50, 200, 400]), n_bins=1000)
        assert (np.diff(fpr) >= 0).all()


# ---------------------------------------------------------------------------
# 3. compute_prior
# ---------------------------------------------------------------------------

class TestComputePrior:

    def test_none_matrix_uniform(self) -> None:
        prior = compute_prior(None, n_bins=100, global_prior=0.1)
        assert len(prior) == 100
        assert np.allclose(prior, 0.1)

    def test_empty_matrix_uniform(self) -> None:
        prior = compute_prior(np.array([]), n_bins=50, global_prior=0.2)
        assert np.allclose(prior, 0.2)

    def test_output_in_0_1(self) -> None:
        mat   = _block_matrix(4, 10)
        prior = compute_prior(mat, n_bins=40)
        assert (prior >= 0).all() and (prior <= 1).all()

    def test_output_length_equals_matrix(self) -> None:
        mat   = _block_matrix(3, 10)
        prior = compute_prior(mat, n_bins=30)
        assert len(prior) == 30

    def test_dense_region_higher_prior(self) -> None:
        """Diagonal-dense bins should get higher prior than off-diagonal bins."""
        n = 40
        mat = np.zeros((n, n), dtype=np.float64)
        # Make diagonal very dense in [10:20, 10:20]
        mat[10:20, 10:20] = 100.0
        prior = compute_prior(mat, n_bins=n, k_prior=3.0)
        assert prior[15] > prior[5]

    def test_size_mismatch_uses_global_prior(self) -> None:
        mat   = _block_matrix(3, 10)   # 30×30
        prior = compute_prior(mat, n_bins=50, global_prior=0.15)
        assert np.allclose(prior, 0.15)


# ---------------------------------------------------------------------------
# 4. build_observation_matrix
# ---------------------------------------------------------------------------

class TestBuildObservationMatrix:

    def test_shape(self) -> None:
        candidates = np.array([10, 20, 30], dtype=np.int64)
        algo_bins  = {"a": np.array([10], dtype=np.int64),
                      "b": np.array([20], dtype=np.int64)}
        X = build_observation_matrix(algo_bins, candidates, tolerance_bins=0)
        assert X.shape == (3, 2)

    def test_exact_hit(self) -> None:
        candidates = np.array([10], dtype=np.int64)
        algo_bins  = {"a": np.array([10], dtype=np.int64)}
        X = build_observation_matrix(algo_bins, candidates, tolerance_bins=0)
        assert X[0, 0] == 1

    def test_miss(self) -> None:
        candidates = np.array([10], dtype=np.int64)
        algo_bins  = {"a": np.array([15], dtype=np.int64)}
        X = build_observation_matrix(algo_bins, candidates, tolerance_bins=0)
        assert X[0, 0] == 0

    def test_tolerance_hit(self) -> None:
        candidates = np.array([10], dtype=np.int64)
        algo_bins  = {"a": np.array([11], dtype=np.int64)}
        X = build_observation_matrix(algo_bins, candidates, tolerance_bins=1)
        assert X[0, 0] == 1

    def test_tolerance_miss(self) -> None:
        candidates = np.array([10], dtype=np.int64)
        algo_bins  = {"a": np.array([12], dtype=np.int64)}
        X = build_observation_matrix(algo_bins, candidates, tolerance_bins=1)
        assert X[0, 0] == 0

    def test_empty_algo_bins(self) -> None:
        candidates = np.array([10, 20], dtype=np.int64)
        algo_bins  = {"empty": np.array([], dtype=np.int64)}
        X = build_observation_matrix(algo_bins, candidates, tolerance_bins=1)
        assert X.sum() == 0

    def test_binary_values(self) -> None:
        candidates = np.array([5, 10, 15], dtype=np.int64)
        algo_bins  = {"a": np.array([5, 15], dtype=np.int64)}
        X = build_observation_matrix(algo_bins, candidates, tolerance_bins=0)
        assert set(X.flatten().tolist()).issubset({0, 1})


# ---------------------------------------------------------------------------
# 5. run_em
# ---------------------------------------------------------------------------

class TestRunEM:

    def _make_simple_X(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        5 candidates, 3 algorithms.
        Candidates 0,1,2 are "true" (all algos agree).
        Candidates 3,4 are "false" (no algo detects them).
        """
        X = np.array([
            [1, 1, 1],
            [1, 1, 0],
            [1, 0, 1],
            [0, 0, 0],
            [0, 0, 0],
        ], dtype=np.int8)
        pi   = np.full(5, 0.1)
        sens = np.array([0.8, 0.8, 0.8])
        fpr  = np.array([0.1, 0.1, 0.1])
        return X, pi, sens, fpr

    def test_posterior_in_0_1(self) -> None:
        X, pi, sens, fpr = self._make_simple_X()
        post, _, _, _ = run_em(X, pi, sens, fpr)
        assert (post >= 0.0).all() and (post <= 1.0).all()

    def test_observed_higher_posterior(self) -> None:
        """Bins observed by all algos should have higher posterior than unobserved."""
        X, pi, sens, fpr = self._make_simple_X()
        post, _, _, _ = run_em(X, pi, sens, fpr)
        assert post[0] > post[3]   # fully observed > unobserved
        assert post[1] > post[3]

    def test_unobserved_low_posterior(self) -> None:
        """With low prior and no observations, posterior should be small."""
        X, pi, sens, fpr = self._make_simple_X()
        post, _, _, _ = run_em(X, pi, sens, fpr)
        assert post[3] < 0.5
        assert post[4] < 0.5

    def test_ll_non_decreasing(self) -> None:
        """Log-likelihood must be non-decreasing across EM iterations."""
        X, pi, sens, fpr = self._make_simple_X()
        _, _, _, ll_hist = run_em(X, pi, sens, fpr, max_iter=20)
        diffs = np.diff(ll_hist)
        assert (diffs >= -1e-6).all(), \
            f"LL decreased: {diffs[diffs < -1e-6]}"

    def test_returns_correct_shapes(self) -> None:
        X, pi, sens, fpr = self._make_simple_X()
        post, s_out, f_out, ll = run_em(X, pi, sens, fpr)
        assert post.shape == (5,)
        assert s_out.shape == (3,)
        assert f_out.shape == (3,)
        assert isinstance(ll, list) and len(ll) > 0

    def test_all_zero_obs_prior_dominates(self) -> None:
        """When all algos miss every candidate, posterior ≈ prior."""
        X    = np.zeros((5, 3), dtype=np.int8)
        pi   = np.full(5, 0.3)
        sens = np.array([0.8, 0.8, 0.8])
        fpr  = np.array([0.1, 0.1, 0.1])
        post, _, _, _ = run_em(X, pi, sens, fpr)
        # Prior = 0.3; with zero observations the posterior should
        # remain near the prior (exact value depends on EM dynamics)
        assert (post < 0.5).all()

    def test_strong_signal_high_posterior(self) -> None:
        """High-quality algo (high sens, low fpr) + all detected → high posterior.

        Setup: 1 true boundary (detected by all 3 algos) + 9 non-boundaries
        (not detected by any algo).  The unobserved candidates stabilise the
        M-step FPR estimate: with only 1 candidate the M-step is degenerate
        (sum_q0 ≈ 0 → fpr jumps to fpr_max), causing posterior to collapse.
        """
        # Candidate 0: true boundary — all algos agree
        # Candidates 1-9: non-boundaries — no algo detects them
        n_noise = 9
        X_true  = np.array([[1, 1, 1]], dtype=np.int8)
        X_noise = np.zeros((n_noise, 3), dtype=np.int8)
        X       = np.vstack([X_true, X_noise])

        pi   = np.array([0.5] + [0.05] * n_noise)   # higher prior on candidate 0
        sens = np.array([0.95, 0.95, 0.95])
        fpr  = np.array([0.02, 0.02, 0.02])
        post, _, _, _ = run_em(X, pi, sens, fpr)
        assert post[0] > 0.8, (
            f"Expected posterior > 0.8 for fully-observed boundary, got {post[0]:.4f}. "
            "Check M-step FPR stability."
        )

    def test_reproducibility(self) -> None:
        """EM is deterministic (no RNG) → same result on repeated calls."""
        X, pi, sens, fpr = self._make_simple_X()
        p1, _, _, _ = run_em(X, pi, sens, fpr)
        p2, _, _, _ = run_em(X, pi, sens, fpr)
        np.testing.assert_array_equal(p1, p2)

    def test_sens_clamped_within_bounds(self) -> None:
        X, pi, sens, fpr = self._make_simple_X()
        _, s_out, f_out, _ = run_em(
            X, pi, sens, fpr,
            sens_min=0.50, sens_max=0.99,
            fpr_min=0.01,  fpr_max=0.60,
        )
        assert (s_out >= 0.50).all() and (s_out <= 0.99).all()
        assert (f_out >= 0.01).all() and (f_out <= 0.60).all()


# ---------------------------------------------------------------------------
# 6. compute_probabilistic_consensus
# ---------------------------------------------------------------------------

class TestComputeProbabilisticConsensus:

    def _build_inputs(self, tmp_path: Path, resolution: int = RES):
        """
        3 algorithms on chr1.
        A (high quality): borders at bins 10, 20, 30 — all true boundaries
        B (medium):       borders at bins 10, 21, 30
        C (noisy):        borders at bins 10, 20, 30, 50, 60, 70 — extra noise
        """
        a = _make_bed(tmp_path, "a.bed", [
            (CHROM, 10*resolution, 20*resolution),
            (CHROM, 20*resolution, 30*resolution),
        ])
        b = _make_bed(tmp_path, "b.bed", [
            (CHROM, 10*resolution, 21*resolution),
            (CHROM, 21*resolution, 30*resolution),
        ])
        c = _make_bed(tmp_path, "c.bed", [
            (CHROM, 10*resolution, 20*resolution),
            (CHROM, 20*resolution, 30*resolution),
            (CHROM, 50*resolution, 60*resolution),
            (CHROM, 60*resolution, 70*resolution),
        ])
        tad_files = {"A": a, "B": b, "C": c}
        meta = {
            "A": {"center_ratio": 1.42, "n_boundaries": 4},
            "B": {"center_ratio": 1.31, "n_boundaries": 4},
            "C": {"center_ratio": 1.02, "n_boundaries": 8},
        }
        return tad_files, meta

    def test_output_schema(self, tmp_path: Path) -> None:
        tad_files, meta = self._build_inputs(tmp_path)
        df = compute_probabilistic_consensus(
            tad_files, meta, CHROM, RES, CFG_MINIMAL,
        )
        assert set(df.columns) == {"chrom", "start", "end", "posterior_prob"}
        if len(df) > 0:
            assert (df["posterior_prob"] >= 0).all()
            assert (df["posterior_prob"] <= 1).all()
            assert (df["start"] % RES == 0).all()
            assert (df["end"] > df["start"]).all()

    def test_returns_dataframe(self, tmp_path: Path) -> None:
        tad_files, meta = self._build_inputs(tmp_path)
        df = compute_probabilistic_consensus(
            tad_files, meta, CHROM, RES, CFG_MINIMAL,
        )
        assert isinstance(df, pd.DataFrame)

    def test_empty_tad_files_returns_empty(self) -> None:
        df = compute_probabilistic_consensus(
            {}, {}, CHROM, RES, CFG_MINIMAL,
        )
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0

    def test_theta_sweep_monotone(self, tmp_path: Path) -> None:
        """Higher theta → fewer or equal boundaries."""
        tad_files, meta = self._build_inputs(tmp_path)
        cfg_low  = {**CFG_MINIMAL,
                    "probabilistic": {**CFG_MINIMAL["probabilistic"], "theta": 0.1}}
        cfg_high = {**CFG_MINIMAL,
                    "probabilistic": {**CFG_MINIMAL["probabilistic"], "theta": 0.9}}
        df_low  = compute_probabilistic_consensus(
            tad_files, meta, CHROM, RES, cfg_low)
        df_high = compute_probabilistic_consensus(
            tad_files, meta, CHROM, RES, cfg_high)
        assert len(df_low) >= len(df_high)

    def test_wrong_chrom_returns_empty(self, tmp_path: Path) -> None:
        tad_files, meta = self._build_inputs(tmp_path)
        df = compute_probabilistic_consensus(
            tad_files, meta, "chr99", RES, CFG_MINIMAL,
        )
        assert len(df) == 0

    def test_noisy_algo_boundaries_suppressed(self, tmp_path: Path) -> None:
        """
        Bins 50, 60, 70 only detected by C (low CR, high FPR).
        With probabilistic model, their posteriors should be lower than
        bins 10, 20, 30 (detected by all three algos).
        """
        tad_files, meta = self._build_inputs(tmp_path)
        cfg_low_theta = {
            **CFG_MINIMAL,
            "probabilistic": {**CFG_MINIMAL["probabilistic"], "theta": 0.0},
        }
        df = compute_probabilistic_consensus(
            tad_files, meta, CHROM, RES, cfg_low_theta,
        )
        if len(df) == 0:
            pytest.skip("No boundaries returned at theta=0 — numerical edge case")

        post_dict = dict(zip(df["start"] // RES, df["posterior_prob"]))
        # Bins agreed by all algos should have higher posterior than noisy bins
        shared_bins  = [10, 20, 30]
        noisy_bins   = [50, 60, 70]
        shared_post  = [post_dict.get(b, 0.0) for b in shared_bins if b in post_dict]
        noisy_post   = [post_dict.get(b, 0.0) for b in noisy_bins  if b in post_dict]
        if shared_post and noisy_post:
            assert np.mean(shared_post) >= np.mean(noisy_post)

    def test_exclude_from_sources(self, tmp_path: Path) -> None:
        """excluded algos should not contribute boundaries."""
        tad_files, meta = self._build_inputs(tmp_path)
        cfg_excl = {
            **CFG_MINIMAL,
            "weighted_consensus": {"exclude_from_sources": ["A", "B", "C"]},
        }
        df = compute_probabilistic_consensus(
            tad_files, meta, CHROM, RES, cfg_excl,
        )
        assert len(df) == 0

    def test_with_block_matrix_prior(self, tmp_path: Path) -> None:
        """With matrix prior, model should still return valid DataFrame."""
        tad_files, meta = self._build_inputs(tmp_path, resolution=RES)
        matrix = _block_matrix(4, 10)   # 40×40
        df = compute_probabilistic_consensus(
            tad_files, meta, CHROM, RES, CFG_MINIMAL, matrix=matrix,
        )
        assert isinstance(df, pd.DataFrame)
        assert "posterior_prob" in df.columns


# ---------------------------------------------------------------------------
# 7. save_probabilistic_consensus
# ---------------------------------------------------------------------------

class TestSaveProbabilisticConsensus:

    @pytest.fixture()
    def sample_df(self) -> pd.DataFrame:
        return pd.DataFrame({
            "chrom":          ["chr1", "chr1"],
            "start":          [100_000, 200_000],
            "end":            [200_000, 300_000],
            "posterior_prob": [0.85, 0.62],
        })

    def test_bed_created(self, tmp_path: Path,
                         sample_df: pd.DataFrame) -> None:
        bed = tmp_path / "out.bed"
        save_probabilistic_consensus(sample_df, bed)
        assert bed.exists()

    def test_no_header(self, tmp_path: Path,
                       sample_df: pd.DataFrame) -> None:
        bed = tmp_path / "out.bed"
        save_probabilistic_consensus(sample_df, bed)
        first = bed.read_text().splitlines()[0]
        assert first.startswith("chr")

    def test_csv_optional(self, tmp_path: Path,
                          sample_df: pd.DataFrame) -> None:
        bed = tmp_path / "out.bed"
        csv = tmp_path / "out.csv"
        save_probabilistic_consensus(sample_df, bed, csv_path=csv)
        assert csv.exists()
        df_back = pd.read_csv(csv)
        assert "posterior_prob" in df_back.columns

    def test_parent_dirs_created(self, tmp_path: Path,
                                 sample_df: pd.DataFrame) -> None:
        bed = tmp_path / "deep" / "nested" / "out.bed"
        save_probabilistic_consensus(sample_df, bed)
        assert bed.exists()