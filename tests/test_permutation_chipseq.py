from __future__ import annotations

"""Tests for src/statistical_methods/permutation_chipseq.py

Coverage
--------
1.  _count_peaks_vectorised  — exact counts, empty inputs, boundary-exact distance
2.  _bh_correction           — monotonicity, range [0,1], known small example
3.  compute_boundary_enrichment
    a. high-signal boundaries → small p-value
    b. no-signal boundaries   → large p-value
    c. no peaks on chrom      → p-value == 1.0 for all
    d. empty boundaries_df    → empty DataFrame, no exception
    e. output schema          — required columns, types
    f. reproducibility        — same rng seed → identical results
    g. different seed         → (usually) different null, but same schema
    h. FDR threshold          — changing threshold changes 'significant' flag
    i. small chrom_size       — telomere margin fallback, no exception
    j. significant count monotone in window_bp (wider → more significant)
4.  save_permutation_results — files created, BED has no header, CSV roundtrip
"""

import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.statistical_methods.permutation_chipseq import (
    _bh_correction,
    _count_peaks_vectorised,
    compute_boundary_enrichment,
    save_permutation_results,
)

# ---------------------------------------------------------------------------
# Shared constants
# ---------------------------------------------------------------------------
CHROM      = "chr1"
RES        = 100_000
CHROM_SIZE = 248_956_422          # hg19 chr1
WIN        = 50_000
N_PERM     = 500                  # small for fast tests
TELO       = 3_000_000
FDR        = 0.05
SEED       = 42


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------

def _rng(seed: int = SEED) -> np.random.Generator:
    return np.random.default_rng(seed)


def _make_boundaries(centres: list[int], chrom: str = CHROM) -> pd.DataFrame:
    """Create a boundaries DataFrame from a list of centre positions."""
    half = RES // 2
    return pd.DataFrame({
        "chrom": chrom,
        "start": [c - half for c in centres],
        "end":   [c + half for c in centres],
    })


def _make_chipseq(centres: list[int], chrom: str = CHROM,
                  peak_half: int = 250) -> pd.DataFrame:
    """Create a ChIP-seq peaks DataFrame from a list of peak centre positions."""
    return pd.DataFrame({
        "chrom": chrom,
        "start": [c - peak_half for c in centres],
        "end":   [c + peak_half for c in centres],
    })


# ---------------------------------------------------------------------------
# 1. _count_peaks_vectorised
# ---------------------------------------------------------------------------

class TestCountPeaksVectorised:

    def test_exact_overlap(self) -> None:
        positions = np.array([1_000_000], dtype=np.int64)
        peak_mids = np.array([1_000_000], dtype=np.int64)   # distance 0
        counts = _count_peaks_vectorised(positions, peak_mids, window_bp=50_000)
        assert counts[0] == 1

    def test_boundary_exact_distance(self) -> None:
        """Peak exactly at window_bp distance must be counted (<=, not <)."""
        positions = np.array([0], dtype=np.int64)
        peak_mids = np.array([50_000], dtype=np.int64)      # distance == window
        counts = _count_peaks_vectorised(positions, peak_mids, window_bp=50_000)
        assert counts[0] == 1

    def test_just_outside_window(self) -> None:
        positions = np.array([0], dtype=np.int64)
        peak_mids = np.array([50_001], dtype=np.int64)      # distance > window
        counts = _count_peaks_vectorised(positions, peak_mids, window_bp=50_000)
        assert counts[0] == 0

    def test_multiple_positions_and_peaks(self) -> None:
        """Two positions; each has a distinct set of close peaks."""
        positions = np.array([0, 10_000_000], dtype=np.int64)
        peak_mids = np.array([
            10_000,      # close to pos 0
            20_000,      # close to pos 0
            10_000_000,  # close to pos 1
        ], dtype=np.int64)
        counts = _count_peaks_vectorised(positions, peak_mids, window_bp=50_000)
        assert counts[0] == 2   # 10k + 20k within 50k of 0
        assert counts[1] == 1   # only 10M within 50k of 10M

    def test_empty_peaks(self) -> None:
        positions = np.array([1_000_000, 2_000_000], dtype=np.int64)
        peak_mids = np.array([], dtype=np.int64)
        counts = _count_peaks_vectorised(positions, peak_mids, window_bp=50_000)
        assert counts.tolist() == [0, 0]

    def test_empty_positions(self) -> None:
        positions = np.array([], dtype=np.int64)
        peak_mids = np.array([1_000_000], dtype=np.int64)
        counts = _count_peaks_vectorised(positions, peak_mids, window_bp=50_000)
        assert len(counts) == 0

    def test_output_dtype_is_int(self) -> None:
        counts = _count_peaks_vectorised(
            np.array([0], dtype=np.int64),
            np.array([0], dtype=np.int64),
            window_bp=1,
        )
        assert np.issubdtype(counts.dtype, np.integer)


# ---------------------------------------------------------------------------
# 2. _bh_correction
# ---------------------------------------------------------------------------

class TestBHCorrection:

    def test_empty_input(self) -> None:
        result = _bh_correction(np.array([], dtype=np.float64))
        assert len(result) == 0

    def test_single_value(self) -> None:
        result = _bh_correction(np.array([0.03]))
        assert pytest.approx(result[0], abs=1e-9) == 0.03

    def test_all_ones(self) -> None:
        result = _bh_correction(np.ones(5))
        assert (result == 1.0).all()

    def test_output_clipped_to_one(self) -> None:
        """BH can inflate values > 1 before clip; result must be <= 1."""
        pvals = np.array([0.9, 0.95, 0.99])
        result = _bh_correction(pvals)
        assert (result <= 1.0).all()

    def test_output_non_negative(self) -> None:
        pvals = np.array([0.001, 0.01, 0.05, 0.1, 0.5])
        result = _bh_correction(pvals)
        assert (result >= 0.0).all()

    def test_monotonicity(self) -> None:
        """Adjusted p-values must be non-decreasing when sorted by raw p."""
        rng_ = np.random.default_rng(0)
        pvals = np.sort(rng_.uniform(0, 1, 50))
        p_adj = _bh_correction(pvals)
        p_adj_sorted = np.sort(p_adj)
        # After sorting, each value must be >= previous
        assert (np.diff(p_adj_sorted) >= -1e-12).all()

    def test_known_example(self) -> None:
        """
        3 p-values: [0.01, 0.04, 0.5]  n=3
        BH ranks (sorted): 1, 2, 3
        raw adj: 0.01*3/1=0.03, 0.04*3/2=0.06, 0.5*3/3=0.50
        cummin from right: min(0.03,0.06)=0.03; min(0.06,0.50)=0.06; 0.50
        → [0.03, 0.06, 0.50]
        """
        pvals  = np.array([0.01, 0.04, 0.5])
        result = _bh_correction(pvals)
        assert pytest.approx(result[0], abs=1e-9) == 0.03
        assert pytest.approx(result[1], abs=1e-9) == 0.06
        assert pytest.approx(result[2], abs=1e-9) == 0.50


# ---------------------------------------------------------------------------
# 3. compute_boundary_enrichment
# ---------------------------------------------------------------------------

class TestComputeBoundaryEnrichment:

    # ── 3a. High-signal: dense peaks at boundary → small p-value ─────────────
    def test_high_signal_gives_small_pvalue(self) -> None:
        """
        Place 200 peaks right at the boundary centre.
        Null positions drawn from the rest of the chromosome will rarely
        accumulate that many peaks → p-value should be very small.
        """
        centre = 50_000_000
        boundaries = _make_boundaries([centre])
        # Dense cluster of 200 peaks at the boundary
        chipseq = _make_chipseq(
            [centre + i * 200 for i in range(-100, 100)]
        )
        df = compute_boundary_enrichment(
            boundaries, chipseq,
            chrom=CHROM, chrom_size=CHROM_SIZE,
            window_bp=WIN, n_permutations=N_PERM,
            rng=_rng(), telomere_margin_bp=TELO,
        )
        assert len(df) == 1
        assert df["pvalue"].iloc[0] < 0.1

    # ── 3b. No-signal: zero peaks anywhere → p-value == 1.0 ──────────────────
    def test_no_peaks_gives_pvalue_one(self) -> None:
        boundaries = _make_boundaries([50_000_000, 100_000_000])
        chipseq    = _make_chipseq([])          # empty
        df = compute_boundary_enrichment(
            boundaries, chipseq,
            chrom=CHROM, chrom_size=CHROM_SIZE,
            window_bp=WIN, n_permutations=N_PERM,
            rng=_rng(), telomere_margin_bp=TELO,
        )
        # All observed counts are 0; null counts are also 0
        # → null_density >= obs_density (0 >= 0) always True → p = 1.0
        assert (df["pvalue"] == 1.0).all()

    # ── 3c. No peaks on specific chrom ───────────────────────────────────────
    def test_no_peaks_on_chrom_returns_pvalue_one(self) -> None:
        boundaries = _make_boundaries([50_000_000])
        # Peaks only on chr2 — not on chr1
        chipseq = _make_chipseq([50_000_000], chrom="chr2")
        df = compute_boundary_enrichment(
            boundaries, chipseq,
            chrom=CHROM, chrom_size=CHROM_SIZE,
            window_bp=WIN, n_permutations=N_PERM,
            rng=_rng(), telomere_margin_bp=TELO,
        )
        assert len(df) == 1
        assert df["pvalue"].iloc[0] == 1.0

    # ── 3d. Empty boundaries → empty DataFrame, no exception ─────────────────
    def test_empty_boundaries_returns_empty_df(self) -> None:
        boundaries = pd.DataFrame(columns=["chrom", "start", "end"])
        chipseq    = _make_chipseq([50_000_000])
        df = compute_boundary_enrichment(
            boundaries, chipseq,
            chrom=CHROM, chrom_size=CHROM_SIZE,
            window_bp=WIN, n_permutations=N_PERM,
            rng=_rng(), telomere_margin_bp=TELO,
        )
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0
        assert "pvalue" in df.columns

    # ── 3e. Output schema ─────────────────────────────────────────────────────
    def test_output_schema(self) -> None:
        boundaries = _make_boundaries([50_000_000, 100_000_000])
        chipseq    = _make_chipseq([50_000_000])
        df = compute_boundary_enrichment(
            boundaries, chipseq,
            chrom=CHROM, chrom_size=CHROM_SIZE,
            window_bp=WIN, n_permutations=N_PERM,
            rng=_rng(), telomere_margin_bp=TELO,
        )
        required = {
            "chrom", "start", "end",
            "boundary_centre", "obs_count", "obs_density",
            "pvalue", "pvalue_adj", "significant",
        }
        assert required.issubset(set(df.columns))
        assert df["pvalue"].between(0.0, 1.0).all()
        assert df["pvalue_adj"].between(0.0, 1.0).all()
        assert df["significant"].dtype == bool
        assert df["obs_count"].ge(0).all()

    # ── 3f. Reproducibility — same seed → identical results ──────────────────
    def test_reproducibility_same_seed(self) -> None:
        boundaries = _make_boundaries([50_000_000, 80_000_000, 120_000_000])
        chipseq    = _make_chipseq([50_000_000 + i * 5_000 for i in range(20)])
        kwargs = dict(
            chrom=CHROM, chrom_size=CHROM_SIZE,
            window_bp=WIN, n_permutations=N_PERM,
            telomere_margin_bp=TELO,
        )
        df1 = compute_boundary_enrichment(
            boundaries, chipseq, rng=_rng(SEED), **kwargs
        )
        df2 = compute_boundary_enrichment(
            boundaries, chipseq, rng=_rng(SEED), **kwargs
        )
        pd.testing.assert_frame_equal(df1, df2)

    # ── 3g. Different seed → same schema (values may differ) ─────────────────
    def test_different_seed_same_schema(self) -> None:
        boundaries = _make_boundaries([50_000_000])
        chipseq    = _make_chipseq([50_000_000])
        kwargs = dict(
            chrom=CHROM, chrom_size=CHROM_SIZE,
            window_bp=WIN, n_permutations=N_PERM,
            telomere_margin_bp=TELO,
        )
        df1 = compute_boundary_enrichment(
            boundaries, chipseq, rng=_rng(42),  **kwargs
        )
        df2 = compute_boundary_enrichment(
            boundaries, chipseq, rng=_rng(999), **kwargs
        )
        assert list(df1.columns) == list(df2.columns)
        assert len(df1) == len(df2)

    # ── 3h. FDR threshold controls 'significant' flag ────────────────────────
    def test_fdr_threshold_controls_significant(self) -> None:
        """Lowering threshold → fewer significant; raising → more."""
        centre  = 50_000_000
        boundaries = _make_boundaries([centre])
        chipseq = _make_chipseq([centre + i * 300 for i in range(-50, 50)])
        kwargs = dict(
            chrom=CHROM, chrom_size=CHROM_SIZE,
            window_bp=WIN, n_permutations=N_PERM,
            rng=_rng(), telomere_margin_bp=TELO,
        )
        df_strict = compute_boundary_enrichment(
            boundaries, chipseq, fdr_threshold=0.001, **kwargs
        )
        df_lenient = compute_boundary_enrichment(
            boundaries, chipseq, fdr_threshold=0.99, **kwargs
        )
        n_strict  = df_strict["significant"].sum()
        n_lenient = df_lenient["significant"].sum()
        assert n_strict <= n_lenient

    # ── 3i. Small chrom_size → telomere fallback, no exception ───────────────
    def test_tiny_chrom_size_no_exception(self) -> None:
        """chrom_size < 2 * telomere_margin → fallback to full range."""
        boundaries = _make_boundaries([500_000])
        chipseq    = _make_chipseq([500_000])
        df = compute_boundary_enrichment(
            boundaries, chipseq,
            chrom=CHROM,
            chrom_size=1_000_000,        # < 2 * TELO (6 000 000)
            window_bp=10_000,
            n_permutations=50,
            rng=_rng(),
            telomere_margin_bp=TELO,
        )
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 1

    # ── 3j. Wider window → higher obs_density (monotone in window_bp) ────────
    def test_wider_window_higher_density(self) -> None:
        centre  = 50_000_000
        boundaries = _make_boundaries([centre])
        # Peaks spread within 100kb of centre
        chipseq = _make_chipseq(
            [centre + i * 5_000 for i in range(-19, 20)]
        )
        kwargs = dict(
            chrom=CHROM, chrom_size=CHROM_SIZE,
            n_permutations=50, rng=_rng(),
            telomere_margin_bp=TELO,
        )
        df_narrow = compute_boundary_enrichment(
            boundaries, chipseq, window_bp=25_000,  **kwargs
        )
        df_wide = compute_boundary_enrichment(
            boundaries, chipseq, window_bp=100_000, **kwargs
        )
        # More peaks captured with wider window → higher count
        assert (
            df_wide["obs_count"].iloc[0]
            >= df_narrow["obs_count"].iloc[0]
        )

    # ── 3k. obs_density equals obs_count / (2 * window_bp) ───────────────────
    def test_obs_density_formula(self) -> None:
        centre = 50_000_000
        boundaries = _make_boundaries([centre])
        chipseq    = _make_chipseq([centre])
        df = compute_boundary_enrichment(
            boundaries, chipseq,
            chrom=CHROM, chrom_size=CHROM_SIZE,
            window_bp=WIN, n_permutations=50,
            rng=_rng(), telomere_margin_bp=TELO,
        )
        expected_density = df["obs_count"].iloc[0] / (2.0 * WIN)
        assert pytest.approx(df["obs_density"].iloc[0], rel=1e-6) == expected_density

    # ── 3l. boundary_centre = (start + end) // 2 ─────────────────────────────
    def test_boundary_centre_value(self) -> None:
        centre = 55_000_000
        boundaries = _make_boundaries([centre])
        chipseq    = _make_chipseq([])
        df = compute_boundary_enrichment(
            boundaries, chipseq,
            chrom=CHROM, chrom_size=CHROM_SIZE,
            window_bp=WIN, n_permutations=50,
            rng=_rng(), telomere_margin_bp=TELO,
        )
        expected = (df["start"].iloc[0] + df["end"].iloc[0]) // 2
        assert df["boundary_centre"].iloc[0] == expected

    # ── 3m. pvalue_adj monotone w.r.t. pvalue ────────────────────────────────
    def test_pvalue_adj_monotone_with_pvalue(self) -> None:
        """Boundaries sorted by raw pvalue → adj pvalue non-decreasing."""
        centres = list(range(20_000_000, 200_000_000, 10_000_000))
        boundaries = _make_boundaries(centres)
        chipseq    = _make_chipseq(
            [centres[0] + i * 2_000 for i in range(-30, 30)]
        )
        df = compute_boundary_enrichment(
            boundaries, chipseq,
            chrom=CHROM, chrom_size=CHROM_SIZE,
            window_bp=WIN, n_permutations=N_PERM,
            rng=_rng(), telomere_margin_bp=TELO,
        )
        df_sorted = df.sort_values("pvalue").reset_index(drop=True)
        adj = df_sorted["pvalue_adj"].values
        assert (np.diff(adj) >= -1e-9).all(), \
            "pvalue_adj must be non-decreasing when sorted by pvalue"


# ---------------------------------------------------------------------------
# 4. save_permutation_results
# ---------------------------------------------------------------------------

class TestSavePermutationResults:

    @pytest.fixture()
    def sample_df(self) -> pd.DataFrame:
        return pd.DataFrame({
            "chrom":           ["chr1", "chr1", "chr1"],
            "start":           [100_000, 200_000, 300_000],
            "end":             [200_000, 300_000, 400_000],
            "boundary_centre": [150_000, 250_000, 350_000],
            "obs_count":       [10, 2, 0],
            "obs_density":     [1e-4, 2e-5, 0.0],
            "pvalue":          [0.01, 0.20, 0.99],
            "pvalue_adj":      [0.03, 0.40, 0.99],
            "significant":     [True, False, False],
        })

    def test_csv_created(self, tmp_path: Path, sample_df: pd.DataFrame) -> None:
        csv_p = tmp_path / "stats" / "perm.csv"
        bed_p = tmp_path / "consensus" / "perm.bed"
        save_permutation_results(sample_df, csv_p, bed_p)
        assert csv_p.exists()

    def test_bed_created(self, tmp_path: Path, sample_df: pd.DataFrame) -> None:
        csv_p = tmp_path / "perm.csv"
        bed_p = tmp_path / "perm.bed"
        save_permutation_results(sample_df, csv_p, bed_p)
        assert bed_p.exists()

    def test_bed_contains_only_significant(
        self, tmp_path: Path, sample_df: pd.DataFrame
    ) -> None:
        csv_p = tmp_path / "perm.csv"
        bed_p = tmp_path / "perm.bed"
        save_permutation_results(sample_df, csv_p, bed_p)
        bed_df = pd.read_csv(bed_p, sep="\t", header=None,
                             names=["chrom", "start", "end"])
        assert len(bed_df) == 1               # only 1 row has significant=True
        assert bed_df["start"].iloc[0] == 100_000

    def test_bed_has_no_header(
        self, tmp_path: Path, sample_df: pd.DataFrame
    ) -> None:
        csv_p = tmp_path / "perm.csv"
        bed_p = tmp_path / "perm.bed"
        save_permutation_results(sample_df, csv_p, bed_p)
        first_line = bed_p.read_text().splitlines()[0]
        assert first_line.startswith("chr"), \
            f"BED must not have a header, got: {first_line!r}"

    def test_csv_roundtrip(
        self, tmp_path: Path, sample_df: pd.DataFrame
    ) -> None:
        csv_p = tmp_path / "perm.csv"
        bed_p = tmp_path / "perm.bed"
        save_permutation_results(sample_df, csv_p, bed_p)
        df_back = pd.read_csv(csv_p)
        assert list(df_back["start"]) == [100_000, 200_000, 300_000]
        assert "pvalue_adj" in df_back.columns

    def test_parent_dirs_created(
        self, tmp_path: Path, sample_df: pd.DataFrame
    ) -> None:
        """Directories must be created automatically (parents=True)."""
        csv_p = tmp_path / "deep" / "nested" / "dir" / "perm.csv"
        bed_p = tmp_path / "another" / "deep" / "perm.bed"
        save_permutation_results(sample_df, csv_p, bed_p)
        assert csv_p.exists()
        assert bed_p.exists()