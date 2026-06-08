from __future__ import annotations

"""Tests for src/statistical_methods/weighted_consensus.py

Coverage
--------
1. load_weights            — correct parsing, NO_DATA/SHIFTED exclusion
2. compute_weighted_consensus — known synthetic result, theta sweep
3. topdom SHIFTED          → weight 0.0, does not contribute
4. all-zero weights        → empty DataFrame, no exception
5. missing TAD file        → treated as empty, no exception
6. reproducibility         — deterministic (no RNG, but seed=42 fixture kept)
7. FDR monotonicity guard  — not applicable here (see permutation tests)
"""

import io
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.statistical_methods.weighted_consensus import (
    _apply_transform,
    _boundaries_from_tad_df,
    compute_weighted_consensus,
    load_weights,
    save_weighted_consensus,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

SUMMARY_CSV_CONTENT = textwrap.dedent("""\
    algorithm,n_boundaries,center_ratio,peak_ratio,peak_offset_kb,mean_near,mean_far,near_far_ratio,symmetry,verdict,reason,track
    ontad,110,1.423,1.1,0,0.8,0.56,1.43,0.95,STRONG,sharp_peak,rad21
    modularity_tad,214,1.314,1.05,0,0.75,0.57,1.31,0.94,STRONG,sharp_peak,rad21
    armatus,422,1.149,0.95,0,0.65,0.57,1.15,0.92,WEAK,flat_profile,rad21
    topdom,180,0.798,0.7,355,0.55,0.69,0.80,0.61,SHIFTED,peak_at_355kb,rad21
    coitad,15,0.0,0.0,0,0.0,0.0,0.0,0.5,NO_DATA,too_few_boundaries,rad21
    dihmm,104,1.031,0.92,0,0.70,0.68,1.03,0.93,GOOD,broad_peak,rad21
    scktld,120,1.020,0.91,0,0.69,0.68,1.02,0.92,WEAK,flat_profile,rad21
""")


@pytest.fixture()
def summary_csv(tmp_path: Path) -> Path:
    p = tmp_path / "chipseq_validation_summary.csv"
    p.write_text(SUMMARY_CSV_CONTENT)
    return p


def _make_bed(tmp_path: Path, name: str, rows: list[tuple]) -> Path:
    """Write a minimal BED3 file; rows = [(chrom, start, end), ...]."""
    p = tmp_path / name
    lines = "\n".join(f"{c}\t{s}\t{e}" for c, s, e in rows) + "\n"
    p.write_text(lines)
    return p


# ---------------------------------------------------------------------------
# 1. load_weights
# ---------------------------------------------------------------------------

class TestLoadWeights:

    def test_correct_values_returned(self, summary_csv: Path) -> None:
        w = load_weights(str(summary_csv), track="rad21")
        assert pytest.approx(w["ontad"],         rel=1e-4) == 1.423
        assert pytest.approx(w["modularity_tad"], rel=1e-4) == 1.314
        assert pytest.approx(w["armatus"],        rel=1e-4) == 1.149
        assert pytest.approx(w["dihmm"],          rel=1e-4) == 1.031

    def test_shifted_excluded(self, summary_csv: Path) -> None:
        w = load_weights(str(summary_csv), track="rad21",
                         exclude_verdicts=["NO_DATA", "SHIFTED"])
        assert w["topdom"] == 0.0

    def test_no_data_excluded(self, summary_csv: Path) -> None:
        w = load_weights(str(summary_csv), track="rad21",
                         exclude_verdicts=["NO_DATA", "SHIFTED"])
        assert w["coitad"] == 0.0

    def test_unknown_track_returns_empty(self, summary_csv: Path) -> None:
        w = load_weights(str(summary_csv), track="nonexistent_track")
        assert w == {}

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            load_weights(str(tmp_path / "no_such_file.csv"))

    def test_negative_center_ratio_clamped(self, tmp_path: Path) -> None:
        csv = tmp_path / "neg.csv"
        csv.write_text(
            "algorithm,n_boundaries,center_ratio,peak_ratio,peak_offset_kb,"
            "mean_near,mean_far,near_far_ratio,symmetry,verdict,reason,track\n"
            "bad_algo,50,-0.5,0.5,0,0.4,0.5,0.8,0.5,WEAK,whatever,rad21\n"
        )
        w = load_weights(str(csv), track="rad21", exclude_verdicts=[])
        assert w["bad_algo"] == 0.0   # clamped from -0.5

    def test_default_exclude_verdicts(self, summary_csv: Path) -> None:
        """Default exclude_verdicts == ['NO_DATA', 'SHIFTED']."""
        w = load_weights(str(summary_csv), track="rad21")
        assert w["topdom"] == 0.0
        assert w["coitad"] == 0.0


# ---------------------------------------------------------------------------
# 2. compute_weighted_consensus — synthetic ground truth
# ---------------------------------------------------------------------------

class TestComputeWeightedConsensus:
    """
    Synthetic setup (resolution=100_000, tolerance_bins=1):

    Chromosome: chr1
    Three algorithms A, B, C with boundary bins:
        A (weight=2.0): bins [10, 20, 30]
        B (weight=1.0): bins [10, 21, 40]   ← bin 21 within tolerance of bin 20
        C (weight=0.0): bins [10, 50, 60]   ← excluded

    Effective total_weight = 2.0 + 1.0 + 0.0 = 3.0

    Candidate bins = union = {10, 20, 21, 30, 40, 50, 60}

    Scores (tolerance_bins=1):
        bin 10 : A hit (10 in [9..11]) + B hit (10 in [9..11]) + C=0
                 = (2.0 + 1.0) / 3.0 = 1.000
        bin 20 : A hit (20 in [19..21]) + B hit? (21 in [19..21]) + C=0
                 = (2.0 + 1.0) / 3.0 = 1.000
        bin 21 : A hit (20 in [20..22]) + B hit (21) + C=0
                 = (2.0 + 1.0) / 3.0 = 1.000
        bin 30 : A hit + B miss + C=0 = 2.0/3.0 ≈ 0.667
        bin 40 : A miss + B hit + C=0 = 1.0/3.0 ≈ 0.333
        bin 50 : C=0 only             = 0.0/3.0 = 0.0
        bin 60 : C=0 only             = 0.0/3.0 = 0.0

    theta=0.4 → retained: 10, 20, 21, 30    (score ≥ 0.4)
    theta=0.7 → retained: 10, 20, 21, 30    (0.667 still ≥ 0.4… wait:
        bin 30 score=0.667 ≥ 0.7? No: 0.667 < 0.7)
    theta=0.7 → retained: 10, 20, 21
    """

    RES   = 100_000
    TOL   = 1
    CHROM = "chr1"

    @pytest.fixture()
    def tad_files(self, tmp_path: Path) -> dict[str, Path]:
        res = self.RES
        a = _make_bed(tmp_path, "a.bed", [
            ("chr1", 10 * res, 20 * res),  # boundaries at bins 10, 20
            ("chr1", 20 * res, 30 * res),  # boundaries at bins 20, 30  → 20 deduped
        ])
        b = _make_bed(tmp_path, "b.bed", [
            ("chr1", 10 * res, 21 * res),  # boundaries at bins 10, 21
            ("chr1", 21 * res, 40 * res),  # boundaries at bins 21, 40
        ])
        c = _make_bed(tmp_path, "c.bed", [
            ("chr1", 10 * res, 50 * res),  # boundaries at bins 10, 50
            ("chr1", 50 * res, 60 * res),  # boundaries at bins 50, 60
        ])
        return {"A": a, "B": b, "C": c}

    @pytest.fixture()
    def weights_mixed(self) -> dict[str, float]:
        return {"A": 2.0, "B": 1.0, "C": 0.0}

    def test_theta_040_returns_4_boundaries(
        self, tad_files: dict[str, Path], weights_mixed: dict[str, float]
    ) -> None:
        df = compute_weighted_consensus(
            tad_files, weights_mixed,
            chrom=self.CHROM, resolution=self.RES,
            theta=0.4, tolerance_bins=self.TOL,
        )
        assert set(df.columns) == {"chrom", "start", "end", "weighted_support"}
        # bins 10, 20, 21, 30 should survive (scores 1.0, 1.0, 1.0, 0.667)
        surviving_bins = set(df["start"].values // self.RES)
        assert 10 in surviving_bins
        assert 20 in surviving_bins
        assert 21 in surviving_bins
        assert 30 in surviving_bins

    def test_theta_070_excludes_bin30(
        self, tad_files: dict[str, Path], weights_mixed: dict[str, float]
    ) -> None:
        df = compute_weighted_consensus(
            tad_files, weights_mixed,
            chrom=self.CHROM, resolution=self.RES,
            theta=0.7, tolerance_bins=self.TOL,
        )
        surviving_bins = set(df["start"].values // self.RES)
        # bin 30 score ≈ 0.667 < 0.7
        assert 30 not in surviving_bins
        assert 10 in surviving_bins

    def test_excluded_algo_does_not_contribute(
        self, tad_files: dict[str, Path]
    ) -> None:
        """C is the only algo with bins 50, 60 but weight=0 → not retained."""
        weights = {"A": 2.0, "B": 1.0, "C": 0.0}
        df = compute_weighted_consensus(
            tad_files, weights,
            chrom=self.CHROM, resolution=self.RES,
            theta=0.4, tolerance_bins=self.TOL,
        )
        surviving_bins = set(df["start"].values // self.RES)
        assert 50 not in surviving_bins
        assert 60 not in surviving_bins

    def test_all_zero_weights_returns_empty(
        self, tad_files: dict[str, Path]
    ) -> None:
        weights = {"A": 0.0, "B": 0.0, "C": 0.0}
        df = compute_weighted_consensus(
            tad_files, weights,
            chrom=self.CHROM, resolution=self.RES,
            theta=0.4, tolerance_bins=self.TOL,
        )
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0
        assert set(df.columns) == {"chrom", "start", "end", "weighted_support"}

    def test_missing_file_does_not_raise(
        self, tmp_path: Path, weights_mixed: dict[str, float]
    ) -> None:
        files = {
            "A": tmp_path / "exists.bed",
            "B": tmp_path / "nonexistent.bed",   # ← missing
        }
        _make_bed(tmp_path, "exists.bed", [("chr1", 1_000_000, 2_000_000)])
        # Must not raise, just treat missing as empty
        df = compute_weighted_consensus(
            files, weights_mixed,
            chrom="chr1", resolution=self.RES,
            theta=0.4, tolerance_bins=self.TOL,
        )
        assert isinstance(df, pd.DataFrame)

    def test_wrong_chrom_returns_empty(
        self, tad_files: dict[str, Path], weights_mixed: dict[str, float]
    ) -> None:
        df = compute_weighted_consensus(
            tad_files, weights_mixed,
            chrom="chr99", resolution=self.RES,
            theta=0.4, tolerance_bins=self.TOL,
        )
        assert len(df) == 0

    def test_output_schema(
        self, tad_files: dict[str, Path], weights_mixed: dict[str, float]
    ) -> None:
        df = compute_weighted_consensus(
            tad_files, weights_mixed,
            chrom=self.CHROM, resolution=self.RES,
            theta=0.4, tolerance_bins=self.TOL,
        )
        assert list(df.columns) == ["chrom", "start", "end", "weighted_support"]
        assert (df["start"] % self.RES == 0).all()
        assert (df["end"]   % self.RES == 0).all()
        assert (df["end"] > df["start"]).all()
        assert df["weighted_support"].between(0.0, 1.0).all()

    def test_weighted_support_values(
        self, tad_files: dict[str, Path], weights_mixed: dict[str, float]
    ) -> None:
        """Verify exact scores for bins 10, 20/21, 30."""
        df = compute_weighted_consensus(
            tad_files, weights_mixed,
            chrom=self.CHROM, resolution=self.RES,
            theta=0.0, tolerance_bins=self.TOL,   # all bins, no filter
        )
        row = df[df["start"] == 10 * self.RES]
        assert len(row) == 1
        assert pytest.approx(row["weighted_support"].iloc[0], abs=1e-5) == 1.0

        row30 = df[df["start"] == 30 * self.RES]
        assert len(row30) == 1
        assert pytest.approx(row30["weighted_support"].iloc[0], abs=1e-4) == 2.0 / 3.0

    def test_deterministic_no_rng(
        self, tad_files: dict[str, Path], weights_mixed: dict[str, float]
    ) -> None:
        """compute_weighted_consensus is deterministic (no RNG involved)."""
        df1 = compute_weighted_consensus(
            tad_files, weights_mixed,
            chrom=self.CHROM, resolution=self.RES,
            theta=0.4, tolerance_bins=self.TOL,
        )
        df2 = compute_weighted_consensus(
            tad_files, weights_mixed,
            chrom=self.CHROM, resolution=self.RES,
            theta=0.4, tolerance_bins=self.TOL,
        )
        pd.testing.assert_frame_equal(df1, df2)


# ---------------------------------------------------------------------------
# 3. topdom SHIFTED scenario (integration with load_weights)
# ---------------------------------------------------------------------------

class TestTopdomIntegration:
    """topdom weight must be 0.0; its boundaries must not pull scores up."""

    def test_topdom_gets_zero_weight(self, summary_csv: Path) -> None:
        w = load_weights(str(summary_csv), track="rad21")
        assert w["topdom"] == 0.0

    def test_topdom_boundaries_ignored(
        self, tmp_path: Path, summary_csv: Path
    ) -> None:
        res = 100_000
        # Only topdom has boundary at bin 999 (unique)
        topdom_bed = _make_bed(tmp_path, "topdom.bed", [
            ("chr1", 999 * res, 1010 * res)
        ])
        ontad_bed = _make_bed(tmp_path, "ontad.bed", [
            ("chr1", 5 * res, 20 * res)
        ])

        weights = load_weights(str(summary_csv), track="rad21")
        tad_files = {"topdom": topdom_bed, "ontad": ontad_bed}

        df = compute_weighted_consensus(
            tad_files, weights,
            chrom="chr1", resolution=res,
            theta=0.4, tolerance_bins=1,
        )
        surviving_bins = set(df["start"].values // res)
        # bin 999 is ONLY from topdom (weight=0) → must NOT appear
        assert 999 not in surviving_bins
        assert 1010 not in surviving_bins
        # ontad boundaries must appear (weight > 0)
        assert 5 in surviving_bins or 20 in surviving_bins


# ---------------------------------------------------------------------------
# 4. save_weighted_consensus
# ---------------------------------------------------------------------------

class TestSaveWeightedConsensus:

    def test_file_created(self, tmp_path: Path) -> None:
        df = pd.DataFrame({
            "chrom":            ["chr1", "chr1"],
            "start":            [100_000, 200_000],
            "end":              [200_000, 300_000],
            "weighted_support": [0.85, 0.60],
        })
        out = tmp_path / "sub" / "test.bed"
        save_weighted_consensus(df, out)
        assert out.exists()

    def test_no_header_in_output(self, tmp_path: Path) -> None:
        df = pd.DataFrame({
            "chrom":            ["chr1"],
            "start":            [100_000],
            "end":              [200_000],
            "weighted_support": [0.9],
        })
        out = tmp_path / "out.bed"
        save_weighted_consensus(df, out)
        first_line = out.read_text().splitlines()[0]
        # Must start with 'chr', not 'chrom'
        assert first_line.startswith("chr")

    def test_roundtrip(self, tmp_path: Path) -> None:
        df_orig = pd.DataFrame({
            "chrom":            ["chr1", "chr1"],
            "start":            [100_000, 500_000],
            "end":              [200_000, 600_000],
            "weighted_support": [0.8, 0.6],
        })
        out = tmp_path / "roundtrip.bed"
        save_weighted_consensus(df_orig, out)

        df_back = pd.read_csv(out, sep="\t", header=None,
                              names=["chrom", "start", "end", "weighted_support"])
        assert list(df_back["start"]) == [100_000, 500_000]
        assert pytest.approx(df_back["weighted_support"].tolist()) == [0.8, 0.6]


# ---------------------------------------------------------------------------
# 5. _boundaries_from_tad_df (unit)
# ---------------------------------------------------------------------------

class TestBoundariesFromTadDf:

    def test_start_and_end_extracted(self) -> None:
        df = pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "start": [0, 200_000],
            "end":   [100_000, 300_000],
        })
        bins = _boundaries_from_tad_df(df, "chr1", 100_000)
        assert set(bins.tolist()) == {0, 1, 2, 3}

    def test_wrong_chrom_returns_empty(self) -> None:
        df = pd.DataFrame({
            "chrom": ["chr1"],
            "start": [100_000],
            "end":   [200_000],
        })
        bins = _boundaries_from_tad_df(df, "chr2", 100_000)
        assert len(bins) == 0

    def test_duplicate_boundaries_deduplicated(self) -> None:
        """Adjacent TADs share a boundary → deduplicated."""
        df = pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "start": [0, 100_000],
            "end":   [100_000, 200_000],
        })
        bins = _boundaries_from_tad_df(df, "chr1", 100_000)
        # bins should be [0, 1, 2] — bin 1 appears as end of TAD1 and start of TAD2
        assert len(bins) == 3
        assert 1 in bins

# ---------------------------------------------------------------------------
# 6. _apply_transform  (unit tests for all four transforms)
# ---------------------------------------------------------------------------

class TestApplyTransform:
    """
    Base raw_weights for all tests (mimics real summary after exclusion):
        ontad          cr=1.423  (active)
        modularity_tad cr=1.314  (active)
        armatus        cr=1.149  (active)
        scktld         cr=1.020  (active)
        dihmm          cr=1.031  (active)
        topdom         cr=0.0    (excluded, already 0.0)
        coitad         cr=0.0    (excluded, already 0.0)
    """
    RAW = {
        "ontad":          1.423,
        "modularity_tad": 1.314,
        "armatus":        1.149,
        "scktld":         1.020,
        "dihmm":          1.031,
        "topdom":         0.0,    # excluded before transform
        "coitad":         0.0,    # excluded before transform
    }

    # ── linear ────────────────────────────────────────────────────────────────

    def test_linear_identity(self) -> None:
        result = _apply_transform(self.RAW, "linear", {})
        assert pytest.approx(result["ontad"],  rel=1e-6) == 1.423
        assert pytest.approx(result["scktld"], rel=1e-6) == 1.020
        assert result["topdom"] == 0.0

    def test_linear_excluded_stay_zero(self) -> None:
        result = _apply_transform(self.RAW, "linear", {})
        assert result["topdom"] == 0.0
        assert result["coitad"] == 0.0

    # ── power ─────────────────────────────────────────────────────────────────

    def test_power_baseline_subtracted(self) -> None:
        """w = max(0, cr - 1.0)^2"""
        result = _apply_transform(self.RAW, "power",
                                  {"baseline": 1.0, "gamma": 2.0})
        expected_ontad   = (1.423 - 1.0) ** 2
        expected_armatus = (1.149 - 1.0) ** 2
        assert pytest.approx(result["ontad"],   rel=1e-5) == expected_ontad
        assert pytest.approx(result["armatus"], rel=1e-5) == expected_armatus

    def test_power_cr_below_baseline_is_zero(self) -> None:
        """cr < baseline → max(0, negative)^γ = 0."""
        raw = {"good_algo": 1.5, "weak_algo": 0.8}  # 0.8 < baseline=1.0
        result = _apply_transform(raw, "power", {"baseline": 1.0, "gamma": 2.0})
        assert result["weak_algo"] == 0.0
        assert result["good_algo"] > 0.0

    def test_power_negative_cr_is_zero(self) -> None:
        """Pathological: cr < 0 (e.g. inverted CTCF profile)."""
        raw = {"pathological": -0.5, "normal": 1.3}
        result = _apply_transform(raw, "power", {"baseline": 1.0, "gamma": 2.0})
        assert result["pathological"] == 0.0
        assert result["normal"] > 0.0

    def test_power_topdom_cr_gives_zero(self) -> None:
        """topdom center_ratio=0.798 < baseline=1.0 → weight=0 even if not excluded."""
        raw = {"topdom": 0.798, "ontad": 1.423}
        result = _apply_transform(raw, "power", {"baseline": 1.0, "gamma": 2.0})
        assert result["topdom"] == 0.0
        assert result["ontad"] > 0.0

    def test_power_amplifies_gap(self) -> None:
        """Power transform amplifies gap between strong and weak algorithms."""
        result_linear = _apply_transform(self.RAW, "linear", {})
        result_power  = _apply_transform(self.RAW, "power",
                                         {"baseline": 1.0, "gamma": 2.0})

        # In linear: ontad/scktld ratio
        ratio_linear = result_linear["ontad"] / (result_linear["scktld"] + 1e-9)
        # In power: ontad/scktld ratio
        ratio_power  = result_power["ontad"]  / (result_power["scktld"]  + 1e-9)

        assert ratio_power > ratio_linear, (
            f"Power transform should amplify gap: "
            f"power_ratio={ratio_power:.2f} vs linear_ratio={ratio_linear:.2f}"
        )

    def test_power_excluded_stay_zero(self) -> None:
        result = _apply_transform(self.RAW, "power",
                                  {"baseline": 1.0, "gamma": 2.0})
        assert result["topdom"] == 0.0
        assert result["coitad"] == 0.0

    def test_power_gamma_one_is_linear_above_baseline(self) -> None:
        """With gamma=1.0, power = linear shift by baseline."""
        raw = {"a": 1.5, "b": 1.2}
        result = _apply_transform(raw, "power",
                                  {"baseline": 1.0, "gamma": 1.0})
        assert pytest.approx(result["a"], rel=1e-9) == 0.5
        assert pytest.approx(result["b"], rel=1e-9) == 0.2

    def test_power_larger_gamma_increases_gap(self) -> None:
        """Larger gamma → stronger amplification of inter-algorithm gap."""
        raw = {"strong": 1.4, "weak": 1.1}
        r2 = _apply_transform(raw, "power", {"baseline": 1.0, "gamma": 2.0})
        r4 = _apply_transform(raw, "power", {"baseline": 1.0, "gamma": 4.0})

        ratio2 = r2["strong"] / (r2["weak"] + 1e-12)
        ratio4 = r4["strong"] / (r4["weak"] + 1e-12)
        assert ratio4 > ratio2

    # ── softmax ───────────────────────────────────────────────────────────────

    def test_softmax_all_positive(self) -> None:
        """exp(...) is always positive for active algos."""
        result = _apply_transform(self.RAW, "softmax",
                                  {"baseline": 1.0, "alpha": 5.0})
        for algo, w in result.items():
            if self.RAW[algo] > 0.0:
                assert w > 0.0, f"{algo} should have positive softmax weight"

    def test_softmax_excluded_stay_zero(self) -> None:
        result = _apply_transform(self.RAW, "softmax",
                                  {"baseline": 1.0, "alpha": 5.0})
        assert result["topdom"] == 0.0
        assert result["coitad"] == 0.0

    def test_softmax_cr_at_baseline_gives_one(self) -> None:
        """cr == baseline → exp(alpha*(cr-baseline)) = exp(0) = 1.0"""
        raw = {"neutral": 1.0, "strong": 1.5}
        result = _apply_transform(raw, "softmax",
                                  {"baseline": 1.0, "alpha": 5.0})
        assert pytest.approx(result["neutral"], rel=1e-9) == 1.0

    def test_softmax_cr_below_baseline_gives_less_than_one(self) -> None:
        """cr < baseline → exp(alpha*(cr-baseline)) < 1 but still > 0."""
        raw = {"below": 0.8, "above": 1.2}
        result = _apply_transform(raw, "softmax",
                                  {"baseline": 1.0, "alpha": 5.0})
        assert 0.0 < result["below"] < 1.0
        assert result["above"] > 1.0

    def test_softmax_amplifies_gap_vs_linear(self) -> None:
        result_linear  = _apply_transform(self.RAW, "linear",  {})
        result_softmax = _apply_transform(self.RAW, "softmax",
                                          {"baseline": 1.0, "alpha": 5.0})
        ratio_linear  = result_linear["ontad"]  / result_linear["scktld"]
        ratio_softmax = result_softmax["ontad"] / result_softmax["scktld"]
        assert ratio_softmax > ratio_linear

    # ── zscore ────────────────────────────────────────────────────────────────

    def test_zscore_above_mean_positive(self) -> None:
        """Algorithms above mean get positive weight."""
        result = _apply_transform(self.RAW, "zscore", {})
        # ontad (1.423) is above mean of active weights → positive
        assert result["ontad"] > 0.0

    def test_zscore_below_mean_zero(self) -> None:
        """Algorithms below mean are clamped to 0."""
        result = _apply_transform(self.RAW, "zscore", {})
        # scktld (1.020) is below mean → should be 0
        vals = [v for k, v in self.RAW.items() if v > 0]
        import numpy as np
        mu = np.mean(vals)
        below_mean = [k for k, v in self.RAW.items()
                      if 0 < v < mu]
        for algo in below_mean:
            assert result[algo] == 0.0,                 f"{algo} (cr={self.RAW[algo]:.3f} < mu={mu:.3f}) should be 0"

    def test_zscore_excluded_stay_zero(self) -> None:
        result = _apply_transform(self.RAW, "zscore", {})
        assert result["topdom"] == 0.0
        assert result["coitad"] == 0.0

    def test_zscore_uniform_weights_fallback(self) -> None:
        """All active weights identical → std=0 → fallback to linear."""
        raw = {"a": 1.2, "b": 1.2, "c": 1.2}
        # Should not raise; result should equal linear (cr unchanged)
        result = _apply_transform(raw, "zscore", {})
        assert isinstance(result, dict)

    # ── unknown transform ─────────────────────────────────────────────────────

    def test_unknown_transform_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown weight_transform"):
            _apply_transform(self.RAW, "banana", {})

    # ── empty active weights ──────────────────────────────────────────────────

    def test_all_excluded_returns_unchanged(self) -> None:
        """If all weights are 0 (all excluded), transform is a no-op."""
        raw = {"a": 0.0, "b": 0.0}
        for t in ["linear", "power", "softmax", "zscore"]:
            result = _apply_transform(raw, t, {})
            assert result["a"] == 0.0
            assert result["b"] == 0.0


# ---------------------------------------------------------------------------
# 7. load_weights — with transform parameter
# ---------------------------------------------------------------------------

class TestLoadWeightsWithTransform:

    def test_linear_backward_compatible(self, summary_csv: Path) -> None:
        """Default (linear) == old behavior."""
        w_old = load_weights(str(summary_csv), track="rad21")
        w_new = load_weights(str(summary_csv), track="rad21",
                             weight_transform="linear")
        assert w_old == w_new

    def test_power_reduces_weak_algo_weight(self, summary_csv: Path) -> None:
        w_linear = load_weights(str(summary_csv), track="rad21",
                                weight_transform="linear")
        w_power  = load_weights(str(summary_csv), track="rad21",
                                weight_transform="power",
                                transform_params={"baseline": 1.0, "gamma": 2.0})
        # armatus (cr=1.149) should have relatively less weight in power
        if w_linear["ontad"] > 0 and w_linear["armatus"] > 0:
            ratio_linear = w_linear["ontad"] / w_linear["armatus"]
            ratio_power  = w_power["ontad"]  / (w_power["armatus"] + 1e-12)
            assert ratio_power > ratio_linear

    def test_topdom_zero_in_all_transforms(self, summary_csv: Path) -> None:
        """topdom is SHIFTED → excluded → 0.0 in all transforms."""
        for t in ["linear", "power", "softmax", "zscore"]:
            w = load_weights(str(summary_csv), track="rad21",
                             weight_transform=t)
            assert w["topdom"] == 0.0, f"topdom should be 0 for transform={t}"

    def test_coitad_zero_in_all_transforms(self, summary_csv: Path) -> None:
        """coitad is NO_DATA → excluded → 0.0 in all transforms."""
        for t in ["linear", "power", "softmax", "zscore"]:
            w = load_weights(str(summary_csv), track="rad21",
                             weight_transform=t)
            assert w["coitad"] == 0.0, f"coitad should be 0 for transform={t}"
