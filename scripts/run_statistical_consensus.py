from __future__ import annotations

"""Orchestrator for statistical TAD boundary consensus methods.

Usage
-----
python scripts/run_statistical_consensus.py \
    --chroms chr1 \
    --resolution 100000 \
    --tracks rad21 smc3 \
    --methods weighted permutation pelt

Steps executed for each (chrom, resolution) pair
-------------------------------------------------
1. weighted   — compute_weighted_consensus()  → BED4
2. permutation— compute_boundary_enrichment() → CSV + BED3  (per track)
3. pelt       — run_pelt_tad()               → BED3  (proof-of-concept)
4. For every produced BED: compute center_ratio via
   run_ctcf_profile_analysis() from src/ctcf_analysis.py
5. Append to results/stats/statistical_consensus_comparison.csv
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import yaml

# ── project root on sys.path (same mechanism as conftest.py) ────────────────
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.statistical_methods.weighted_consensus import (
    compute_weighted_consensus,
    load_weights,
    save_weighted_consensus,
)
from src.statistical_methods.permutation_chipseq import (
    compute_boundary_enrichment,
    save_permutation_results,
)
from src.statistical_methods.changepoint_tad import run_pelt_tad

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s — %(message)s",
    datefmt="%H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger("run_statistical_consensus")


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------

def load_config(config_path: str = "config/config.yaml") -> dict:
    p = Path(config_path)
    if not p.exists():
        raise FileNotFoundError(f"Config not found: {p}")
    with p.open() as fh:
        return yaml.safe_load(fh)


def _path(cfg: dict, key: str, fallback: str) -> Path:
    """Resolve a path from cfg['paths'][key], falling back to fallback."""
    return Path(cfg.get("paths", {}).get(key, fallback))


def _chrom_size(cfg: dict, chrom: str) -> int:
    """Return chromosome size from cfg['chrom_sizes'] or hg19 defaults."""
    sizes = cfg.get("chrom_sizes", {})
    if chrom in sizes:
        return int(sizes[chrom])
    # hg19 defaults for the chromosomes used in the project
    HG19 = {
        "chr1": 249_250_621, "chr2": 243_199_373, "chr3": 198_022_430,
        "chr4": 191_154_276, "chr5": 180_915_260, "chr6": 171_115_067,
        "chr7": 159_138_663, "chr8": 146_364_022, "chr9": 141_213_431,
        "chr10": 135_534_747, "chr11": 135_006_516, "chr12": 133_851_895,
        "chr13": 115_169_878, "chr14": 107_349_540, "chr15": 102_531_392,
        "chr16": 90_354_753,  "chr17": 81_195_210,  "chr18": 78_077_248,
        "chr19": 59_128_983,  "chr20": 63_025_520,  "chr21": 48_129_895,
        "chr22": 51_304_566,  "chrX": 155_270_560,
    }
    if chrom in HG19:
        return HG19[chrom]
    raise ValueError(
        f"Chromosome size for {chrom} not found in config or hg19 defaults"
    )


# ---------------------------------------------------------------------------
# ChIP-seq profile helper
# ---------------------------------------------------------------------------

def _center_ratio_for_bed(
    bed_path: Path,
    track_name: str,
    chrom: str,
    resolution: int,
    cfg: dict,
) -> Optional[float]:
    """Compute center_ratio for a BED file using run_ctcf_profile_analysis.

    Matches the real signature of run_ctcf_profile_analysis():
        run_ctcf_profile_analysis(
            algo_results: Dict[str, pd.DataFrame],
            ctcf_df: pd.DataFrame,
            chrom, resolution, ...
        )

    Returns None if analysis fails or BED is empty / < 30 boundaries.
    """
    try:
        from src.ctcf_analysis import run_ctcf_profile_analysis  # noqa: PLC0415
    except ImportError:
        logger.warning("src.ctcf_analysis not available — skipping enrichment")
        return None

    # ── Read TAD/boundary BED (first 3 columns only) ─────────────────────────
    if not bed_path.exists():
        logger.warning("BED not found: %s — skipping center_ratio", bed_path)
        return None

    try:
        tad_df = pd.read_csv(
            str(bed_path), sep="\t", header=None,
            usecols=[0, 1, 2], names=["chrom", "start", "end"],
            dtype={"chrom": str, "start": int, "end": int},
        )
        tad_df = tad_df[tad_df["chrom"] == chrom].reset_index(drop=True)
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to read BED %s: %s", bed_path, exc)
        return None

    if tad_df.empty:
        logger.warning("No rows for %s in %s — skipping", chrom, bed_path.name)
        return None

    # ── Read ChIP-seq BED ────────────────────────────────────────────────────
    track_key    = f"{track_name}_bed"
    chipseq_path = _path(
        cfg, track_key,
        f"data/reference/GM12878_{track_name.upper()}_peaks_hg19.bed",
    )
    if not chipseq_path.exists():
        logger.warning("ChIP-seq file not found: %s — skipping", chipseq_path)
        return None

    try:
        ctcf_df = pd.read_csv(
            str(chipseq_path), sep="\t", header=None,
            usecols=[0, 1, 2], names=["chrom", "start", "end"],
            comment="#",
        )
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to read ChIP-seq %s: %s", chipseq_path, exc)
        return None

    # ── Call with real signature ─────────────────────────────────────────────
    method_label = bed_path.stem   # e.g. "weighted_consensus_chr1_100000bp"
    try:
        result_df = run_ctcf_profile_analysis(
            algo_results={method_label: tad_df},
            ctcf_df=ctcf_df,
            chrom=chrom,
            resolution=resolution,
        )
    except Exception as exc:  # noqa: BLE001
        logger.warning(
            "center_ratio computation failed for %s / %s: %s",
            bed_path.name, track_name, exc,
        )
        return None

    if result_df is None or result_df.empty:
        return None

    # Extract the row for our method (first row if label not found)
    row = result_df[result_df["algorithm"] == method_label]
    if row.empty:
        row = result_df.iloc[[0]]

    val = row["center_ratio"].iloc[0]
    return float(val) if pd.notna(val) else None


# ---------------------------------------------------------------------------
# Method runners
# ---------------------------------------------------------------------------

def run_weighted(
    chrom: str,
    resolution: int,
    cfg: dict,
    tracks: list[str],
) -> list[dict]:
    """Run Method 3: weighted consensus. Returns list of result-row dicts."""
    rows: list[dict] = []

    wc_cfg   = cfg.get("weighted_consensus", {})
    theta    = float(wc_cfg.get("theta", 0.4))
    w_track  = str(wc_cfg.get("weight_track", "rad21"))
    excl     = list(wc_cfg.get("exclude_verdicts", ["NO_DATA", "SHIFTED"]))
    tol_bins = int(cfg.get("consensus", {}).get("tolerance_bins", 1))

    summary_csv = str(_path(cfg, "chipseq_validation_summary",
                            "results/stats/chipseq_validation_summary.csv"))

    # Load weights (with nonlinear transform from config)
    w_transform = str(wc_cfg.get("weight_transform", "linear"))
    t_params    = dict(wc_cfg.get("transform_params", {}))
    try:
        weights = load_weights(
            summary_csv,
            track=w_track,
            exclude_verdicts=excl,
            weight_transform=w_transform,
            transform_params=t_params,
        )
    except FileNotFoundError:
        logger.warning(
            "[weighted] summary CSV not found: %s — skipping weighted method",
            summary_csv,
        )
        return rows

    # Collect TAD BED files
    tad_dir   = Path(cfg.get("paths", {}).get("tads_out", "results/tads"))
    algo_list = list(cfg.get("algorithms", {}).keys())
    if not algo_list:
        # Fallback: glob existing BED files for this chrom/res
        pattern = f"*_{chrom}_{resolution}bp.bed"
        algo_list = [
            p.name.replace(f"_{chrom}_{resolution}bp.bed", "")
            for p in tad_dir.glob(pattern)
        ]

    tad_files: dict[str, Path] = {}
    for algo in algo_list:
        bed = tad_dir / f"{algo}_{chrom}_{resolution}bp.bed"
        if bed.exists():
            tad_files[algo] = bed

    if not tad_files:
        logger.warning(
            "[weighted] no TAD BED files found in %s for %s@%dbp",
            tad_dir, chrom, resolution,
        )
        return rows

    df = compute_weighted_consensus(
        tad_files=tad_files,
        weights=weights,
        chrom=chrom,
        resolution=resolution,
        theta=theta,
        tolerance_bins=tol_bins,
    )

    out_dir  = Path(cfg.get("paths", {}).get("consensus_out", "results/consensus"))
    out_path = out_dir / f"weighted_consensus_{chrom}_{resolution}bp.bed"
    save_weighted_consensus(df, out_path)

    logger.info(
        "[weighted_consensus] %s@%dbp: %d boundaries, theta=%.2f "
        "(transform=%s)",
        chrom, resolution, len(df), theta, w_transform,
    )

    # Compute enrichment for each requested track
    for track in tracks:
        cr = _center_ratio_for_bed(out_path, track, chrom, resolution, cfg)
        rows.append({
            "method":        "weighted_consensus",
            "chrom":         chrom,
            "resolution_bp": resolution,
            "track":         track,
            "n_boundaries":  len(df),
            "theta":         theta,
            "weight_track":  w_track,
            "center_ratio":  cr,
            "output_file":   str(out_path),
        })

    return rows


def run_permutation(
    chrom: str,
    resolution: int,
    cfg: dict,
    tracks: list[str],
    source_bed: Optional[Path] = None,
) -> list[dict]:
    """Run Method 4: permutation test.

    source_bed: boundaries to test (default: weak consensus BED).
    Returns list of result-row dicts.
    """
    rows: list[dict] = []
    perm_cfg = cfg.get("permutation", {})
    window_bp    = int(perm_cfg.get("window_bp",          50_000))
    n_perm       = int(perm_cfg.get("n_permutations",     1_000))
    fdr          = float(perm_cfg.get("fdr_threshold",    0.05))
    telo         = int(perm_cfg.get("telomere_margin_bp", 3_000_000))
    rng          = np.random.default_rng(42)
    chrom_size   = _chrom_size(cfg, chrom)

    # Default source: boundary consensus (BED of boundary positions, not domains).
    # consensus_{chrom}_{res}bp.bed  — жадная кластеризация границ (§6 rules.md)
    # НЕ tad_consensus_...bed (домены) — их центры находятся внутри TAD,
    # а не на границах, что даёт нулевое обогащение ChIP-seq.
    if source_bed is None:
        cons_dir   = Path(cfg.get("paths", {}).get("consensus_out", "results/consensus"))
        source_bed = cons_dir / f"consensus_{chrom}_{resolution}bp.bed"

    if not source_bed.exists():
        logger.warning(
            "[permutation] source BED not found: %s — skipping", source_bed
        )
        return rows

    # Read source boundaries
    try:
        bnd_df = pd.read_csv(source_bed, sep="\t", header=None,
                             usecols=[0, 1, 2],
                             names=["chrom", "start", "end"])
        bnd_df = bnd_df[bnd_df["chrom"] == chrom].copy()
    except Exception as exc:  # noqa: BLE001
        logger.warning("[permutation] failed to read %s: %s", source_bed, exc)
        return rows

    for track in tracks:
        track_key    = f"{track}_bed"
        chipseq_path = _path(cfg, track_key,
                             f"data/reference/GM12878_{track.upper()}_peaks_hg19.bed")

        if not chipseq_path.exists():
            logger.warning(
                "[permutation] ChIP-seq file not found: %s — skip track %s",
                chipseq_path, track,
            )
            continue

        try:
            cs_df = pd.read_csv(str(chipseq_path), sep="\t", header=None,
                                usecols=[0, 1, 2],
                                names=["chrom", "start", "end"],
                                comment="#")
        except Exception as exc:  # noqa: BLE001
            logger.warning("[permutation] failed to read ChIP-seq %s: %s",
                           chipseq_path, exc)
            continue

        result_df = compute_boundary_enrichment(
            boundaries_df=bnd_df,
            chipseq_df=cs_df,
            chrom=chrom,
            chrom_size=chrom_size,
            window_bp=window_bp,
            n_permutations=n_perm,
            rng=rng,
            telomere_margin_bp=telo,
            fdr_threshold=fdr,
        )

        stats_dir = Path(cfg.get("paths", {}).get("stats_out", "results/stats"))
        cons_dir  = Path(cfg.get("paths", {}).get("consensus_out", "results/consensus"))

        csv_path = stats_dir / f"permutation_{track}_{chrom}_{resolution}bp.csv"
        bed_path = cons_dir  / f"permutation_filtered_{track}_{chrom}_{resolution}bp.bed"
        save_permutation_results(result_df, csv_path, bed_path)

        n_sig = int(result_df["significant"].sum()) if len(result_df) > 0 else 0
        logger.info(
            "[permutation] %s@%dbp: %d significant boundaries "
            "(FDR<%.2f, %s)",
            chrom, resolution, n_sig, fdr, track,
        )

        cr = _center_ratio_for_bed(bed_path, track, chrom, resolution, cfg)
        rows.append({
            "method":        f"permutation_{track}",
            "chrom":         chrom,
            "resolution_bp": resolution,
            "track":         track,
            "n_boundaries":  n_sig,
            "theta":         fdr,
            "weight_track":  "n/a",
            "center_ratio":  cr,
            "output_file":   str(bed_path),
        })

    return rows


def run_pelt(
    chrom: str,
    resolution: int,
    cfg: dict,
    tracks: list[str],
) -> list[dict]:
    """Run Method 1: PELT change-point detection. Returns list of result-row dicts."""
    rows: list[dict] = []

    df = run_pelt_tad(
        chrom=chrom,
        resolution=resolution,
        data_path=str(_path(cfg, "hic_data", "data/raw")),
        cfg=cfg,
    )

    n_bnd = len(df)
    logger.info("[pelt] %s@%dbp: %d boundaries (proof-of-concept)", chrom, resolution, n_bnd)

    if n_bnd == 0:
        for track in tracks:
            rows.append({
                "method":        "pelt",
                "chrom":         chrom,
                "resolution_bp": resolution,
                "track":         track,
                "n_boundaries":  0,
                "theta":         float("nan"),
                "weight_track":  "n/a",
                "center_ratio":  float("nan"),
                "output_file":   "",
            })
        return rows

    cons_dir = Path(cfg.get("paths", {}).get("consensus_out", "results/consensus"))
    out_path = cons_dir / f"pelt_consensus_{chrom}_{resolution}bp.bed"
    cons_dir.mkdir(parents=True, exist_ok=True)
    df[["chrom", "start", "end"]].to_csv(
        out_path, sep="\t", header=False, index=False
    )
    logger.info("Saved PELT BED → %s", out_path)

    for track in tracks:
        cr = _center_ratio_for_bed(out_path, track, chrom, resolution, cfg)
        rows.append({
            "method":        "pelt",
            "chrom":         chrom,
            "resolution_bp": resolution,
            "track":         track,
            "n_boundaries":  n_bnd,
            "theta":         float("nan"),
            "weight_track":  "n/a",
            "center_ratio":  cr,
            "output_file":   str(out_path),
        })

    return rows


# ---------------------------------------------------------------------------
# Summary I/O
# ---------------------------------------------------------------------------

def _update_summary(new_rows: list[dict], summary_path: Path) -> pd.DataFrame:
    """Append new_rows to summary CSV, deduplicating by (method, chrom, resolution_bp, track)."""
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    new_df = pd.DataFrame(new_rows)

    if summary_path.exists():
        existing = pd.read_csv(summary_path)
        key_cols = ["method", "chrom", "resolution_bp", "track"]
        # Drop old rows that are being overwritten
        mask = pd.Series([True] * len(existing))
        for _, row in new_df.iterrows():
            match = pd.Series([True] * len(existing))
            for col in key_cols:
                if col in existing.columns:
                    match &= (existing[col] == row[col])
            mask &= ~match
        combined = pd.concat([existing[mask], new_df], ignore_index=True)
    else:
        combined = new_df

    combined.to_csv(summary_path, index=False)
    return combined


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run statistical TAD boundary consensus methods.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
# All three methods, chr1, 100kb, two tracks
python scripts/run_statistical_consensus.py \\
    --chroms chr1 --resolution 100000 \\
    --tracks rad21 smc3 --methods weighted permutation pelt

# Only weighted + permutation, multiple chromosomes
python scripts/run_statistical_consensus.py \\
    --chroms chr1 chr17 chr18 --resolution 50000 \\
    --tracks rad21 --methods weighted permutation
""",
    )
    p.add_argument(
        "--chroms", nargs="+", required=True,
        help="Chromosomes to process (e.g. chr1 chr17)",
    )
    p.add_argument(
        "--resolution", type=int, default=100_000,
        help="Resolution in bp (default: 100000)",
    )
    p.add_argument(
        "--tracks", nargs="+", default=["rad21"],
        choices=["rad21", "smc3", "h3k4me3", "h3k27ac", "ctcf"],
        help="ChIP-seq tracks for enrichment validation (default: rad21)",
    )
    p.add_argument(
        "--methods", nargs="+", default=["weighted", "permutation", "pelt"],
        choices=["weighted", "permutation", "pelt"],
        help="Methods to run (default: all three)",
    )
    p.add_argument(
        "--config", default="config/config.yaml",
        help="Path to config YAML (default: config/config.yaml)",
    )
    p.add_argument(
        "--summary-csv", default="results/stats/statistical_consensus_comparison.csv",
        help="Output summary CSV path",
    )
    return p.parse_args()


def main() -> None:
    args   = _parse_args()
    cfg    = load_config(args.config)
    all_rows: list[dict] = []

    logger.info(
        "Starting statistical consensus | chroms=%s | res=%dbp | "
        "methods=%s | tracks=%s",
        args.chroms, args.resolution, args.methods, args.tracks,
    )

    for chrom in args.chroms:
        logger.info("─── %s @ %dbp ───", chrom, args.resolution)

        if "weighted" in args.methods:
            try:
                rows = run_weighted(chrom, args.resolution, cfg, args.tracks)
                all_rows.extend(rows)
            except Exception as exc:  # noqa: BLE001
                logger.error("[weighted] %s: unexpected error — %s", chrom, exc)

        if "permutation" in args.methods:
            try:
                rows = run_permutation(chrom, args.resolution, cfg, args.tracks)
                all_rows.extend(rows)
            except Exception as exc:  # noqa: BLE001
                logger.error("[permutation] %s: unexpected error — %s", chrom, exc)

        if "pelt" in args.methods:
            try:
                rows = run_pelt(chrom, args.resolution, cfg, args.tracks)
                all_rows.extend(rows)
            except Exception as exc:  # noqa: BLE001
                logger.error("[pelt] %s: unexpected error — %s", chrom, exc)

    # ── Summary ──────────────────────────────────────────────────────────────
    if all_rows:
        summary_path = Path(args.summary_csv)
        combined = _update_summary(all_rows, summary_path)
        logger.info(
            "Summary saved → %s (%d rows total)",
            summary_path, len(combined),
        )
        # Pretty-print to stdout
        display_cols = [
            "method", "chrom", "resolution_bp", "track",
            "n_boundaries", "center_ratio",
        ]
        cols = [c for c in display_cols if c in combined.columns]
        print("\n" + combined[cols].to_string(index=False))
    else:
        logger.warning("No results produced — summary not written")

    logger.info("Done.")


if __name__ == "__main__":
    main()