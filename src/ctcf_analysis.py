"""
ctcf_analysis.py  v2
====================
Количественный анализ CTCF-профилей TAD-границ.

Метрики:
  center_ratio   — mean(density в ±50kb) / mean(всего профиля)
                   > 1.0 = обогащение у границы
  peak_ratio     — max / mean  (>1.5 = есть структура)
  peak_offset_kb — смещение пика от 0 (идеал = 0)
  near_far_ratio — mean_density(±100kb) / mean_density(200–500kb)
                   > 1.0 = концентрация у границы  [НОРМИРОВАНО по зоне]
  symmetry       — min(left_auc,right_auc)/max(...)  (1.0 = симметрично)
  n_boundaries   — число границ (< 30 → NO_DATA)

Вердикт:
  STRONG  — убедительный пик у 0 (near_far>1.2 И center>1.1 И |offset|<100kb)
  WEAK    — слабый пик у 0 (center>1.05 ИЛИ near_far>1.05)
  SHIFTED — пик есть, но далеко от 0
  FLAT    — нет структуры
  NO_DATA — мало данных
"""

from __future__ import annotations
import logging, os
from typing import Dict, Optional
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def analyse_ctcf_profile(
    bins: np.ndarray,
    density: np.ndarray,
    algo: str,
    n_boundaries: int,
    profile_range_bp: int = 500_000,
    near_window_bp:   int = 100_000,
    center_window_bp: int =  50_000,
    min_boundaries:   int = 30,
) -> dict:
    """Рассчитать метрики одного CTCF-профиля."""

    dens = density.astype(float)
    mean_d = dens.mean()
    if mean_d < 1e-12:
        return _empty(algo, n_boundaries, "zero_density")

    # ── Нормировка на среднее ─────────────────────────────────────────────────
    dens_norm = dens / mean_d   # 1.0 = случайный уровень

    # ── center_ratio: среднее в ±center_window вместо одной точки ───────────
    center_mask = np.abs(bins) <= center_window_bp
    if center_mask.sum() == 0:
        center_mask = np.array([np.argmin(np.abs(bins))])  # fallback
    center_ratio = float(dens_norm[center_mask].mean())

    # ── peak_ratio и peak_offset ─────────────────────────────────────────────
    peak_idx       = int(np.argmax(dens_norm))
    peak_ratio     = float(dens_norm[peak_idx])
    peak_offset_kb = float(bins[peak_idx] / 1000)

    # ── near_far_ratio: СРЕДНЯЯ плотность (не интеграл) ─────────────────────
    near_mask = np.abs(bins) <= near_window_bp
    far_mask  = (np.abs(bins) > 2 * near_window_bp) & \
                (np.abs(bins) <= profile_range_bp)

    mean_near = float(dens_norm[near_mask].mean()) if near_mask.sum() > 0 else 0.0
    mean_far  = float(dens_norm[far_mask].mean())  if far_mask.sum()  > 0 else 1.0
    near_far_ratio = mean_near / (mean_far + 1e-9)

    # ── Симметрия ─────────────────────────────────────────────────────────────
    left_mask  = (bins >= -near_window_bp) & (bins < 0)
    right_mask = (bins > 0) & (bins <= near_window_bp)
    ml = float(dens_norm[left_mask].mean())  if left_mask.sum()  > 0 else 0.0
    mr = float(dens_norm[right_mask].mean()) if right_mask.sum() > 0 else 0.0
    symmetry = min(ml, mr) / (max(ml, mr) + 1e-9) if max(ml, mr) > 1e-9 else 0.0

    # ── Вердикт (калиброван по Arrowhead chr1@100kb) ─────────────────────────
    # Arrowhead-эталон: center_ratio≈1.08, near_far≈1.04
    # STRONG  ≥ Arrowhead (center≥1.07 AND near_far≥1.03 AND |peak|≤100kb)
    # GOOD    ≥ 90% Arrowhead (center≥0.97 AND near_far≥0.93)
    # WEAK    слабее Arrowhead но есть сигнал (center≥1.03 OR near_far≥1.02)
    # SHIFTED структура есть, но пик далеко от границы
    # FLAT    нет сигнала
    if n_boundaries < min_boundaries:
        verdict, reason = "NO_DATA", f"n_bnd={n_boundaries}<{min_boundaries}"
    elif center_ratio >= 1.07 and near_far_ratio >= 1.03 and abs(peak_offset_kb) <= 100:
        verdict, reason = "STRONG", f"≥Arrowhead (ctr={center_ratio:.2f}, n/f={near_far_ratio:.2f})"
    elif center_ratio >= 0.97 and near_far_ratio >= 0.93 and abs(peak_offset_kb) <= 150:
        verdict, reason = "GOOD",   f"~Arrowhead уровень (ctr={center_ratio:.2f})"
    elif center_ratio >= 1.03 or near_far_ratio >= 1.02:
        verdict, reason = "WEAK",   f"слабый сигнал (ctr={center_ratio:.2f})"
    elif peak_ratio >= 1.3 and abs(peak_offset_kb) > 150:
        verdict, reason = "SHIFTED", f"пик смещён на {peak_offset_kb:+.0f} kb"
    else:
        verdict, reason = "FLAT", f"нет структуры (ctr={center_ratio:.2f})"

    return {
        "algorithm":      algo,
        "n_boundaries":   n_boundaries,
        "center_ratio":   round(center_ratio,    3),
        "peak_ratio":     round(peak_ratio,      3),
        "peak_offset_kb": round(peak_offset_kb,  1),
        "mean_near":      round(mean_near,        3),
        "mean_far":       round(mean_far,         3),
        "near_far_ratio": round(near_far_ratio,   3),
        "symmetry":       round(symmetry,         3),
        "verdict":        verdict,
        "reason":         reason,
    }


def _empty(algo: str, n: int, reason: str) -> dict:
    return dict(algorithm=algo, n_boundaries=n,
                center_ratio=np.nan, peak_ratio=np.nan,
                peak_offset_kb=np.nan, mean_near=np.nan,
                mean_far=np.nan, near_far_ratio=np.nan,
                symmetry=np.nan, verdict="NO_DATA", reason=reason)


def run_ctcf_profile_analysis(
    algo_results: Dict[str, pd.DataFrame],
    ctcf_df: pd.DataFrame,
    chrom: str,
    resolution: int,
    arrowhead_path: Optional[str] = None,
    profile_range_bp: int = 500_000,
    profile_bin_bp:   int = 10_000,
    near_window_bp:   int = 100_000,
    center_window_bp: int =  50_000,
    min_boundaries:   int = 30,
    out_csv: Optional[str] = None,
) -> pd.DataFrame:
    """Анализ для всех алгоритмов + Arrowhead."""
    from src.validation import compute_ctcf_profile
    from src.consensus import extract_boundaries

    rows = []

    # ── Arrowhead ─────────────────────────────────────────────────────────────
    if arrowhead_path and os.path.exists(arrowhead_path):
        try:
            ah = pd.read_csv(arrowhead_path, sep="\t", comment="#",
                             header=0, dtype=str).iloc[:, :3]
            ah.columns = ["chrom", "start", "end"]
            ah["start"] = pd.to_numeric(ah["start"], errors="coerce")
            ah["end"]   = pd.to_numeric(ah["end"],   errors="coerce")
            ah = ah.dropna(subset=["start", "end"])
            if not str(ah["chrom"].iloc[0]).startswith("chr"):
                ah["chrom"] = "chr" + ah["chrom"].astype(str)
            ah = ah[ah["chrom"] == chrom].reset_index(drop=True)
            if not ah.empty:
                bins, dens = compute_ctcf_profile(
                    ah, ctcf_df, chrom, resolution,
                    profile_range_bp=profile_range_bp,
                    profile_bin_bp=profile_bin_bp,
                )
                bnd = extract_boundaries(ah, resolution)
                rows.append(analyse_ctcf_profile(
                    bins, dens, "Arrowhead", len(bnd),
                    profile_range_bp, near_window_bp, center_window_bp, min_boundaries,
                ))
                logger.info("Arrowhead: %d доменов на %s", len(ah), chrom)
        except Exception as exc:
            logger.warning("Arrowhead пропущен: %s", exc)

    # ── Алгоритмы ─────────────────────────────────────────────────────────────
    for algo, df in algo_results.items():
        if df is None or df.empty:
            rows.append(_empty(algo, 0, "empty_df"))
            continue
        try:
            bins, dens = compute_ctcf_profile(
                df, ctcf_df, chrom, resolution,
                profile_range_bp=profile_range_bp,
                profile_bin_bp=profile_bin_bp,
            )
            bnd = extract_boundaries(df, resolution)
            rows.append(analyse_ctcf_profile(
                bins, dens, algo, len(bnd),
                profile_range_bp, near_window_bp, center_window_bp, min_boundaries,
            ))
        except Exception as exc:
            logger.warning("[%s] пропущен: %s", algo, exc)
            rows.append(_empty(algo, 0, str(exc)))

    df_out = pd.DataFrame(rows)
    _print_report(df_out, chrom, resolution)

    if out_csv:
        os.makedirs(os.path.dirname(out_csv), exist_ok=True)
        df_out.to_csv(out_csv, index=False)
        logger.info("Сохранено: %s", out_csv)

    return df_out


def _print_report(df: pd.DataFrame, chrom: str, res: int) -> None:
    ICON = {"STRONG": "🟢", "GOOD": "🟢", "WEAK": "🟡", "SHIFTED": "🟠",
            "FLAT": "🔴", "NO_DATA": "⚪"}
    print(f"\n{'═'*82}")
    print(f"  CTCF-профиль | {chrom} @ {res//1000} kb")
    print(f"{'═'*82}")
    print(f"{'Algorithm':18s} {'n_bnd':>5} {'ctr_r':>6} {'pk_r':>6} "
          f"{'peak@kb':>8} {'near':>6} {'far':>6} {'n/f':>6} {'sym':>5}  Verdict")
    print(f"{'-'*82}")

    order = {"STRONG":0,"GOOD":1,"WEAK":2,"SHIFTED":3,"FLAT":4,"NO_DATA":5}
    df_s = df.copy()
    df_s["_o"] = df_s["verdict"].map(order).fillna(5)
    df_s = df_s.sort_values(["_o","center_ratio"], ascending=[True,False])

    for _, r in df_s.iterrows():
        def f(v, fmt=".3f"):
            return format(v, fmt) if pd.notna(v) else " — "
        print(f"{r['algorithm']:18s} {int(r['n_boundaries']):5d} "
              f"{f(r['center_ratio']):>6} {f(r['peak_ratio']):>6} "
              f"{f(r['peak_offset_kb'],'+.0f'):>8} "
              f"{f(r['mean_near']):>6} {f(r['mean_far']):>6} "
              f"{f(r['near_far_ratio']):>6} {f(r['symmetry']):>5}  "
              f"{ICON.get(r['verdict'],'?')} {r['verdict']} ({r['reason']})")

    print(f"{'═'*82}")
    print("near/far = средняя_плотность(±100kb) / средняя_плотность(±200-500kb)")
    print("ctr_r    = средняя_плотность(±50kb) / средняя_по_профилю\n")

    for v, label in [("STRONG","Сильный сигнал (≥Arrowhead)"),("GOOD","Хороший (~Arrowhead)"),
                     ("WEAK","Слабый сигнал"),("SHIFTED","Пик смещён"),("FLAT","Случайные границы")]:
        algos = df[df["verdict"]==v]["algorithm"].tolist()
        if algos:
            icons = {"STRONG":"✅","GOOD":"✅","WEAK":"⚠️","SHIFTED":"🔶","FLAT":"❌"}
            print(f"  {icons[v]} {label}: {', '.join(algos)}")
    print()
