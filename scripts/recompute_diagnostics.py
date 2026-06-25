#!/usr/bin/env python
"""scripts/recompute_diagnostics.py — recompute the manuscript-facing diagnostics
from per-seed severity-detail h5ads.

This is Phase 1 of the training-aware reproducibility pipeline. It takes the
already-evaluated severity-detail artefacts (produced by Phase 2,
scripts/evaluate_severity_panel.py) and produces exactly the eval, table, and
figure-input CSVs that the figure notebooks (01-05) consume.

Input contract
--------------
For each (model, cell_type, split_type, seed), one h5ad at:

    {predictions-dir}/{model}_{cell_type}_{split_type}_seed{seed}.severity.h5ad

with .obs columns:

    predicted_mean_abs_delta : float, per-perturbation predicted magnitude
    leverage_score           : float, per-perturbation severity reference value
    perturbation_target      : str,   gene symbol (used for LOO mode-driver)

The full cell-type-specific severity reference is read separately for the
winsorisation thresholds:

    {severity-refs-dir}/replogle_{K562|RPE1}_severity.csv

Output contract
---------------
data/eval/diag_loo_sensitivity_n100.csv
data/eval/diag_loo_sensitivity_gears.csv
data/eval/diag_winsorise_n100.csv
data/eval/diag_winsorise_n100_summary.csv
data/eval/diag_spearman_n100.csv             (§3.3 metric robustness)
data/eval/diag_spearman_n100_summary.csv     (§3.3 metric robustness)
data/eval/diag_alttargets_n100.csv           (§3.3 + Appendix B target robustness)
data/eval/diag_alttargets_n100_summary.csv   (§3.3 + Appendix B target robustness)
data/eval/stage5_comparison_n100.csv
data/tables/table1_mechanism_summary.csv
data/tables/table2_gears_matched_n.csv
data/tables/appendix_a_threshold_sensitivity.csv
data/figure_inputs/bootstrap_ci_summary.csv
data/figure_inputs/figure1_panel_summary.csv

Plus a manifest JSON recording inputs, outputs, configuration, git commit,
row counts, and timestamps.

Validation
----------
`--validate` compares every output against the corresponding file in
precomputed/ and reports any differences beyond a configurable tolerance.
Bootstrap CIs are compared with looser tolerance because of Monte Carlo noise.

Usage
-----
    python scripts/recompute_diagnostics.py
    python scripts/recompute_diagnostics.py --predictions-dir custom/path
    python scripts/recompute_diagnostics.py --validate --precomputed-dir precomputed
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import warnings
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

warnings.filterwarnings("ignore", category=RuntimeWarning)  # pearsonr/spearmanr on constant slices

# ── Configuration constants (Methods §2.5) ───────────────────────────────────
DRIVER_THRESH = 0.10            # threshold for driver_count_above_0.10
LOO_MAX_THRESH = 0.15           # operational single-driver threshold
LOW_LOO_THRESH = 0.10           # below which "distributed" is considered
RECURRENCE_FRAC = 0.30          # mode-driver recurrence threshold for single-driver
APPENDIX_THRESHOLDS = [0.125, 0.150, 0.175]

CPA_SEEDS = list(range(1, 101))
GEARS_SEEDS = list(range(1, 8))
CONDITIONS = [("K562", "random"), ("K562", "stratified"),
              ("RPE1", "random"), ("RPE1", "stratified")]


# ── Helpers ──────────────────────────────────────────────────────────────────
def severity_h5ad_path(predictions_dir: Path, model: str,
                       cell_type: str, split_type: str, seed: int) -> Path:
    fname = f"{model}_{cell_type}_{split_type}_seed{seed}.severity.h5ad"
    return predictions_dir / fname


def load_severity_detail(path: Path) -> dict | None:
    """Read one severity-detail h5ad; return None if missing."""
    if not path.exists():
        return None
    import anndata as ad
    obs = ad.read_h5ad(path).obs
    return {
        "mag": obs["predicted_mean_abs_delta"].to_numpy(float),
        "lev": obs["leverage_score"].to_numpy(float),
        "perts": (obs["perturbation_target"].to_numpy()
                  if "perturbation_target" in obs.columns else None),
    }


def loo_scores(mag: np.ndarray, lev: np.ndarray,
               perts: np.ndarray | None) -> dict:
    """LOO sensitivity scores for one holdout (Methods §2.5)."""
    raw = pearsonr(mag, lev)[0]
    n = len(mag)
    absdelta = np.empty(n)
    for i in range(n):
        keep = np.arange(n) != i
        absdelta[i] = abs(raw - pearsonr(mag[keep], lev[keep])[0])
    order_idx = np.argsort(absdelta)[::-1]
    loo_max = float(absdelta[order_idx[0]])
    top3 = float(absdelta[order_idx[:3]].sum())
    driver_count = int((absdelta > DRIVER_THRESH).sum())
    top1_frac = float(loo_max / top3) if top3 > 0 else float("nan")
    out = {
        "severity_pearson_raw": round(float(raw), 4),
        "LOO_max_abs_delta": round(loo_max, 4),
        "LOO_top3_sum_abs_delta": round(top3, 4),
        "driver_count_above_0.10": driver_count,
        "top1_fraction": round(top1_frac, 4),
    }
    if perts is not None:
        out["top_loo_driver"] = str(np.asarray(perts)[order_idx[0]])
    return out


def panel_thresholds(severity_ref_path: Path) -> tuple[float, float, int]:
    """(p5, p95, n_reference_genes) from the full cell-type-specific reference."""
    sev = pd.read_csv(severity_ref_path)
    lev = sev["leverage_score"].dropna().to_numpy(float)
    return (float(np.percentile(lev, 5)), float(np.percentile(lev, 95)), int(len(lev)))


def _mad(x: np.ndarray) -> float:
    x = np.asarray(x, float)
    return float(np.median(np.abs(x - np.median(x))))


def _pct_change(raw_v: float, wins_v: float) -> float:
    if raw_v == 0:
        return float("nan")
    return (wins_v - raw_v) / abs(raw_v) * 100.0


def mechanism_label_and_mode(slice_df: pd.DataFrame) -> tuple[str, str | None, int, int, float]:
    """Apply the three-label rule to a (condition × seeds) slice.
    Returns (label, mode_driver, mode_count, n_high_loo, loo_max_median)."""
    loo_med = float(slice_df["LOO_max_abs_delta"].median())
    high = slice_df[slice_df["LOO_max_abs_delta"] > LOO_MAX_THRESH]
    n_high = len(high)
    if n_high == 0 or "top_loo_driver" not in slice_df.columns:
        mode_driver, mode_count = None, 0
    else:
        c = Counter(high["top_loo_driver"].dropna())
        mode_driver, mode_count = (c.most_common(1)[0] if c else (None, 0))
    mode_frac = (mode_count / n_high) if n_high > 0 else 0.0
    if loo_med > LOO_MAX_THRESH and mode_frac > RECURRENCE_FRAC:
        label = "single-driver"
    elif loo_med < LOW_LOO_THRESH and mode_frac <= RECURRENCE_FRAC:
        label = "distributed"
    else:
        label = "mixed"
    return label, mode_driver, mode_count, n_high, loo_med


def bootstrap_ci_mean(x: np.ndarray, n_resamples: int = 10_000,
                      ci: float = 0.95, seed: int = 0) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    boots = np.empty(n_resamples)
    for i in range(n_resamples):
        boots[i] = float(rng.choice(x, size=len(x), replace=True).mean())
    alpha = (1 - ci) / 2
    return (float(np.percentile(boots, 100 * alpha)),
            float(np.percentile(boots, 100 * (1 - alpha))))


def _git_sha(root: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unknown"


# ── Stage computations ───────────────────────────────────────────────────────
def compute_loo_panel(predictions_dir: Path) -> tuple[pd.DataFrame, list[str]]:
    """Pass 1: per-seed LOO scores for CPA n=100 + GEARS n=7."""
    rows, missing = [], []
    for model, seeds in [("CPA", CPA_SEEDS), ("GEARS", GEARS_SEEDS)]:
        for cell, split in CONDITIONS:
            for seed in seeds:
                d = load_severity_detail(severity_h5ad_path(
                    predictions_dir, model, cell, split, seed))
                if d is None:
                    missing.append(f"{model} {cell} {split} seed{seed}")
                    continue
                rows.append({"model_id": model, "cell_type": cell,
                             "split_type": split, "seed": seed,
                             **loo_scores(d["mag"], d["lev"], d["perts"])})
    return pd.DataFrame(rows), missing


def compute_winsorise_panel(predictions_dir: Path,
                            refs: dict[str, tuple[float, float, int]]
                            ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Pass 2: per-seed CPA winsorisation + 4-row aggregate summary."""
    rows = []
    for cell, split in CONDITIONS:
        p5, p95, _ = refs[cell]
        for seed in CPA_SEEDS:
            d = load_severity_detail(severity_h5ad_path(
                predictions_dir, "CPA", cell, split, seed))
            if d is None:
                continue
            lev_c = np.clip(d["lev"], p5, p95)
            rows.append({
                "condition": f"{cell} {split}",
                "cell_type": cell, "split_type": split, "seed": seed,
                "r_raw": float(pearsonr(d["mag"], d["lev"])[0]),
                "r_wins": float(pearsonr(d["mag"], lev_c)[0]),
                "n_perts_capped": int(((d["lev"] < p5) | (d["lev"] > p95)).sum()),
            })
    per_seed = pd.DataFrame(rows)

    summ_rows = []
    for cell, split in CONDITIONS:
        g = per_seed[(per_seed.cell_type == cell) & (per_seed.split_type == split)]
        if g.empty:
            continue
        raw = g["r_raw"].to_numpy(float)
        wins = g["r_wins"].to_numpy(float)
        rr, wr = float(raw.max() - raw.min()), float(wins.max() - wins.min())
        ram, wam = _mad(raw), _mad(wins)
        rw = float(np.percentile(raw, 95) - np.percentile(raw, 5))
        ww = float(np.percentile(wins, 95) - np.percentile(wins, 5))
        summ_rows.append({
            "condition": f"{cell} {split}",
            "n_present": int(len(g)),
            "raw_median": round(float(np.median(raw)), 4),
            "wins_median": round(float(np.median(wins)), 4),
            "raw_range": round(rr, 4),
            "wins_range": round(wr, 4),
            "range_pct_change": round(_pct_change(rr, wr), 1),
            "raw_MAD": round(ram, 4),
            "wins_MAD": round(wam, 4),
            "MAD_pct_change": round(_pct_change(ram, wam), 1),
            "raw_p5_p95_width": round(rw, 4),
            "wins_p5_p95_width": round(ww, 4),
            "p5_p95_pct_change": round(_pct_change(rw, ww), 1),
            "raw_sign_flips": int((raw < 0).sum()),
            "wins_sign_flips": int((wins < 0).sum()),
            "median_n_perts_capped": float(np.median(g["n_perts_capped"].to_numpy(float))),
        })
    return per_seed, pd.DataFrame(summ_rows)


def compute_stage5(loo_df: pd.DataFrame) -> pd.DataFrame:
    """Pass 3: per-condition CPA mechanism summary (Table 1 source)."""
    cpa = loo_df[loo_df.model_id == "CPA"]
    rows = []
    for cell, split in CONDITIONS:
        g = cpa[(cpa.cell_type == cell) & (cpa.split_type == split)]
        if g.empty:
            continue
        sev = g["severity_pearson_raw"].to_numpy(float)
        label, mode_drv, mode_cnt, n_high, _ = mechanism_label_and_mode(g)
        rows.append({
            "condition": f"{cell} {split}",
            "n": int(len(g)),
            "median": round(float(np.median(sev)), 4),
            "MAD": round(_mad(sev), 4),
            "IQR": round(float(np.percentile(sev, 75) - np.percentile(sev, 25)), 4),
            "p5_p95_width": round(float(np.percentile(sev, 95) - np.percentile(sev, 5)), 4),
            "sign_flips": int((sev < 0).sum()),
            "LOO_max_median": round(float(g.LOO_max_abs_delta.median()), 4),
            "top1_frac_median": round(float(g.top1_fraction.median()), 4),
            "n_LOOmax_gt_0.15": int((g.LOO_max_abs_delta > LOO_MAX_THRESH).sum()),
            "n_high_LOO_seeds": n_high,
            "top_driver_mode": mode_drv,
            "top_driver_mode_count_among_high_LOO": mode_cnt,
        })
    return pd.DataFrame(rows)


def build_table1(stage5: pd.DataFrame, wins_summary: pd.DataFrame) -> pd.DataFrame:
    """Table 1: 8-column mechanism summary, formatted for display.

    Columns match the manuscript Table 1 cell-by-cell: Raw median, Sign flips,
    LOO_max median, tau_1 median, Top driver (recurrence in high-LOO seeds),
    and winsorisation range percent change.
    """
    ORDER = [f"{c} {s}" for c, s in CONDITIONS]
    m = stage5.merge(wins_summary[["condition", "range_pct_change"]],
                     on="condition", how="left")
    m["condition"] = pd.Categorical(m["condition"], categories=ORDER, ordered=True)
    m = m.sort_values("condition").reset_index(drop=True)

    def fmt_sign(x): return f"{x:+.3f}"
    def fmt_pct(x): return f"{x:+.0f}%"
    return pd.DataFrame({
        "Cell type":   [c.split()[0] for c in m["condition"]],
        "Split":       [c.split()[1] for c in m["condition"]],
        "Raw median":  [fmt_sign(v) for v in m["median"]],
        "Sign flips":  [f"{int(v)}/100" for v in m["sign_flips"]],
        "LOO_max median": [f"{v:.3f}" for v in m["LOO_max_median"]],
        "tau_1 median":  [f"{v:.3f}" for v in m["top1_frac_median"]],
        "Top driver in high-LOO seeds": [
            f"{drv} ({int(cnt)}/{int(tot)})"
            for drv, cnt, tot in zip(m["top_driver_mode"],
                                     m["top_driver_mode_count_among_high_LOO"],
                                     m["n_high_LOO_seeds"])
        ],
        "Wins. range Delta": [fmt_pct(v) for v in m["range_pct_change"]],
    })


def build_table2(loo_df: pd.DataFrame) -> pd.DataFrame:
    """Table 2: GEARS vs CPA matched-n verdict."""
    cpa_all = loo_df[loo_df.model_id == "CPA"]
    gears = loo_df[loo_df.model_id == "GEARS"]

    def render_cell(label, driver, loo_med):
        drv_str = driver if driver is not None else "-"
        return f"{label} ({drv_str}; {loo_med:.3f})"

    rows = []
    for cell, split in CONDITIONS:
        gs = gears[(gears.cell_type == cell) & (gears.split_type == split)]
        cm = cpa_all[(cpa_all.cell_type == cell) & (cpa_all.split_type == split)
                     & (cpa_all.seed.between(1, 7))]
        cn = cpa_all[(cpa_all.cell_type == cell) & (cpa_all.split_type == split)]

        g_lbl, g_drv, _, _, g_loo = mechanism_label_and_mode(gs)
        m_lbl, m_drv, _, _, m_loo = mechanism_label_and_mode(cm)
        n_lbl, n_drv, _, _, n_loo = mechanism_label_and_mode(cn)
        verdict = "match" if g_lbl == m_lbl else "differ"
        rows.append({
            "Condition":         f"{cell} {split}",
            "GEARS (n=7)":       render_cell(g_lbl, g_drv, g_loo),
            "CPA matched (n=7)": render_cell(m_lbl, m_drv, m_loo),
            "Verdict":           verdict,
            "CPA n=100":         render_cell(n_lbl, n_drv, n_loo),
        })
    return pd.DataFrame(rows)


def build_appendix_a(loo_df: pd.DataFrame) -> pd.DataFrame:
    """Appendix A: threshold sensitivity counts (CPA only)."""
    cpa = loo_df[loo_df.model_id == "CPA"]
    rows = []
    for thr in APPENDIX_THRESHOLDS:
        row = {"Threshold": f"LOO_max > {thr:.3f}"}
        for cell, split in CONDITIONS:
            sub = cpa[(cpa.cell_type == cell) & (cpa.split_type == split)]
            row[f"{cell} {split}"] = f"{int((sub.LOO_max_abs_delta > thr).sum())}/{len(sub)}"
        rows.append(row)
    return pd.DataFrame(rows)


def build_bootstrap_ci(loo_df: pd.DataFrame, n_resamples: int, seed: int) -> pd.DataFrame:
    cpa = loo_df[loo_df.model_id == "CPA"]
    rows = []
    for cell, split in CONDITIONS:
        g = cpa[(cpa.cell_type == cell) & (cpa.split_type == split)]
        vals = g["severity_pearson_raw"].to_numpy(float)
        if len(vals) == 0:
            continue
        lo, hi = bootstrap_ci_mean(vals, n_resamples=n_resamples,
                                   ci=0.95, seed=seed)
        rows.append({
            "condition": f"{cell} {split}",
            "n": int(len(vals)),
            "mean": round(float(np.mean(vals)), 4),
            "ci95_low": round(lo, 4),
            "ci95_high": round(hi, 4),
        })
    return pd.DataFrame(rows)


def compute_spearman_panel(predictions_dir: Path
                           ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Pass 3: per-seed CPA Pearson + Spearman + 4-row aggregate summary.

    Parallel to compute_winsorise_panel: same per-seed shape, same summary
    metrics, but the contrast is metric choice (Pearson vs rank-based Spearman)
    rather than target choice (raw vs winsorised leverage). Used to test whether
    severity-correlation instability is anchored on Pearson's tail sensitivity
    or persists under rank-based scoring (Manuscript §3.3 metric-robustness
    paragraph).
    """
    rows = []
    for cell, split in CONDITIONS:
        for seed in CPA_SEEDS:
            d = load_severity_detail(severity_h5ad_path(
                predictions_dir, "CPA", cell, split, seed))
            if d is None:
                continue
            rows.append({
                "condition": f"{cell} {split}",
                "cell_type": cell, "split_type": split, "seed": seed,
                "r_pearson": float(pearsonr(d["mag"], d["lev"])[0]),
                "r_spearman": float(spearmanr(d["mag"], d["lev"])[0]),
            })
    per_seed = pd.DataFrame(rows)

    summ_rows = []
    for cell, split in CONDITIONS:
        g = per_seed[(per_seed.cell_type == cell) & (per_seed.split_type == split)]
        if g.empty:
            continue
        pear = g["r_pearson"].to_numpy(float)
        spea = g["r_spearman"].to_numpy(float)
        pr, sr = float(pear.max() - pear.min()), float(spea.max() - spea.min())
        pam, sam = _mad(pear), _mad(spea)
        pw = float(np.percentile(pear, 95) - np.percentile(pear, 5))
        sw = float(np.percentile(spea, 95) - np.percentile(spea, 5))
        summ_rows.append({
            "condition": f"{cell} {split}",
            "n_present": int(len(g)),
            "pearson_median": round(float(np.median(pear)), 4),
            "spearman_median": round(float(np.median(spea)), 4),
            "pearson_range": round(pr, 4),
            "spearman_range": round(sr, 4),
            "range_pct_change": round(_pct_change(pr, sr), 1),
            "pearson_MAD": round(pam, 4),
            "spearman_MAD": round(sam, 4),
            "MAD_pct_change": round(_pct_change(pam, sam), 1),
            "pearson_p5_p95_width": round(pw, 4),
            "spearman_p5_p95_width": round(sw, 4),
            "p5_p95_pct_change": round(_pct_change(pw, sw), 1),
            "pearson_sign_flips": int((pear < 0).sum()),
            "spearman_sign_flips": int((spea < 0).sum()),
        })
    return per_seed, pd.DataFrame(summ_rows)


def compute_alttarget_panel(predictions_dir: Path, severity_refs_dir: Path
                            ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Robustness pass (alternative targets): per-seed CPA correlations against
    deg_count, log1p(deg_count), and |knockdown_pct|, alongside the canonical
    leverage target.

    Joins each held-out perturbation_target to the cell-type-specific severity
    reference. Each holdout records the usable-N for each alternative target,
    plus the count of knockdown_pct = -1.0 sentinel/flag rows in the Replogle
    reference (interpretation of the -1 value is not relied on here).

    Manuscript-facing: leverage / DEG count / log(1+DEG) columns underpin
    the §3.3 alt-target paragraph and the Appendix B sign-flip table.
    Knockdown columns are retained as exploratory output only and are not
    reported in the manuscript; the Replogle reference contains a large
    point mass at the -1 sentinel/flag value, which makes |knockdown|
    correlations hard to interpret without a defensible filtering rule.
    """
    refs = {cell: pd.read_csv(severity_refs_dir / f"replogle_{cell}_severity.csv")
                  .set_index("perturbation_target")
            for cell in ("K562", "RPE1")}
    target_names = ["leverage", "deg", "logdeg", "abskd"]

    rows = []
    for cell, split in CONDITIONS:
        ref = refs[cell]
        for seed in CPA_SEEDS:
            d = load_severity_detail(severity_h5ad_path(
                predictions_dir, "CPA", cell, split, seed))
            if d is None or d["perts"] is None:
                continue
            mag = np.asarray(d["mag"], float)
            j = pd.DataFrame({
                "perturbation_target": d["perts"],
                "mag": mag,
                "leverage": d["lev"],
            }).merge(ref[["deg_count", "knockdown_pct"]],
                     left_on="perturbation_target", right_index=True, how="left")
            row = {
                "condition": f"{cell} {split}",
                "cell_type": cell, "split_type": split, "seed": seed,
                "n_perts": int(len(j)),
                "n_kd_sentinel": int((j["knockdown_pct"] == -1.0).sum()),
                "n_deg_missing": int(j["deg_count"].isna().sum()),
                "n_kd_missing": int(j["knockdown_pct"].isna().sum()),
            }
            target_vals = {
                "leverage": j["leverage"].to_numpy(float),
                "deg": j["deg_count"].to_numpy(float),
                "logdeg": np.log1p(j["deg_count"].to_numpy(float)),
                "abskd": j["knockdown_pct"].abs().to_numpy(float),
            }
            mag_arr = j["mag"].to_numpy(float)
            for tname, tv in target_vals.items():
                ok = np.isfinite(tv) & np.isfinite(mag_arr)
                n_ok = int(ok.sum())
                row[f"n_{tname}"] = n_ok
                if n_ok < 3 or np.std(tv[ok]) == 0 or np.std(mag_arr[ok]) == 0:
                    row[f"r_pearson_{tname}"] = float("nan")
                    row[f"r_spearman_{tname}"] = float("nan")
                else:
                    row[f"r_pearson_{tname}"] = float(pearsonr(mag_arr[ok], tv[ok])[0])
                    row[f"r_spearman_{tname}"] = float(spearmanr(mag_arr[ok], tv[ok])[0])
            rows.append(row)
    per_seed = pd.DataFrame(rows)

    summ_rows = []
    for cell, split in CONDITIONS:
        g = per_seed[(per_seed.cell_type == cell) & (per_seed.split_type == split)]
        if g.empty:
            continue
        row = {
            "condition": f"{cell} {split}",
            "n_seeds": int(len(g)),
            "median_n_perts": float(np.median(g["n_perts"])),
            "median_n_kd_sentinel": float(np.median(g["n_kd_sentinel"])),
            "median_n_kd_missing": float(np.median(g["n_kd_missing"])),
            "median_n_deg_missing": float(np.median(g["n_deg_missing"])),
        }
        for t in target_names:
            for metric in ("pearson", "spearman"):
                col = f"r_{metric}_{t}"
                v = g[col].dropna().to_numpy(float)
                if len(v) == 0:
                    continue
                row[f"{metric}_{t}_median"] = round(float(np.median(v)), 4)
                row[f"{metric}_{t}_MAD"] = round(_mad(v), 4)
                row[f"{metric}_{t}_sign_flips"] = int((v < 0).sum())
                row[f"{metric}_{t}_n_seeds_ok"] = int(len(v))
        summ_rows.append(row)
    return per_seed, pd.DataFrame(summ_rows)


def build_figure1_panel_summary(loo_df: pd.DataFrame,
                                wins_summary: pd.DataFrame) -> pd.DataFrame:
    """Per-condition values plotted in Figure 1 Panel A + B."""
    cpa = loo_df[loo_df.model_id == "CPA"]
    rows = []
    for cell, split in CONDITIONS:
        g = cpa[(cpa.cell_type == cell) & (cpa.split_type == split)]
        if g.empty:
            continue
        sev = g["severity_pearson_raw"].to_numpy(float)
        rows.append({
            "condition": f"{cell} {split}",
            "median": round(float(np.median(sev)), 4),
            "mad": round(_mad(sev), 4),
            "sign_flips": int((sev < 0).sum()),
            "raw_range": round(float(sev.max() - sev.min()), 4),
        })
    panel_a = pd.DataFrame(rows)
    panel_b = wins_summary[["condition", "raw_range", "wins_range",
                            "range_pct_change", "raw_MAD", "wins_MAD",
                            "raw_p5_p95_width", "wins_p5_p95_width"]]
    # Use pandas default merge suffixes (_x, _y) to match the precomputed file's
    # column structure: raw_range_x (from Panel A), raw_range_y (from Panel B).
    return panel_a.merge(panel_b, on="condition", how="left")


# ── Validation ───────────────────────────────────────────────────────────────
def validate_against(out_files: dict[str, Path], precomputed_dir: Path,
                     bootstrap_tolerance: float = 0.01) -> int:
    """Compare each output CSV against precomputed/ and report differences.
    Returns 0 on full match, 2 on any non-tolerated difference."""
    failures = 0
    for name, out_path in out_files.items():
        # Map data/ output paths to precomputed/ paths
        rel = str(out_path).split("data/", 1)[-1]
        ref_path = precomputed_dir / rel
        print(f"\n[validate] {out_path.name}")
        print(f"  out: {out_path}")
        print(f"  ref: {ref_path}")
        if not ref_path.exists():
            print(f"  SKIP — no reference at {ref_path}")
            continue
        if not out_path.exists():
            print(f"  FAIL — generated output missing")
            failures += 1
            continue

        out_df = pd.read_csv(out_path)
        ref_df = pd.read_csv(ref_path)

        if out_df.shape != ref_df.shape:
            print(f"  FAIL — shape differs: out {out_df.shape} vs ref {ref_df.shape}")
            failures += 1
            continue

        is_bootstrap = "bootstrap_ci_summary" in str(out_path)
        tol = bootstrap_tolerance if is_bootstrap else 1e-3
        diffs = 0
        for col in out_df.columns:
            if out_df[col].dtype.kind in "iuf":
                d = (out_df[col].astype(float) - ref_df[col].astype(float)).abs()
                col_diffs = int((d > tol).sum())
                if col_diffs:
                    print(f"  WARN — {col}: {col_diffs} rows differ beyond tol={tol}")
                    diffs += col_diffs
            else:
                col_diffs = int((out_df[col].astype(str) != ref_df[col].astype(str)).sum())
                if col_diffs:
                    print(f"  WARN — {col}: {col_diffs} string rows differ")
                    diffs += col_diffs
        if diffs == 0:
            print(f"  OK — matches precomputed reference{' (within MC tol)' if is_bootstrap else ''}")
        else:
            failures += 1
    return 0 if failures == 0 else 2


# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions-dir", type=Path,
                        default=Path("data/predictions/severity_details"))
    parser.add_argument("--severity-refs-dir", type=Path,
                        default=Path("data/severity_refs"))
    parser.add_argument("--out-eval-dir", type=Path, default=Path("data/eval"))
    parser.add_argument("--out-tables-dir", type=Path, default=Path("data/tables"))
    parser.add_argument("--out-figure-inputs-dir", type=Path,
                        default=Path("data/figure_inputs"))
    parser.add_argument("--manifest", type=Path,
                        default=Path("data/manifest_recompute_diagnostics.json"))
    parser.add_argument("--bootstrap-resamples", type=int, default=10_000)
    parser.add_argument("--bootstrap-seed", type=int, default=0)
    parser.add_argument("--validate", action="store_true",
                        help="After computing outputs, compare against precomputed/")
    parser.add_argument("--precomputed-dir", type=Path,
                        default=Path("precomputed"),
                        help="Reference directory for --validate")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    started = datetime.now(timezone.utc)

    # ── Check inputs ─────────────────────────────────────────────────────────
    if not args.predictions_dir.exists() or not any(args.predictions_dir.iterdir()):
        print(f"ERROR: predictions directory is missing or empty: {args.predictions_dir}",
              file=sys.stderr)
        print("Run scripts/train_cpa_panel.py and scripts/train_gears_panel.py first,",
              file=sys.stderr)
        print("then scripts/evaluate_severity_panel.py to produce severity-detail h5ads.",
              file=sys.stderr)
        return 1
    if not args.severity_refs_dir.exists():
        print(f"ERROR: severity references directory missing: {args.severity_refs_dir}",
              file=sys.stderr)
        return 1

    # ── Compute ──────────────────────────────────────────────────────────────
    print("=== Pass 1: LOO sensitivity panel ===", flush=True)
    loo_df, missing = compute_loo_panel(args.predictions_dir)
    print(f"  {len(loo_df)} rows computed; {len(missing)} h5ads missing")

    print("\n=== Pass 2: winsorisation panel ===", flush=True)
    refs = {}
    for cell in ("K562", "RPE1"):
        path = args.severity_refs_dir / f"replogle_{cell}_severity.csv"
        if not path.exists():
            print(f"ERROR: missing {path}", file=sys.stderr)
            return 1
        refs[cell] = panel_thresholds(path)
        print(f"  {cell}: p5={refs[cell][0]:.4f}  p95={refs[cell][1]:.4f}  "
              f"n_reference_genes={refs[cell][2]}")
    wins_df, wins_summary = compute_winsorise_panel(args.predictions_dir, refs)
    print(f"  {len(wins_df)} per-seed rows; {len(wins_summary)} summary rows")

    print("\n=== Robustness pass: Pearson vs Spearman (CPA n=100) ===", flush=True)
    spear_df, spear_summary = compute_spearman_panel(args.predictions_dir)
    print(f"  {len(spear_df)} per-seed rows; {len(spear_summary)} summary rows")
    print(spear_summary[["condition", "pearson_median", "spearman_median",
                          "pearson_sign_flips", "spearman_sign_flips",
                          "MAD_pct_change", "range_pct_change"]].to_string(index=False))

    print("\n=== Robustness pass: alternative severity targets (CPA n=100) ===", flush=True)
    alt_df, alt_summary = compute_alttarget_panel(args.predictions_dir,
                                                  args.severity_refs_dir)
    print(f"  {len(alt_df)} per-seed rows; {len(alt_summary)} summary rows")
    cols = ["condition", "median_n_perts", "median_n_kd_sentinel",
            "pearson_leverage_median", "spearman_leverage_median",
            "pearson_deg_median",      "spearman_deg_median",
            "pearson_logdeg_median",   "spearman_logdeg_median",
            "pearson_abskd_median",    "spearman_abskd_median"]
    print(alt_summary[cols].to_string(index=False))
    print("\n  Sign-flip counts per target (Pearson | Spearman):")
    for _, r in alt_summary.iterrows():
        print(f"    {r['condition']:18s}  "
              f"lev: {int(r['pearson_leverage_sign_flips'])}/{int(r['pearson_leverage_n_seeds_ok'])} | "
              f"{int(r['spearman_leverage_sign_flips'])}/{int(r['spearman_leverage_n_seeds_ok'])}   "
              f"deg: {int(r['pearson_deg_sign_flips'])}/{int(r['pearson_deg_n_seeds_ok'])} | "
              f"{int(r['spearman_deg_sign_flips'])}/{int(r['spearman_deg_n_seeds_ok'])}   "
              f"logdeg: {int(r['pearson_logdeg_sign_flips'])}/{int(r['pearson_logdeg_n_seeds_ok'])} | "
              f"{int(r['spearman_logdeg_sign_flips'])}/{int(r['spearman_logdeg_n_seeds_ok'])}   "
              f"abskd: {int(r['pearson_abskd_sign_flips'])}/{int(r['pearson_abskd_n_seeds_ok'])} | "
              f"{int(r['spearman_abskd_sign_flips'])}/{int(r['spearman_abskd_n_seeds_ok'])}")

    print("\n=== Pass 3: stage5 per-condition aggregates ===", flush=True)
    stage5 = compute_stage5(loo_df)
    print(stage5.to_string(index=False))

    print("\n=== Build Table 1 ===", flush=True)
    table1 = build_table1(stage5, wins_summary)
    print(table1.to_string(index=False))

    print("\n=== Build Table 2 ===", flush=True)
    table2 = build_table2(loo_df)
    print(table2.to_string(index=False))

    print("\n=== Build Appendix A ===", flush=True)
    appendix_a = build_appendix_a(loo_df)
    print(appendix_a.to_string(index=False))

    print("\n=== Build bootstrap CI summary ===", flush=True)
    boot_ci = build_bootstrap_ci(loo_df, args.bootstrap_resamples, args.bootstrap_seed)
    print(boot_ci.to_string(index=False))

    print("\n=== Build Figure 1 panel summary ===", flush=True)
    fig1 = build_figure1_panel_summary(loo_df, wins_summary)
    print(fig1.to_string(index=False))

    # ── Write outputs ────────────────────────────────────────────────────────
    args.out_eval_dir.mkdir(parents=True, exist_ok=True)
    args.out_tables_dir.mkdir(parents=True, exist_ok=True)
    args.out_figure_inputs_dir.mkdir(parents=True, exist_ok=True)

    # Split LOO output by model: precomputed convention has CPA in
    # diag_loo_sensitivity_n100.csv (400 rows) and GEARS in
    # diag_loo_sensitivity_gears.csv (28 rows) as separate files.
    loo_cpa = loo_df[loo_df.model_id == "CPA"].copy()
    loo_gears = loo_df[loo_df.model_id == "GEARS"].copy()

    out_files = {
        "diag_loo_sensitivity_n100":     args.out_eval_dir / "diag_loo_sensitivity_n100.csv",
        "diag_loo_sensitivity_gears":    args.out_eval_dir / "diag_loo_sensitivity_gears.csv",
        "diag_winsorise_n100":           args.out_eval_dir / "diag_winsorise_n100.csv",
        "diag_winsorise_n100_summary":   args.out_eval_dir / "diag_winsorise_n100_summary.csv",
        "diag_spearman_n100":            args.out_eval_dir / "diag_spearman_n100.csv",
        "diag_spearman_n100_summary":    args.out_eval_dir / "diag_spearman_n100_summary.csv",
        "diag_alttargets_n100":          args.out_eval_dir / "diag_alttargets_n100.csv",
        "diag_alttargets_n100_summary":  args.out_eval_dir / "diag_alttargets_n100_summary.csv",
        "stage5_comparison_n100":        args.out_eval_dir / "stage5_comparison_n100.csv",
        "table1_mechanism_summary":      args.out_tables_dir / "table1_mechanism_summary.csv",
        "table2_gears_matched_n":        args.out_tables_dir / "table2_gears_matched_n.csv",
        "appendix_a_threshold_sensitivity": args.out_tables_dir / "appendix_a_threshold_sensitivity.csv",
        "bootstrap_ci_summary":          args.out_figure_inputs_dir / "bootstrap_ci_summary.csv",
        "figure1_panel_summary":         args.out_figure_inputs_dir / "figure1_panel_summary.csv",
    }
    loo_cpa.to_csv(out_files["diag_loo_sensitivity_n100"], index=False)
    loo_gears.to_csv(out_files["diag_loo_sensitivity_gears"], index=False)
    wins_df.to_csv(out_files["diag_winsorise_n100"], index=False)
    wins_summary.to_csv(out_files["diag_winsorise_n100_summary"], index=False)
    spear_df.to_csv(out_files["diag_spearman_n100"], index=False)
    spear_summary.to_csv(out_files["diag_spearman_n100_summary"], index=False)
    alt_df.to_csv(out_files["diag_alttargets_n100"], index=False)
    alt_summary.to_csv(out_files["diag_alttargets_n100_summary"], index=False)
    stage5.to_csv(out_files["stage5_comparison_n100"], index=False)
    table1.to_csv(out_files["table1_mechanism_summary"], index=False)
    table2.to_csv(out_files["table2_gears_matched_n"], index=False)
    appendix_a.to_csv(out_files["appendix_a_threshold_sensitivity"], index=False)
    boot_ci.to_csv(out_files["bootstrap_ci_summary"], index=False)
    fig1.to_csv(out_files["figure1_panel_summary"], index=False)

    print("\n=== Outputs ===")
    for name, path in out_files.items():
        size = path.stat().st_size if path.exists() else 0
        print(f"  {name:<40s}  {path}  ({size:,} bytes)")

    # ── Manifest ─────────────────────────────────────────────────────────────
    finished = datetime.now(timezone.utc)
    manifest = {
        "command": " ".join(sys.argv),
        "git_commit": _git_sha(repo_root),
        "started_utc": started.isoformat(),
        "finished_utc": finished.isoformat(),
        "duration_seconds": (finished - started).total_seconds(),
        "inputs": {
            "predictions_dir": str(args.predictions_dir),
            "severity_refs_dir": str(args.severity_refs_dir),
            "panel_thresholds": {
                cell: {"p5": p5, "p95": p95, "n_reference_genes": n}
                for cell, (p5, p95, n) in refs.items()
            },
            "missing_h5ads": missing,
        },
        "configuration": {
            "DRIVER_THRESH": DRIVER_THRESH,
            "LOO_MAX_THRESH": LOO_MAX_THRESH,
            "LOW_LOO_THRESH": LOW_LOO_THRESH,
            "RECURRENCE_FRAC": RECURRENCE_FRAC,
            "APPENDIX_THRESHOLDS": APPENDIX_THRESHOLDS,
            "bootstrap_resamples": args.bootstrap_resamples,
            "bootstrap_seed": args.bootstrap_seed,
        },
        "outputs": {
            name: {"path": str(path), "size_bytes": path.stat().st_size,
                   "row_count": int(pd.read_csv(path).shape[0])}
            for name, path in out_files.items() if path.exists()
        },
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2))
    print(f"\nmanifest: {args.manifest}")

    # ── Validate ─────────────────────────────────────────────────────────────
    if args.validate:
        print("\n=== Validation against precomputed/ ===")
        return validate_against(out_files, args.precomputed_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
