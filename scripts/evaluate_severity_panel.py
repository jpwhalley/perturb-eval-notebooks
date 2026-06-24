#!/usr/bin/env python
"""scripts/evaluate_severity_panel.py — adapter layer between raw model
predictions and the standardised severity-detail h5ads consumed by
scripts/recompute_diagnostics.py.

This is Phase 2 of the training-aware reproducibility pipeline. It takes the
model-specific raw prediction artefacts produced by P01 (CPA) and P02 (GEARS)
and emits severity-detail h5ads with the contract Phase 1 expects:

    .obs columns:
        perturbation_target        : gene symbol (one row per held-out perturbation)
        predicted_mean_abs_delta   : per-perturbation mean absolute delta from
                                     control mean (on log1p_normalized scale)
        leverage_score             : per-perturbation severity (from severity ref)
        knockdown_pct, deg_count   : optional metadata from severity ref
        predicted_l2_delta         : per-perturbation L2 of delta
        predicted_max_abs_delta    : per-perturbation max abs delta
        predicted_n_genes_above_thresh : count of |delta| > MAGNITUDE_THRESHOLD
        rank_predicted, rank_leverage  : descending rank of magnitude / leverage

This is the "contract layer". After Phase 2, everything downstream is
model-agnostic.

Input contract
--------------
Raw prediction h5ad at `{predictions-raw-dir}/{model}_{cell_type}_{split_type}_seed{seed}.h5ad`

with .obs columns:
    perturbation_target : gene symbol per cell (or per perturbation for GEARS)
    split               : "test" / "train" (used to filter held-out rows)
    is_control          : bool
and .uns:
    source_dataset_id   : key for the ground-truth atlas (atlas resolution)
    expression_scale    : source scale, e.g. "raw_counts_predicted"
    control_label       : how controls are labelled in this prediction file

The corresponding ground-truth atlas is read from
`{atlases-dir}/{cell_type}_essential.h5ad` (D01 convention) for accurate
control-mean computation. The cell type is inferred from the filename.

The per-cell-type severity reference is read from
`{severity-refs-dir}/replogle_{cell_type}_severity.csv` (D02 convention).

Output contract
---------------
`{out-dir}/{model}_{cell_type}_{split_type}_seed{seed}.severity.h5ad`

The Phase 1 diagnostics script (scripts/recompute_diagnostics.py) consumes
exactly the three columns predicted_mean_abs_delta, leverage_score, and
perturbation_target. The other columns above are written for completeness
and to match the existing a_perturb_AI detail h5ad schema.

Validation
----------
`--validate-existing {dir}` reads each existing severity-detail h5ad in {dir}
and compares predicted_mean_abs_delta + leverage_score against the freshly-
derived values. This is the test that proves Phase 2 is a faithful adapter:
running it on the user's a_perturb_AI severity details (symlinked into the
public repo) should yield zero-or-tiny-tolerance differences.

Usage
-----
    python scripts/evaluate_severity_panel.py --dry-run
    python scripts/evaluate_severity_panel.py --model CPA --cell-type K562 --seed-start 1 --seed-end 2
    python scripts/evaluate_severity_panel.py --validate-existing data/predictions/severity_details
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

MAGNITUDE_THRESHOLD = 0.1


# ── Canonical scale (ported from tools/canonical_scale.py) ──────────────────
def to_canonical_scale(adata: ad.AnnData, source_scale: str) -> np.ndarray:
    """Convert .X to log1p_normalized scale, returning a (n_obs, n_vars) float64 array."""
    import scanpy as sc
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X, dtype=np.float64)
    X = X.astype(np.float64)
    if source_scale in ("raw_counts", "raw_counts_predicted"):
        # Clip negatives to zero before log1p
        X = np.clip(X, 0, None)
        tmp = ad.AnnData(X=X)
        sc.pp.normalize_total(tmp, target_sum=1e4)
        X = np.log1p(tmp.X)
    elif source_scale == "log1p_normalized":
        pass
    elif source_scale == "normalized_linear":
        X = np.log1p(X)
    else:
        raise ValueError(f"Unknown source_scale: {source_scale!r}")
    return np.asarray(X, dtype=np.float64)


# ── Helpers ──────────────────────────────────────────────────────────────────
def parse_pred_filename(path: Path) -> tuple[str, str, str, int] | None:
    """Parse {MODEL}_{CELL}_{SPLIT}_seed{SEED}.h5ad → (model, cell, split, seed)."""
    m = re.match(r"^(CPA|GEARS)_(K562|RPE1)_(random|stratified)_seed(\d+)\.h5ad$",
                 path.name)
    if not m:
        return None
    return m.group(1), m.group(2), m.group(3), int(m.group(4))


def list_raw_predictions(predictions_raw_dir: Path,
                         model: str | None, cell: str | None,
                         split: str | None,
                         seed_range: tuple[int, int] | None) -> list[tuple[Path, str, str, str, int]]:
    """List raw prediction files matching the filter criteria."""
    out = []
    for p in sorted(predictions_raw_dir.glob("*.h5ad")):
        parsed = parse_pred_filename(p)
        if parsed is None:
            continue
        m, c, s, sd = parsed
        if model and model != m:
            continue
        if cell and cell != c:
            continue
        if split and split != s:
            continue
        if seed_range and not (seed_range[0] <= sd <= seed_range[1]):
            continue
        out.append((p, m, c, s, sd))
    return out


def compute_severity_detail(pred_path: Path, atlas_path: Path,
                            severity_ref_path: Path) -> tuple[ad.AnnData, dict]:
    """Compute severity-detail h5ad for a single (model, cell, split, seed)."""
    # ── Load ─────────────────────────────────────────────────────────────────
    pred = ad.read_h5ad(pred_path)
    atlas = ad.read_h5ad(atlas_path)
    severity_df = pd.read_csv(severity_ref_path)

    # ── Filter prediction rows: test + non-control ───────────────────────────
    pred_control_label = pred.uns.get("control_label", "control")
    test_mask = (pred.obs.get("split", pd.Series(index=pred.obs_names)) == "test")
    if "is_control" in pred.obs.columns:
        test_mask = test_mask & (~pred.obs["is_control"].astype(bool))
    pred_filtered = pred[test_mask].copy()
    test_perts = set(pred_filtered.obs["perturbation_target"].unique())

    # ── Filter atlas: controls + test perturbations ──────────────────────────
    atlas_pert_col = "perturbation" if "perturbation" in atlas.obs.columns else "perturbation_target"
    atlas_control_label = "control"  # convention across both atlases
    atlas_pert_values = atlas.obs[atlas_pert_col].values
    ctrl_mask = atlas_pert_values == atlas_control_label
    atlas_ctrl = atlas[ctrl_mask].copy()

    # ── Intersect genes + canonical scale ────────────────────────────────────
    shared_genes = sorted(set(pred_filtered.var_names) & set(atlas.var_names))
    if not shared_genes:
        raise ValueError(f"No shared genes between {pred_path.name} and atlas {atlas_path.name}")
    pred_sub = pred_filtered[:, shared_genes].copy()
    atlas_ctrl_sub = atlas_ctrl[:, shared_genes].copy()

    pred_scale = pred.uns.get("expression_scale", "raw_counts_predicted")
    atlas_scale = atlas.uns.get("expression_scale", "raw_counts")

    pred_X = to_canonical_scale(pred_sub, pred_scale)
    atlas_ctrl_X = to_canonical_scale(atlas_ctrl_sub, atlas_scale)
    control_mean = atlas_ctrl_X.mean(axis=0)

    # ── Per-perturbation predicted magnitude ─────────────────────────────────
    pred_perts = pred_filtered.obs["perturbation_target"].values
    unique_test_perts = sorted(test_perts - {pred_control_label})

    records = []
    for p in unique_test_perts:
        pmask = pred_perts == p
        if pmask.sum() == 0:
            continue
        pred_post = pred_X[pmask].mean(axis=0)
        delta = pred_post - control_mean
        records.append({
            "perturbation_target": p,
            "predicted_mean_abs_delta": float(np.mean(np.abs(delta))),
            "predicted_l2_delta":      float(np.sqrt(np.mean(delta ** 2))),
            "predicted_max_abs_delta": float(np.max(np.abs(delta))),
            "predicted_n_genes_above_thresh": int(np.sum(np.abs(delta) > MAGNITUDE_THRESHOLD)),
        })
    pred_mag_df = pd.DataFrame(records)

    # ── Join with severity reference ─────────────────────────────────────────
    keep_sev_cols = ["perturbation_target", "leverage_score"]
    for extra in ("knockdown_pct", "deg_count", "n_source_rows", "leverage_score_std"):
        if extra in severity_df.columns:
            keep_sev_cols.append(extra)
    joined = pred_mag_df.merge(severity_df[keep_sev_cols], on="perturbation_target", how="inner")
    dropped = sorted(set(pred_mag_df["perturbation_target"]) - set(joined["perturbation_target"]))

    if len(joined) == 0:
        raise ValueError(f"No perturbations survived severity join for {pred_path.name}")

    # ── Ranks ────────────────────────────────────────────────────────────────
    joined["rank_predicted"] = joined["predicted_mean_abs_delta"].rank(ascending=False, method="min").astype(int)
    joined["rank_leverage"] = joined["leverage_score"].rank(ascending=False, method="min").astype(int)

    # ── Build detail h5ad ────────────────────────────────────────────────────
    obs_df = joined.reset_index(drop=True)
    X = joined["predicted_mean_abs_delta"].values.reshape(-1, 1).astype(np.float32)
    var_df = pd.DataFrame(index=["predicted_mean_abs_delta"])
    detail = ad.AnnData(X=X, obs=obs_df, var=var_df)
    detail.uns.update({
        "canonical_scale": "log1p_normalized",
        "control_rule": "all_atlas_controls_mean",
        "predicted_magnitude_definition": "mean_abs_delta",
        "source_prediction": pred_path.name,
        "source_atlas": atlas_path.name,
        "source_severity_reference": severity_ref_path.name,
        "n_perts_joined": int(len(joined)),
        "n_perts_dropped_unjoined": int(len(dropped)),
        "dropped_perturbations": dropped,
    })

    return detail, {
        "n_perts_joined": int(len(joined)),
        "n_perts_dropped_unjoined": int(len(dropped)),
        "dropped_perturbations": dropped,
    }


def _git_sha(root: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unknown"


def validate_existing(existing_dir: Path, predictions_raw_dir: Path,
                      atlases_dir: Path, severity_refs_dir: Path,
                      tol: float = 1e-4) -> int:
    """Compare each existing severity-detail h5ad against a freshly-derived one.
    Returns 0 on full match, 2 on any non-tolerated difference."""
    failures = 0
    checked = 0
    for existing_path in sorted(existing_dir.glob("*.severity.h5ad")):
        parsed = parse_pred_filename(existing_path.with_name(existing_path.name.replace(".severity.h5ad", ".h5ad")))
        if parsed is None:
            continue
        model, cell, split, seed = parsed
        raw_path = predictions_raw_dir / f"{model}_{cell}_{split}_seed{seed}.h5ad"
        if not raw_path.exists():
            continue
        sev_path = severity_refs_dir / f"replogle_{cell}_severity.csv"
        atlas_path = atlases_dir / f"{cell.lower()}_essential.h5ad"
        if not (sev_path.exists() and atlas_path.exists()):
            continue

        checked += 1
        ex_obs = ad.read_h5ad(existing_path).obs
        new_detail, _ = compute_severity_detail(raw_path, atlas_path, sev_path)
        new_obs = new_detail.obs
        # Align on perturbation_target
        m = ex_obs.merge(new_obs[["perturbation_target", "predicted_mean_abs_delta", "leverage_score"]],
                         on="perturbation_target", how="inner", suffixes=("_existing", "_new"))
        d_mag = (m["predicted_mean_abs_delta_existing"].astype(float)
                 - m["predicted_mean_abs_delta_new"].astype(float)).abs()
        d_lev = (m["leverage_score_existing"].astype(float)
                 - m["leverage_score_new"].astype(float)).abs()
        bad_mag = int((d_mag > tol).sum())
        bad_lev = int((d_lev > tol).sum())
        if bad_mag or bad_lev:
            print(f"[validate-existing] FAIL {existing_path.name}: "
                  f"{bad_mag} mag diffs, {bad_lev} lev diffs > {tol}")
            failures += 1
        else:
            print(f"[validate-existing] OK   {existing_path.name}: "
                  f"{len(m)} perts match within tol={tol}")
    print(f"\n=== Checked {checked} files; {failures} failed ===")
    return 0 if failures == 0 else 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions-raw-dir", type=Path,
                        default=Path("data/predictions/raw"),
                        help="Directory of raw model prediction h5ads.")
    parser.add_argument("--atlases-dir", type=Path,
                        default=Path("data/replogle"),
                        help="Directory of ground-truth atlas h5ads "
                             "({k562,rpe1}_essential.h5ad).")
    parser.add_argument("--severity-refs-dir", type=Path,
                        default=Path("data/severity_refs"))
    parser.add_argument("--out-dir", type=Path,
                        default=Path("data/predictions/severity_details"))
    parser.add_argument("--model", choices=["CPA", "GEARS"], default=None,
                        help="Filter: only process this model.")
    parser.add_argument("--cell-type", choices=["K562", "RPE1"], default=None,
                        help="Filter: only process this cell type.")
    parser.add_argument("--split-type", choices=["random", "stratified"], default=None,
                        help="Filter: only process this split type.")
    parser.add_argument("--seed-start", type=int, default=None,
                        help="Filter: lower seed bound (inclusive).")
    parser.add_argument("--seed-end", type=int, default=None,
                        help="Filter: upper seed bound (inclusive).")
    parser.add_argument("--dry-run", action="store_true",
                        help="List matching raw predictions and exit; do not compute.")
    parser.add_argument("--validate-existing", type=Path, default=None,
                        help="Compare existing severity-detail h5ads in DIR against "
                             "freshly-derived ones; report any non-tolerated diffs.")
    parser.add_argument("--manifest", type=Path,
                        default=Path("data/manifest_evaluate_severity.json"))
    parser.add_argument("--skip-existing", action="store_true",
                        help="If the output severity-detail h5ad already exists, skip it.")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    started = datetime.now(timezone.utc)

    # ── Validation mode ──────────────────────────────────────────────────────
    if args.validate_existing is not None:
        return validate_existing(args.validate_existing, args.predictions_raw_dir,
                                 args.atlases_dir, args.severity_refs_dir)

    # ── Discover raw predictions ─────────────────────────────────────────────
    if not args.predictions_raw_dir.exists():
        print(f"ERROR: predictions-raw-dir missing: {args.predictions_raw_dir}", file=sys.stderr)
        return 1
    seed_range = ((args.seed_start or 1), (args.seed_end or 10**9))
    targets = list_raw_predictions(args.predictions_raw_dir, args.model,
                                   args.cell_type, args.split_type,
                                   seed_range if args.seed_start or args.seed_end else None)
    if not targets:
        print(f"No matching raw predictions in {args.predictions_raw_dir}")
        return 1

    print(f"=== {len(targets)} raw predictions to process ===")
    if args.dry_run:
        for p, m, c, s, sd in targets[:20]:
            print(f"  {p.name}  ->  {args.out_dir / (p.stem + '.severity.h5ad')}")
        if len(targets) > 20:
            print(f"  ... and {len(targets) - 20} more")
        return 0

    args.out_dir.mkdir(parents=True, exist_ok=True)

    # ── Process each ─────────────────────────────────────────────────────────
    results = []
    failures = []
    for i, (raw_path, model, cell, split, seed) in enumerate(targets, start=1):
        out_path = args.out_dir / f"{model}_{cell}_{split}_seed{seed}.severity.h5ad"
        if args.skip_existing and out_path.exists():
            print(f"  [{i}/{len(targets)}] SKIP {out_path.name} (already exists)")
            continue

        atlas_path = args.atlases_dir / f"{cell.lower()}_essential.h5ad"
        sev_path = args.severity_refs_dir / f"replogle_{cell}_severity.csv"
        for required, label in [(atlas_path, "atlas"), (sev_path, "severity ref")]:
            if not required.exists():
                msg = f"  [{i}/{len(targets)}] FAIL {raw_path.name}: missing {label} {required}"
                print(msg)
                failures.append({"raw": str(raw_path), "missing": str(required)})
                continue

        try:
            detail, summary = compute_severity_detail(raw_path, atlas_path, sev_path)
            detail.write_h5ad(out_path)
            results.append({
                "raw": str(raw_path),
                "out": str(out_path),
                **summary,
            })
            print(f"  [{i}/{len(targets)}] OK   {out_path.name}  "
                  f"(joined {summary['n_perts_joined']}, dropped {summary['n_perts_dropped_unjoined']})")
        except Exception as e:
            failures.append({"raw": str(raw_path), "error": f"{type(e).__name__}: {e}"})
            print(f"  [{i}/{len(targets)}] FAIL {raw_path.name}: {e}")

    # ── Manifest ─────────────────────────────────────────────────────────────
    finished = datetime.now(timezone.utc)
    manifest = {
        "command": " ".join(sys.argv),
        "git_commit": _git_sha(repo_root),
        "started_utc": started.isoformat(),
        "finished_utc": finished.isoformat(),
        "duration_seconds": (finished - started).total_seconds(),
        "inputs": {
            "predictions_raw_dir": str(args.predictions_raw_dir),
            "atlases_dir": str(args.atlases_dir),
            "severity_refs_dir": str(args.severity_refs_dir),
            "filters": {
                "model": args.model, "cell_type": args.cell_type,
                "split_type": args.split_type,
                "seed_start": args.seed_start, "seed_end": args.seed_end,
            },
        },
        "configuration": {
            "MAGNITUDE_THRESHOLD": MAGNITUDE_THRESHOLD,
            "canonical_scale": "log1p_normalized",
            "control_rule": "all_atlas_controls_mean",
        },
        "results": results,
        "failures": failures,
        "n_processed": len(results),
        "n_failed": len(failures),
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2))
    print(f"\nmanifest: {args.manifest}")
    print(f"=== Done: {len(results)} processed, {len(failures)} failed ===")
    return 0 if not failures else 2


if __name__ == "__main__":
    sys.exit(main())
