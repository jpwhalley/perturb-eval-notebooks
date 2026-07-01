# Training Provenance

This document records how the per-seed prediction h5ads that underpin the
manuscript were generated. It is intentionally separate from `README.md`
because the public repository's default reproducibility path **does not
retrain the perturbation models**; it reproduces the manuscript from shipped
evaluation artefacts and validates the evaluation layer against stored
prediction outputs.

This file explains where the training step lives, why it is documented rather
than packaged as the default workflow, and what guarantees the public
reproducibility layers do provide.

---

## What the public repository ships

The public repository is organised as three layers, each with a defined
reproducibility commitment:

| Layer | Default-path runtime | Validation |
|---|---|---|
| `precomputed/` (canonical evaluation CSVs and table outputs) | n/a — shipped | Used as the reference for every other layer |
| `scripts/recompute_diagnostics.py` (severity-detail h5ads → evaluation CSVs) | ~10 s on CPU | `--validate` compares all 14 outputs against `precomputed/`; verified 14/14 match |
| `scripts/evaluate_severity_panel.py` (raw prediction h5ads → severity-detail h5ads) | ~1 min for 428 files on CPU | `--validate-existing` compares freshly-derived detail h5ads against the existing 428; verified 0 failed, 0 skipped |

These three layers reproduce every manuscript figure, table, and in-text
statistic from the shipped artefacts in about two minutes, without external
data, model dependencies, or training time.

---

## What the public repository does not ship

The substantive CPA and GEARS training loops that produced the 428 raw
prediction h5ads are not included as runnable code in this repository. The
`P01_train_cpa.ipynb` and `P02_train_gears.ipynb` notebooks document the
training contract (inputs, outputs, expected runtime, hyperparameters) but
delegate the substantive training step to private scripts described below.

This separation is deliberate. Three reasons:

1. **Numerical stability.** CPA and GEARS training are stochastic in ways that
   are not fully controlled by a fixed random seed: floating-point ordering
   inside CUDA/MPS kernels, scvi-tools and torch minor-version drift, and
   batch-ordering differences across hardware can all produce per-seed
   predictions that differ at the 3rd or 4th decimal place from the manuscript
   reference. Shipping a "rerun training from scratch" path would create a
   second reference that does not exactly match the validated 428 prediction
   h5ads, and would erode the byte-exact validation chain documented in
   `MANUSCRIPT_TRACEABILITY.md`.

2. **Dependency weight.** A reviewer running the full retraining pipeline
   would need the `cpa` and `gears` dependency groups
   (`uv sync --group cpa --group gears --group data`), GPU hardware
   (Apple Silicon MPS or CUDA), ~20 GB of disk for atlases plus predictions
   plus logs, and several hours of compute per panel. The default 2-minute
   reproducibility path is a substantially lower commitment than the
   retraining path would require.

3. **Scope discipline.** This manuscript's contribution is an audit of the
   severity-correlation evaluation statistic, not a contribution to CPA or
   GEARS as model architectures. The training code is therefore framed as
   "where the inputs to the eval layer come from" rather than as a
   first-class deliverable of this paper.

---

## Scripts that produced the 428 prediction h5ads

The actual training and prediction was carried out by the following private
scripts in the author's working repository (not part of this public release):

| Script | Role |
|---|---|
| `cpa_20_train.py` | Trains CPA from scratch on the K562 or RPE1 atlas for one (cell type, split type, seed). Forces CPU due to a known scvi-tools float64 BatchNorm crash on MPS. |
| `gears_20_train.py` | Trains GEARS on the K562 or RPE1 atlas for one (cell type, split type, seed). Uses MPS when available, CPU fallback. |
| `tools/expand_cpa_seeds.py` | Orchestrates per-seed CPA training across many seeds in one condition. |
| `gears_rpe1_scout_train.py` | GEARS training variant used for the matched-n cross-architecture scope check. |
| `results/overnight/orchestrate_cpa.sh` | Overnight orchestrator covering CPA seeds 8–30 across all four conditions, sequential, resilient to per-seed failures. |
| `results/overnight/orchestrate_cpa_31_100.sh` | Continuation orchestrator covering CPA seeds 31–100. |
| `results/overnight/orchestrate_gears_rpe1.sh` | GEARS RPE1 matched-n orchestrator (seeds 1–7). |
| `results/overnight/orchestrate_gears_rpe1_67.sh` | GEARS continuation orchestrator. |

Hyperparameters (already documented in Methods §2.2 of the manuscript):

- **CPA**: 30 epochs, batch size 128, Adam optimiser using the CPA 0.8.8 TrainingPlan default learning rate ($5 \times 10^{-4}$; the `--lr` CLI flag in `cpa_20_train.py` is recorded in the run manifest but not forwarded to `model.train`, which uses the TrainingPlan default). The `--training-seed` CLI flag (default 1) is applied only to the validation-split RNG inside the training script; the CPA model initialisation itself uses the library default seed. The input gene set is the top 2,000 highly variable genes plus perturbation target genes. CPU. ≈2 minutes per (cell type, split, seed) on a modern laptop.
- **GEARS**: 30 epochs, default GEARS hyperparameters, fixed training seed 1,
  same input gene set as CPA. MPS preferred, CPU fallback. ≈5–10 minutes per
  (cell type, split, seed).

Output naming convention:
`{model}_{cell_type}_mid_hvg_holdout30_seed{N}_{run_tag}.severity.eval.h5ad`

These files are read by `scripts/evaluate_severity_panel.py`, which projects
them into the standardised `data/predictions/severity_details/` layout used
by the public evaluation layer.

---

## What the public validation chain actually verifies

A reviewer cannot retrain the panels from this repository alone. They can,
however, verify the following with the shipped artefacts:

1. **Evaluation-CSV layer**: `scripts/recompute_diagnostics.py --validate`
   recomputes every manuscript-canonical CSV from the severity-detail h5ads
   and confirms a byte-exact match against `precomputed/` (14/14 outputs;
   bootstrap CI summary within Monte Carlo tolerance).

2. **Severity-detail-layer**:
   `scripts/evaluate_severity_panel.py --validate-existing data/predictions/severity_details`
   re-derives every severity-detail h5ad from raw prediction h5ads and confirms
   each one matches the shipped reference (428 files checked, 0 failed,
   0 skipped). This validates that any prediction h5ad following the documented
   `.obs` contract round-trips correctly through the evaluation layer.

3. **Figure/table layer**: `./run_all.sh --figures` executes the five
   visualisation notebooks against shipped CSVs and regenerates every
   manuscript figure, table, and Appendix B in approximately two minutes,
   without external data or model dependencies.

The validation chain therefore runs:

```
raw prediction h5ads  (privately generated; not shipped, but contract documented)
        ↓ scripts/evaluate_severity_panel.py            (validated round-trip on 428 files)
severity-detail h5ads  (shipped under data/predictions/severity_details/)
        ↓ scripts/recompute_diagnostics.py              (validated byte-exact on 14 outputs)
precomputed/ CSVs      (shipped, canonical manuscript reference)
        ↓ notebooks 01–05                                (deterministic reproduction)
Figure 1, Tables 1–2, Appendix A, Appendix B, in-text §3.3 numbers
```

Every arrow except the first is validated end-to-end and reproducible in the
public repository. The first arrow (training → raw predictions) is documented
here for transparency about provenance.

---

## If you want to retrain anyway

You will need:

- Atlases from Replogle 2022 K562 and RPE1, harmonised via scverse/pertpy
  (D-series notebooks document the acquisition; D01 fetches, D02 builds the
  severity references, D03 produces the deterministic holdout specifications).
- Training environment: `uv sync --group cpa --group gears --group data`,
  plus working PyTorch with MPS (Apple Silicon) or CUDA. CPU-only is
  technically supported but ~10× slower per seed.
- Approximate compute budget for the full manuscript panel: 100 CPA seeds ×
  4 conditions ≈ 13 hours on a single Apple Silicon laptop (CPU-bound);
  7 GEARS seeds × 4 conditions ≈ 2–4 hours on MPS.

Any predictions you produce should be written to
`data/predictions/severity_details/{MODEL}_{cell_type}_{split_type}_seed{N}.severity.h5ad`
with the `.obs` columns documented in
`scripts/evaluate_severity_panel.py`. You can then validate your panel
against the shipped reference with
`scripts/evaluate_severity_panel.py --validate-existing data/predictions/severity_details`
and (within Monte Carlo tolerance for the bootstrap CI) reproduce every
manuscript-facing number via `scripts/recompute_diagnostics.py --validate`.

Per-seed agreement with the shipped predictions at the 3rd–4th decimal place
is the realistic expectation, not byte-exact match.
