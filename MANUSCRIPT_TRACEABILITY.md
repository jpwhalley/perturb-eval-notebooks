# Manuscript–Notebook Traceability Matrix

> **AI-assisted documentation.** This traceability matrix and the cross-references
> between manuscript claims and notebook cells were assembled with the assistance
> of Claude (Anthropic) via Cowork. All scientific analyses, interpretations, and
> editorial decisions were made by the author. The AI tools were used to
> systematically audit the mapping between code outputs and manuscript text.
> Readers can use this document to independently verify any claim in the
> manuscript against the corresponding notebook cell.

This document maps every figure panel, table entry, and in-text statistic in
the manuscript to the notebook and cell that produces it.

**Pipeline:** D-series (acquire) → P-series (analyse) → Numbered (visualise).
See `README.md` for full pipeline documentation.

**Cell-numbering convention.** Cells are 1-indexed including markdown cells, so
"cell 4" of `01_figure1.ipynb` is the fourth cell of that notebook (the data
load + sanity print). Each notebook is small enough (8–10 cells) that this
1-indexed scheme is unambiguous on inspection. The figure/table outputs are
embedded in the executed notebooks for direct inspection without re-running.

**Notebook quick reference.**

| Notebook | Cells of interest | Produces |
|---|---|---|
| `01_figure1` | cell 4 (load), cell 6 (panel renderers), cell 8 (compose + save), cell 10 (bootstrap CIs) | Figure 1 (PNG + PDF); panel summary CSV; bootstrap CI summary |
| `02_table1_mechanism_summary` | cell 4 (load), cell 6 (build), cell 8 (save + display) | Table 1 CSV |
| `03_table2_gears_scope` | cell 4 (load), cell 6 (mechanism label fn), cell 8 (build matched-n) | Table 2 CSV |
| `04_appendix_threshold_sensitivity` | cell 4 (load), cell 6 (compute) | Appendix A CSV |

---

## Script-level provenance

Every manuscript claim has two equally valid provenance paths: the notebook cell tabulated below, and the CLI script that produces the underlying `precomputed/` artefact. The two paths exist in parallel and are validated against one another. A reviewer auditing a specific value can locate it via whichever path matches their workflow.

**Script quick reference.**

| Script | Phase | Inputs | Outputs | Validation |
|---|---|---|---|---|
| `scripts/evaluate_severity_panel.py` | 2 | Raw model prediction h5ads in `data/predictions/raw/`; atlases in `data/replogle/`; severity references in `data/severity_refs/` | Per-seed severity-detail h5ads in `data/predictions/severity_details/` (the `.obs`-column contract that Phase 1 consumes) | `--validate-existing` compares freshly-derived detail h5ads against the existing reference set; verified 428 of 428 match within tolerance |
| `scripts/recompute_diagnostics.py` | 1 | Per-seed severity-detail h5ads + per-cell-type severity references | All 10 manuscript-canonical CSVs in `precomputed/eval/`, `precomputed/tables/`, and `precomputed/figure_inputs/` | `--validate` diffs every output against `precomputed/`; verified exact match (CI summary within Monte Carlo tolerance) |

**Provenance chain for every figure/table/in-text claim:**

```
raw atlas + holdout spec  →  raw prediction h5ad
                                    ↓ (Phase 2: scripts/evaluate_severity_panel.py)
                              severity-detail h5ad (.obs columns: predicted_mean_abs_delta,
                                                                  leverage_score,
                                                                  perturbation_target)
                                    ↓ (Phase 1: scripts/recompute_diagnostics.py)
                              precomputed/{eval, tables, figure_inputs}/*.csv
                                    ↓ (notebooks 01–04)
                              Figure 1, Table 1, Table 2, Appendix A
```

**Output → producing script:**

| Output file | Produced by | Notebook that consumes it |
|---|---|---|
| `precomputed/eval/diag_loo_sensitivity_n100.csv` | `scripts/recompute_diagnostics.py` (Pass 1) | 01 (cell 4), 03 (cell 4), 04 (cell 4) |
| `precomputed/eval/diag_loo_sensitivity_gears.csv` | `scripts/recompute_diagnostics.py` (Pass 1) | 03 (cell 4) |
| `precomputed/eval/diag_winsorise_n100.csv` | `scripts/recompute_diagnostics.py` (Pass 2) | — (analysis intermediate) |
| `precomputed/eval/diag_winsorise_n100_summary.csv` | `scripts/recompute_diagnostics.py` (Pass 2) | 01 (cell 4), 02 (cell 4) |
| `precomputed/eval/stage5_comparison_n100.csv` | `scripts/recompute_diagnostics.py` (Pass 3) | 02 (cell 4) |
| `precomputed/tables/table1_mechanism_summary.csv` | `scripts/recompute_diagnostics.py` | 02 (cell 6) |
| `precomputed/tables/table2_gears_matched_n.csv` | `scripts/recompute_diagnostics.py` | 03 (cell 8) |
| `precomputed/tables/appendix_a_threshold_sensitivity.csv` | `scripts/recompute_diagnostics.py` | 04 (cell 6) |
| `precomputed/figure_inputs/bootstrap_ci_summary.csv` | `scripts/recompute_diagnostics.py` | 01 (cell 10) |
| `precomputed/figure_inputs/figure1_panel_summary.csv` | `scripts/recompute_diagnostics.py` | 01 (cell 8) |

The per-seed severity-detail h5ads that feed `scripts/recompute_diagnostics.py` are produced by `scripts/evaluate_severity_panel.py` from raw CPA and GEARS prediction outputs. Raw predictions in turn are produced by the P01 (CPA) and P02 (GEARS) training notebooks, whose substantive training loops remain documented as the next-engineering-session deliverable.

---

## Figure 1: Cell-type-specific instability and tail sensitivity of the CPA severity metric (n = 100 seeds)

**Notebook:** `01_figure1.ipynb` ← `precomputed/eval/diag_loo_sensitivity_n100.csv`, `precomputed/eval/diag_winsorise_n100_summary.csv`

### Figure panels

| Panel | Description | Notebook | Cell |
|-------|-------------|----------|------|
| 1A | Per-seed severity Pearson across 4 conditions (n=100) | 01 | 8 |
| 1B | Winsorisation effect (raw vs winsorised range) | 01 | 8 |
| 1C | LOO mechanism signature (LOO_max vs τ₁) | 01 | 8 |

### In-text statistics (Results §3.1)

| Claim | Value | Notebook | Cell |
|-------|-------|----------|------|
| K562 random range −0.487 to +0.607 | span 1.094 | 01 | 4 |
| K562 random median +0.060 | +0.060 | 02 | 4 |
| K562 random sign flips | 42/100 | 02 | 4 |
| K562 random bootstrap 95% mean CI | [−0.004, +0.089] | 01 | 10 |
| RPE1 random range −0.135 to +0.620 | span 0.755 | 01 | 4 |
| RPE1 random median +0.290 | +0.290 | 02 | 4 |
| RPE1 random sign flips | 5/100 | 02 | 4 |
| RPE1 random bootstrap 95% mean CI | [+0.238, +0.300] | 01 | 10 |

### In-text statistics (Results §3.2)

| Claim | Value | Notebook | Cell |
|-------|-------|----------|------|
| K562 random LOO_max median | 0.163 | 02 | 4 |
| K562 random τ₁ median | 0.509 | 02 | 4 |
| K562 random high-LOO_max holdouts | 56/100 | 02 | 4 |
| K562 random MED12 mode count | 23/56 | 02 | 6 |
| RPE1 random LOO_max median | 0.088 | 02 | 4 |
| RPE1 random τ₁ median | 0.429 | 02 | 4 |
| RPE1 random high-LOO_max holdouts | 13/100 | 02 | 4 |
| RPE1 random SF3B2 mode count | 8/13 | 02 | 6 |
| RPE1 stratified SF3B2 mode count | 8/9 | 02 | 6 |
| K562 random POLR3A at seeds 1–7 (label-shift vignette) | mixed (POLR3A) | 03 | 8 |

### In-text statistics (Results §3.3)

| Claim | Value | Notebook | Cell |
|-------|-------|----------|------|
| K562 random range drops 22% (1.094 → 0.849) | −22% | 02 | 4 |
| K562 random MAD drops 21% (0.174 → 0.137) | −21% | 02 | 4 |
| K562 random 5–95 width drops 28% (0.763 → 0.547) | −28% | 02 | 4 |
| K562 random sign flips 42 → 22 under winsorisation | 22 | 02 | 4 |
| K562 stratified range drops 16% | −16% | 02 | 4 |
| K562 random median +0.060 → +0.148 under winsorisation | +0.148 | 02 | 4 |
| K562 stratified median +0.101 → +0.152 under winsorisation | +0.152 | 02 | 4 |
| RPE1 random range +4%, MAD +10%, width +9% | marginal widening | 02 | 4 |
| K562 median capped perturbations per holdout | 6 of 30 | 02 | 4 |
| RPE1 median capped perturbations per holdout | 1 of 30 | 02 | 4 |

### In-text statistics (Results §3.4)

| Claim | Value | Notebook | Cell |
|-------|-------|----------|------|
| K562 LOO_max median identical under both splits | 0.163 / 0.163 | 02 | 4 |
| K562 sign flips under random vs stratified | 42 vs 40 | 02 | 4 |
| K562 stratified MED17 mode count | 21/53 | 02 | 6 |
| K562 MAD random vs stratified | 0.174 vs 0.152 | 02 | 6 |
| RPE1 MAD random vs stratified | 0.100 vs 0.112 | 02 | 6 |

---

## Table 1: Mechanism summary across the four CPA n = 100 conditions

**Notebook:** `02_table1_mechanism_summary.ipynb` ← `precomputed/eval/stage5_comparison_n100.csv`, `precomputed/eval/diag_winsorise_n100_summary.csv`

Table 1 is built in **cell 6** (`02`) and saved to `precomputed/tables/table1_mechanism_summary.csv` in cell 8.

| Cell type | Split | Raw median | Sign flips | LOO_max median | Top driver (count) | Wins. range Δ | Cell |
|-----------|-------|-----------|-----------|----------------|-------------------|--------------|------|
| K562 | random | +0.060 | 42/100 | 0.163 | MED12 (23/56) | −22% | 02 | 6 |
| K562 | stratified | +0.101 | 40/100 | 0.163 | MED17 (21/53) | −16% | 02 | 6 |
| RPE1 | random | +0.290 | 5/100 | 0.088 | SF3B2 (8/13) | +4% | 02 | 6 |
| RPE1 | stratified | +0.263 | 8/100 | 0.092 | SF3B2 (8/9) | +5% | 02 | 6 |

---

## Table 2: Cross-architecture mechanism comparison at matched seed count

**Notebook:** `03_table2_gears_scope.ipynb` ← `precomputed/eval/diag_loo_sensitivity_n100.csv` (filtered to CPA seeds 1–7), `precomputed/eval/diag_loo_sensitivity_gears.csv`

Table 2 is **derived from raw data** by applying the mechanism label rule (cell 6) to the three slices (GEARS n=7, CPA matched n=7, CPA n=100) in cell 8, and saved to `precomputed/tables/table2_gears_matched_n.csv` in cell 10.

| Condition | GEARS (n=7) | CPA matched (n=7) | Verdict | CPA n=100 | Cell |
|-----------|-------------|-------------------|---------|-----------|------|
| K562 random | mixed (RPSA; 0.073) | mixed (POLR3A; 0.138) | match | single-driver (MED12; 0.163) | 03 | 8 |
| K562 stratified | mixed (SMC1A; 0.129) | single-driver (MED17; 0.217) | differ | single-driver (MED17; 0.163) | 03 | 8 |
| RPE1 random | mixed (BET1; 0.092) | mixed (—; 0.102) | match | mixed (SF3B2; 0.088) | 03 | 8 |
| RPE1 stratified | mixed (SF3B2; 0.076) | mixed (SF3B2; 0.067) | match | mixed (SF3B2; 0.092) | 03 | 8 |

### In-text statistics (Results §3.5)

| Claim | Value | Notebook | Cell |
|-------|-------|----------|------|
| Mechanism label agreement | 3 of 4 | 03 | 8 |
| RPE1 stratified cross-architecture driver | SF3B2 (both) | 03 | 8 |
| K562 stratified architecture-specific finding | CPA MED17 vs GEARS SMC1A | 03 | 8 |
| K562 random label shift between matched-n and n=100 | mixed → single-driver | 03 | 8 |

---

## Appendix A: LOO threshold sensitivity

**Notebook:** `04_appendix_threshold_sensitivity.ipynb` ← `precomputed/eval/diag_loo_sensitivity_n100.csv`

Appendix A counts are computed in **cell 6** (`04`) and saved to `precomputed/tables/appendix_a_threshold_sensitivity.csv` in cell 8.

| Threshold | K562 random | K562 stratified | RPE1 random | RPE1 stratified | Cell |
|-----------|-------------|-----------------|-------------|-----------------|------|
| LOO_max > 0.125 | 69/100 | 65/100 | 21/100 | 20/100 | 04 | 6 |
| LOO_max > 0.150 | 56/100 | 53/100 | 13/100 | 9/100 | 04 | 6 |
| LOO_max > 0.175 | 49/100 | 44/100 | 6/100 | 2/100 | 04 | 6 |

---

## Abstract claims

| Claim | Value | Notebook | Cell |
|-------|-------|----------|------|
| K562 span 1.09 | 1.0943 | 01 | 4 |
| RPE1 span 0.76 | 0.7552 | 01 | 4 |
| K562 sign flips | 42/100 | 02 | 4 |
| RPE1 sign flips | 5/100 | 02 | 4 |
| K562 winsorised median shift | +0.060 → +0.148 | 02 | 4 |
| GEARS matched-n agreement | 3 of 4 conditions | 03 | 8 |

---

## Manuscript figure provenance

The manuscript-bundled Figure 1 PNG at `/manuscripts/perturbation_ai/figs/figure1_methodology.png` is byte-identical to the repository-generated `figs/figure1_methodology.png` (SHA256 `8ff8697775b9f14bcc288edc6fc7c5d3830e81ffde132da6bd2c92988c3614dd`). The bundled PDF differs only by Matplotlib's `CreationDate` metadata; the rendered content is identical. The manuscript figure is updated by copying the repository output after running `./run_all.sh --figures`.

---

## Notebook implementation status

All eleven notebooks (D01–D03, P01–P04, 01–04) are implemented and in
`NOTEBOOKS_READY` in `run_all.sh`. The D-series and training notebooks
(P01, P02) gracefully exit with prerequisite-check messages when their
external inputs (atlases, Replogle supplementary Excel, holdout specs,
predictions) are not present, so partial setups never produce cryptic
stack traces.

One substantive gap remains: the CPA and GEARS training loops in
`src/perturb_eval/train_cpa.py` and `src/perturb_eval/train_gears.py` are
documented stubs. When P01 and P02 are invoked with all prerequisites
present, they print their planned per-seed flow rather than performing
actual training. The mirror from the original `a_perturb_AI` orchestrators
is queued for a future session and does not affect Levels 1 or 2 of the
reproducibility hierarchy described in `README.md`.

---

## Summary

| Category | Claims traced |
|----------|--------------|
| Figure 1 panels | 3 |
| Table 1 entries | 4 rows × 5 data cols = 20 |
| Table 2 entries | 4 rows × 4 data cols = 16 |
| Appendix A entries | 3 rows × 4 data cols = 12 |
| In-text statistics | ~30 |
| Abstract claims | 6 |
| **Total** | **~87** |

All values were verified against the executed notebook outputs and the
underlying CSVs in `precomputed/eval/`. The bootstrap 95% confidence intervals
reported in §3.1 are reproduced in `01_figure1.ipynb` cell 10 and saved to
`precomputed/figure_inputs/bootstrap_ci_summary.csv`. The reproduced bounds
match the manuscript-quoted values to within bootstrap Monte-Carlo precision
(~0.002 at 10,000 resamples with `seed=0`).
