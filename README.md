# Perturbation Evaluation Audit — Notebooks

Reproducibility notebooks for:

> **Single-shot perturbation benchmarks can mislead: diagnostics for split- and metric-driven instability in single-cell model evaluation**
> Whalley, J.P. (2026). Submitted to MLCB 2026.

---

## Default path: reproduce all manuscript figures and tables in two minutes

```bash
git clone <repo-url>
cd repository/
uv sync
./run_all.sh --figures
```

This reproduces **Figure 1, Table 1, Table 2, and Appendix A** from the shipped CSVs in `precomputed/eval/`. No external data, no model training, no GPU. Approximate runtime: 2 minutes on a laptop.

Everything beyond this point is optional.

---

## Overview

The repository contains 11 notebooks organised in three layers, plus a precomputed CSV layer that lets a reader reproduce every manuscript-facing output without first running the heavy training pipeline:

```
D01–D03  →  P01–P04  →  01–04
(acquire)   (analyse)   (visualise)
```

| Layer | Notebooks | Purpose |
|---|---|---|
| Data acquisition | D01–D03 | Download Replogle atlases (D01), build per-gene severity references from the Replogle supplementary tables (D02), generate the deterministic holdout specifications (D03). |
| Analysis pipeline | P01–P04 | Train CPA (P01) and the matched-n GEARS audit (P02); compute per-seed severity Pearson (P03); compute leave-one-out sensitivity and the winsorisation panel (P04). |
| Figure / table notebooks | 01–04 | Consume CSVs and produce Figure 1, Table 1, Table 2, and Appendix A. |

The **public reader path** uses only the figure notebooks against the shipped `precomputed/eval/` CSVs. The D-series and P-series are needed only by readers who want to regenerate the analysis layer from scratch.

---

## Run modes

```bash
./run_all.sh                  # print help (does not start any compute)
./run_all.sh --figures        # public default; figures + tables from precomputed CSVs
./run_all.sh --demo           # training-path smoke test (see below)
./run_all.sh --analysis       # recompute eval CSVs from existing predictions
./run_all.sh --train-full     # full n=100 CPA + matched-n GEARS (~13 hours)
```

**Mode semantics:**

- **`--figures`** — public default. Runs notebooks 01–04 against the precomputed CSVs in `precomputed/eval/`. No external data, no model checkpoints, no GPU. Approximate runtime: 2 minutes.
- **`--demo`** — training-path smoke test. Runs the full chain (D-series → P-series → figures) with `QUICK_DEMO=1` set. Each notebook *checks its prerequisites and exits gracefully with installation guidance if anything is missing*. To actually train CPA on 2 seeds per condition, you need to have completed (a) `uv sync --group cpa --group gears --group data`, (b) the manual Replogle Excel download documented in D02, and (c) a successful D01 atlas download. The smoke test does not perform any of these for you; it surfaces what is missing.
- **`--analysis`** — recompute eval CSVs from existing per-seed predictions in `data/predictions/severity_details/`. P03 and P04 produce the CSVs that the figure notebooks consume. Useful after a successful `--demo` or `--train-full` run.
- **`--train-full`** — full reproduction: 100 CPA seeds per condition plus the matched-n GEARS audit. Approximately 13 hours on Apple Silicon CPU. Requires all dependency groups (`cpa`, `gears`, `data`).

---

## Reproducibility levels

The repository supports three increasingly heavy levels of reproduction:

| Level | Command | Time | What you reproduce | What you need |
|---|---|---|---|---|
| **1** | `./run_all.sh --figures` | ~2 min | Figure 1, Table 1, Table 2, Appendix A | Cloned repo + `uv sync` |
| **2** | `./run_all.sh --analysis` | minutes (per condition) | Level 1 plus the per-seed evaluation CSVs in `data/eval/` | Level 1 plus a prior `--demo` or `--train-full` run (which writes `data/predictions/`) |
| **3** | `./run_all.sh --train-full` | ~13 hours | Level 2 plus the per-seed CPA and GEARS prediction artefacts in `data/predictions/` | Level 2 prerequisites plus model dependency groups (`cpa`, `gears`, `data`) and a manual Replogle Excel download (see D02) |

Level 1 is the default and the path the manuscript primarily expects readers to use. Level 3 is the path the manuscript's reported results were originally generated through.

---

## Prerequisites

- Python 3.11
- [uv](https://docs.astral.sh/uv/) for dependency management

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Level 1 requires only the base dependency group (installed automatically by `uv sync`). Levels 2 and 3 additionally require:

```bash
uv sync --group cpa --group gears --group data
```

For Level 3 (full retraining) also recommended:

- Apple Silicon (MPS) or CUDA GPU
- ~20 GB disk space for atlases, predictions, and logs

---

## Repository structure

```
repository/
├── pyproject.toml
├── uv.lock
├── README.md
├── MANUSCRIPT_TRACEABILITY.md
├── perturb_style.py
├── run_all.sh
├── LICENSE
├── src/perturb_eval/         # Reusable analysis modules
├── precomputed/              # Tracked lightweight CSVs (Level 1 reproducibility layer)
│   ├── eval/                 # Per-seed evaluation tables consumed by figure notebooks
│   ├── tables/               # Table 1, Table 2, Appendix A as CSVs
│   └── figure_inputs/        # Aggregated per-condition values plotted in Figure 1
├── data/                     # gitignored; regenerated by D-series and P-series
├── D01_replogle_atlases.ipynb
├── D02_severity_references.ipynb
├── D03_holdout_specifications.ipynb
├── P01_train_cpa.ipynb
├── P02_train_gears.ipynb
├── P03_evaluate_severity.ipynb
├── P04_loo_and_winsorisation.ipynb
├── 01_figure1.ipynb
├── 02_table1_mechanism_summary.ipynb
├── 03_table2_gears_scope.ipynb
└── 04_appendix_threshold_sensitivity.ipynb
```

The `precomputed/` directory is the Level 1 reproducibility layer. Its CSVs are intentionally tracked so the figure notebooks can run without first reproducing the analysis pipeline. See `precomputed/README.md` for the curation rationale.

---

## Status

All 11 notebooks are implemented and gated by `NOTEBOOKS_READY` in `run_all.sh`. The public reader path (`--figures`) is verified end-to-end and reproduces every manuscript-facing output. The D-series and P-series notebooks handle missing prerequisites gracefully with installation guidance, so partial setups do not produce cryptic stack traces.

**Substantive gap (Level 3 only):** the CPA and GEARS training loops in `src/perturb_eval/train_cpa.py` and `src/perturb_eval/train_gears.py` are currently documented stubs. P01 and P02 will print their planned per-seed flow and exit cleanly when invoked in `--demo` or `--train-full` mode with all prerequisites present, but they do not yet actually train. The mirror from the original `a_perturb_AI` orchestrators is the last remaining piece for end-to-end Level 3 reproduction. Level 1 and the structure of Level 2 are unaffected.

---

## Manuscript traceability

`MANUSCRIPT_TRACEABILITY.md` maps every figure panel, table entry, in-text statistic, and abstract claim in the manuscript to the specific notebook and cell that produces it. Approximately 87 claims are traced. The document was assembled with the assistance of Claude (Anthropic) via Cowork.

## License

MIT. See `LICENSE`.
