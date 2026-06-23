"""
perturb_eval — Reusable analysis modules for the perturbation-evaluation audit.

Submodules:
    data         — Load Replogle atlases and harmonised severity references.
    holdouts     — Generate and read random / leverage-quintile-stratified holdout specifications.
    train_cpa    — CPA training wrapper, parameterised by holdout and seed.
    train_gears  — GEARS training wrapper, parameterised by holdout and seed.
    severity     — Compute predicted-magnitude and severity Pearson correlation.
    loo          — Leave-one-out sensitivity scoring and condition-level mechanism labels.
    winsorise    — Cell-type-specific p5/p95 leverage clipping and the winsorisation panel.
    plotting     — Figure-building helpers used by the figure notebooks.

The notebooks parameterise and call these modules; the modules carry the
analysis logic. This keeps notebooks readable and the implementations
testable in isolation.
"""

__version__ = "0.1.0"
