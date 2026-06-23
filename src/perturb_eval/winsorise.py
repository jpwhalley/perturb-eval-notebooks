"""
winsorise — Cell-type-specific p5/p95 leverage clipping and the winsorisation panel.

Winsorising the leverage target clips the per-perturbation leverage at the
cell-type-specific 5th and 95th percentile of the full Replogle severity
reference before computing the Pearson correlation. This isolates whether the
per-seed instability is anchored on leverage-tail perturbations.

Functions to add:
    panel_thresholds(cell_type: str) -> tuple[float, float]
        Return (p5, p95) for the cell-type-specific severity reference.
    winsorise_severity_one(magnitudes, leverage, p5, p95) -> float
        Recompute the severity Pearson with leverage clipped at p5/p95.
    winsorise_panel(predictions_dir: Path) -> pd.DataFrame
        Per-seed raw vs winsorised Pearson across all 400 CPA n=100 conditions
        × seeds. Produces the input for Figure 1B and Table 1.

Pending implementation. See P04_loo_and_winsorisation.ipynb.
"""

# TODO: implement panel_thresholds, winsorise_severity_one, winsorise_panel
