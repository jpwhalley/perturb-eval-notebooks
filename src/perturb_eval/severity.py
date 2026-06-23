"""
severity — Compute predicted-magnitude and severity Pearson correlation.

For each held-out perturbation, the model produces per-cell predicted
post-perturbation expression. The predicted magnitude is the mean absolute
delta from the unperturbed control baseline, averaged across the evaluated
gene set. The severity Pearson is the Pearson correlation between predicted
magnitudes and the per-perturbation leverage score across the 30 held-out
perturbations in each holdout.

Functions to add:
    predicted_magnitude(adata_pred, adata_ctrl, gene_set) -> dict[str, float]
        Per-perturbation mean absolute delta from control.
    severity_correlation(magnitudes, leverage, method='pearson') -> float
        Pearson (or Spearman) correlation between predicted magnitude and
        per-perturbation leverage across the 30 held-out perturbations.
    evaluate_seed(prediction_path: Path, severity_ref: pd.DataFrame) -> dict
        Combine the above into a per-seed evaluation record.

Pending implementation. See P03_evaluate_severity.ipynb.
"""

# TODO: implement predicted_magnitude, severity_correlation, evaluate_seed
