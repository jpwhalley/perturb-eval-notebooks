"""
loo — Leave-one-out sensitivity scoring and condition-level mechanism labels.

The leave-one-out sensitivity score quantifies how much any single held-out
perturbation drives the severity correlation in a given holdout. For each
holdout with N=30 perturbations and observed severity correlation r, the score
recomputes r with each perturbation removed in turn, giving N leave-one-out
correlations r_{-i}. The per-perturbation impact is d_i = |r - r_{-i}| and the
maximum absolute change is LOO_max = max_i d_i. The single-driver concentration
τ_1 is d_(1) / (d_(1) + d_(2) + d_(3)).

Three condition-level mechanism labels are assigned across the seeds for a
condition (cell type × split × model):
    single-driver : median LOO_max > 0.15 AND the same perturbation appears as
                    the top driver in > 30% of high-LOO_max holdouts
    distributed   : median LOO_max < 0.10 AND no perturbation recurs above
                    that threshold
    mixed         : otherwise

Functions to add:
    loo_sensitivity_one(magnitudes, leverage) -> dict
        Per-holdout LOO_max, τ_1, and per-perturbation impacts.
    loo_sensitivity_panel(eval_csv: Path) -> pd.DataFrame
        Per-seed LOO scoring across all (condition, seed) rows.
    mechanism_label(loo_panel: pd.DataFrame, condition: str) -> str
        Apply the three-label rule to a condition × seeds slice.

Pending implementation. See P04_loo_and_winsorisation.ipynb.
"""

# TODO: implement loo_sensitivity_one, loo_sensitivity_panel, mechanism_label
