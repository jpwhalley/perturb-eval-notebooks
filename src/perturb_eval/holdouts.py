"""
holdouts — Generate and read holdout specifications.

Two split strategies are supported:
    random      — uniform random draw of 30 perturbations from the 150-pool
    stratified  — leverage-quintile stratified draw (6 perturbations per quintile)

Functions to add:
    build_holdout_spec(cell_type: str, split: str, n_seeds: int) -> pd.DataFrame
        Generate the holdout specification table with one row per (seed,
        perturbation) and a deterministic random_state per seed.
    load_holdout_spec(cell_type: str, split: str) -> pd.DataFrame
        Read the pre-generated holdout specification from data/holdout_specs/.
    held_out_perturbations(spec: pd.DataFrame, seed: int) -> list[str]
        Look up the 30 held-out perturbations for a given seed.

Pending implementation. See D03_holdout_specifications.ipynb for usage.
"""

# TODO: implement build_holdout_spec, load_holdout_spec, held_out_perturbations
