"""
train_gears — GEARS training wrapper.

Functions to add:
    train_one_seed(cell_type: str, split: str, seed: int, n_epochs: int = 30,
                   output_dir: Path = Path('data/predictions')) -> Path
        Train GEARS for a single (cell_type, split, seed) combination using the
        deterministic holdout specification from D03. Used for the matched-n
        cross-architecture scope check; only seeds 1 through 7 are evaluated.
    train_panel(cell_types: list[str], splits: list[str], seeds: list[int],
                **kwargs) -> list[Path]
        Convenience wrapper.

Pending implementation. See P02_train_gears.ipynb for the orchestration.

Note: this module requires the `gears` dependency group:
    uv sync --group gears
"""

# TODO: implement train_one_seed, train_panel
