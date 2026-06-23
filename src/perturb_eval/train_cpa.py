"""
train_cpa — CPA training wrapper.

Functions to add:
    train_one_seed(cell_type: str, split: str, seed: int, n_epochs: int = 30,
                   output_dir: Path = Path('data/predictions')) -> Path
        Train CPA for a single (cell_type, split, seed) combination using the
        deterministic holdout specification from D03. Returns the path to the
        per-seed prediction AnnData. Approximately 2 minutes per call on CPU.
    train_panel(cell_types: list[str], splits: list[str], seeds: list[int],
                **kwargs) -> list[Path]
        Convenience wrapper that calls train_one_seed for each combination
        and returns the list of output paths.

Pending implementation. See P01_train_cpa.ipynb for the orchestration.

Note: this module requires the `cpa` dependency group:
    uv sync --group cpa
"""

# TODO: implement train_one_seed, train_panel
