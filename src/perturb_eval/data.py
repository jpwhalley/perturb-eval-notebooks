"""
data — Load Replogle atlases and severity references.

Functions to add:
    load_replogle_atlas(cell_type: str) -> AnnData
        Load the K562 or RPE1 harmonised Perturb-seq atlas from data/replogle/.
    load_severity_reference(cell_type: str) -> pd.DataFrame
        Load the per-perturbation severity reference (leverage, knockdown,
        DEG count) from data/severity_refs/.

Pending implementation. See D01_replogle_atlases.ipynb and
D02_severity_references.ipynb for the data acquisition flow that produces
the on-disk artefacts these functions read.
"""

# TODO: implement load_replogle_atlas, load_severity_reference
