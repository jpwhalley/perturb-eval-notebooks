"""
plotting — Figure-building helpers for the figure notebooks.

Functions to add:
    figure1_panel_a(ax, eval_df) — Per-seed severity Pearson across 4 conditions.
    figure1_panel_b(ax, wins_summary) — Raw vs winsorised range, paired bars.
    figure1_panel_c(ax, loo_df) — LOO_max vs τ_1 mechanism scatter.
    build_figure1(eval_df, wins_summary, loo_df) — Assemble all three panels.
    table1(eval_df, loo_df, wins_summary) — Mechanism summary CSV.
    table2(gears_eval, cpa_eval_matched, cpa_eval_n100) — GEARS matched-n verdict.
    appendix_a_threshold_counts(loo_df) — Threshold sensitivity table.

All figure helpers respect the `perturb_style` palette and sizing. See
`perturb_style.py` for the colour and marker scheme.

Pending implementation. See 01_figure1.ipynb through 04_appendix_threshold_sensitivity.ipynb.
"""

# TODO: implement figure1_*, table1, table2, appendix_a_threshold_counts
