# Manuscript–Notebook Traceability Matrix

> **AI-assisted documentation.** This traceability matrix and the cross-references
> between manuscript claims and notebook cells were assembled with the assistance
> of Claude (Anthropic) via Cowork. All scientific analyses, interpretations, and
> editorial decisions were made by the author. The AI tools were used to
> systematically audit the mapping between code outputs and manuscript text.
> Readers can use this document to independently verify any claim in the
> manuscript against the corresponding notebook cell.

This document maps every figure panel, table entry, and in-text statistic in
the manuscript to the notebook and cell that produces it.

**Pipeline:** D-series (acquire) → P-series (analyse) → Numbered (visualise).
See `README.md` for full pipeline documentation.

**Status:** scaffolded; notebook cell references will be filled in as each
notebook is implemented.

---

## Figure 1: Cell-type-specific instability and tail sensitivity of the CPA severity metric (n = 100 seeds)

**Notebook:** `01_figure1.ipynb` ← `P03_evaluate_severity.ipynb`, `P04_loo_and_winsorisation.ipynb`

### Figure panels

| Panel | Description | Notebook | Cell |
|-------|-------------|----------|------|
| 1A | Per-seed severity Pearson across 4 conditions (n=100) | 01 | TODO |
| 1B | Winsorisation effect (raw vs winsorised range) | 01 | TODO |
| 1C | LOO mechanism signature (LOO_max vs τ₁) | 01 | TODO |

### In-text statistics (Results §3.1)

| Claim | Value | Notebook | Cell |
|-------|-------|----------|------|
| K562 random range −0.487 to +0.607 | span 1.094 | 01 | TODO |
| K562 random median +0.060 | +0.060 | 01 | TODO |
| K562 random sign flips | 42/100 | 01 | TODO |
| K562 random bootstrap 95% mean CI | [−0.002, +0.090] | 01 | TODO |
| RPE1 random range −0.135 to +0.620 | span 0.755 | 01 | TODO |
| RPE1 random median +0.290 | +0.290 | 01 | TODO |
| RPE1 random sign flips | 5/100 | 01 | TODO |
| RPE1 random bootstrap 95% mean CI | [+0.238, +0.301] | 01 | TODO |

### In-text statistics (Results §3.2)

| Claim | Value | Notebook | Cell |
|-------|-------|----------|------|
| K562 random LOO_max median | 0.163 | 02 | TODO |
| K562 random τ₁ median | 0.509 | 02 | TODO |
| K562 random high-LOO_max holdouts | 56/100 | 02 | TODO |
| K562 random MED12 mode count | 23/56 | 02 | TODO |
| RPE1 random LOO_max median | 0.088 | 02 | TODO |
| RPE1 random τ₁ median | 0.429 | 02 | TODO |
| RPE1 random high-LOO_max holdouts | 13/100 | 02 | TODO |
| RPE1 random SF3B2 mode count | 8/13 | 02 | TODO |
| RPE1 stratified SF3B2 mode count | 8/9 | 02 | TODO |
| K562 random POLR3A at seeds 1–7 (label-shift vignette) | mixed (POLR3A) | 02 | TODO |

### In-text statistics (Results §3.3)

| Claim | Value | Notebook | Cell |
|-------|-------|----------|------|
| K562 random range drops 22% (1.094 → 0.849) | −22% | 01 | TODO |
| K562 random MAD drops 21% (0.174 → 0.137) | −21% | 01 | TODO |
| K562 random 5–95 width drops 28% (0.763 → 0.547) | −28% | 01 | TODO |
| K562 random sign flips 42 → 22 under winsorisation | 22 | 01 | TODO |
| K562 stratified range drops 16% | −16% | 01 | TODO |
| K562 random median +0.060 → +0.148 under winsorisation | +0.148 | 01 | TODO |
| K562 stratified median +0.101 → +0.152 under winsorisation | +0.152 | 01 | TODO |
| RPE1 random range +4%, MAD +10%, width +9% | marginal widening | 01 | TODO |
| K562 median capped perturbations per holdout | 6 of 30 | 01 | TODO |
| RPE1 median capped perturbations per holdout | 1 of 30 | 01 | TODO |

### In-text statistics (Results §3.4)

| Claim | Value | Notebook | Cell |
|-------|-------|----------|------|
| K562 LOO_max median identical under both splits | 0.163 / 0.163 | 02 | TODO |
| K562 sign flips under random vs stratified | 42 vs 40 | 02 | TODO |
| K562 stratified MED17 mode count | 21/53 | 02 | TODO |
| K562 MAD random vs stratified | 0.174 vs 0.152 | 02 | TODO |
| RPE1 MAD random vs stratified | 0.100 vs 0.112 | 02 | TODO |

---

## Table 1: Mechanism summary across the four CPA n = 100 conditions

**Notebook:** `02_table1_mechanism_summary.ipynb` ← `P03_evaluate_severity.ipynb`, `P04_loo_and_winsorisation.ipynb`

| Cell type | Split | Raw median | Sign flips | LOO_max median | Top driver (count) | Wins. range Δ | Cell |
|-----------|-------|-----------|-----------|----------------|-------------------|--------------|------|
| K562 | random | +0.060 | 42/100 | 0.163 | MED12 (23/56) | −22% | TODO |
| K562 | stratified | +0.101 | 40/100 | 0.163 | MED17 (21/53) | −16% | TODO |
| RPE1 | random | +0.290 | 5/100 | 0.088 | SF3B2 (8/13) | +4% | TODO |
| RPE1 | stratified | +0.263 | 8/100 | 0.092 | SF3B2 (8/9) | +5% | TODO |

---

## Table 2: Cross-architecture mechanism comparison at matched seed count

**Notebook:** `03_table2_gears_scope.ipynb` ← `P02_train_gears.ipynb`, `P03_evaluate_severity.ipynb`, `P04_loo_and_winsorisation.ipynb`

| Condition | GEARS (n=7) | CPA matched (n=7) | Verdict | CPA n=100 | Cell |
|-----------|-------------|-------------------|---------|-----------|------|
| K562 random | mixed (RPSA; 0.073) | mixed (POLR3A; 0.138) | match | single-driver (MED12; 0.163) | TODO |
| K562 stratified | mixed (SMC1A; 0.129) | single-driver (MED17; 0.217) | differ | single-driver (MED17; 0.163) | TODO |
| RPE1 random | mixed (BET1; 0.092) | mixed (—; 0.102) | match | mixed (SF3B2; 0.088) | TODO |
| RPE1 stratified | mixed (SF3B2; 0.076) | mixed (SF3B2; 0.067) | match | mixed (SF3B2; 0.092) | TODO |

### In-text statistics (Results §3.5)

| Claim | Value | Notebook | Cell |
|-------|-------|----------|------|
| Mechanism label agreement | 3 of 4 | 03 | TODO |
| RPE1 stratified cross-architecture driver | SF3B2 (both) | 03 | TODO |
| K562 stratified architecture-specific finding | CPA MED17 vs GEARS SMC1A | 03 | TODO |
| K562 random label shift between matched-n and n=100 | mixed → single-driver | 03 | TODO |

---

## Appendix A: LOO threshold sensitivity

**Notebook:** `04_appendix_threshold_sensitivity.ipynb` ← `P04_loo_and_winsorisation.ipynb`

| Threshold | K562 random | K562 stratified | RPE1 random | RPE1 stratified | Cell |
|-----------|-------------|-----------------|-------------|-----------------|------|
| LOO_max > 0.125 | 69/100 | 65/100 | 21/100 | 20/100 | TODO |
| LOO_max > 0.150 | 56/100 | 53/100 | 13/100 | 9/100 | TODO |
| LOO_max > 0.175 | 49/100 | 44/100 | 6/100 | 2/100 | TODO |

---

## Abstract claims

| Claim | Value | Notebook | Cell |
|-------|-------|----------|------|
| K562 span 1.09 | 1.0943 | 01 | TODO |
| RPE1 span 0.76 | 0.7552 | 01 | TODO |
| K562 sign flips | 42/100 | 01 | TODO |
| RPE1 sign flips | 5/100 | 01 | TODO |
| K562 winsorised median shift | +0.060 → +0.148 | 01 | TODO |
| GEARS matched-n agreement | 3 of 4 conditions | 03 | TODO |

---

## Summary

| Category | Claims traced |
|----------|--------------|
| Figure 1 panels | 3 |
| Table 1 entries | 4 rows × 5 data cols = 20 |
| Table 2 entries | 4 rows × 4 data cols = 16 |
| Appendix A entries | 3 rows × 4 data cols = 12 |
| In-text statistics | ~30 |
| **Total** | **~80** |

All values were verified against the evaluation CSVs in `precomputed/eval/` at
the manuscript's final pre-submission state. Cell references will be added
as each notebook is implemented.
