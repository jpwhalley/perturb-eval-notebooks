# G2 — CPA vs GEARS LOO mechanism comparison (matched-n headline)

## HEADLINE — matched-n verdict (GEARS vs CPA on the SAME seeds; identical holdouts per seed)

| condition | GEARS (n=7) | CPA (matched seeds) | verdict | CPA (n=100, context) |
|---|---|---|---|---|
| K562 random | mixed [RPSA, LOO 0.073] | mixed [POLR3A, n=7, LOO 0.138] | MATCH | single-driver [MED12, LOO 0.163] |
| K562 stratified | mixed [SMC1A, LOO 0.129] | single-driver [MED17, n=7, LOO 0.217] | differ | single-driver [MED17, LOO 0.163] |
| RPE1 random | mixed [BET1, LOO 0.092] | mixed [(none), n=7, LOO 0.102] | MATCH | mixed [SF3B2, LOO 0.088] |
| RPE1 stratified | mixed [SF3B2, LOO 0.076] | mixed [SF3B2, n=7, LOO 0.067] | MATCH | mixed [SF3B2, LOO 0.092] |

**Matched-n agreement: 3/4 conditions MATCH.** CPA n=100 is context only (its labels can differ from matched-n because the mechanism_label rule is sample-size sensitive).

mechanism_label rule (both models): single-driver if LOO_max_median>0.15 AND modal high-LOO driver in >30% of high-LOO seeds; distributed if LOO_max_median<0.10 AND no driver >30%; else mixed.

## Context — full-n table (CPA n=100 vs GEARS n=7)

| condition | model | n | LOO_max_median | top_driver_mode | sign_flip_rate | mechanism_label |
|---|---|---|---|---|---|---|
| K562 random | CPA | 100 | 0.163 | MED12 | 42/100 | single-driver |
| K562 random | GEARS | 7 | 0.073 | RPSA | 2/7 | mixed |
| K562 stratified | CPA | 100 | 0.163 | MED17 | 40/100 | single-driver |
| K562 stratified | GEARS | 7 | 0.129 | SMC1A | 2/7 | mixed |
| RPE1 random | CPA | 100 | 0.088 | SF3B2 | 5/100 | mixed |
| RPE1 random | GEARS | 7 | 0.092 | BET1 | 3/7 | mixed |
| RPE1 stratified | CPA | 100 | 0.092 | SF3B2 | 8/100 | mixed |
| RPE1 stratified | GEARS | 7 | 0.076 | SF3B2 | 2/7 | mixed |

## Context — pre-registered reading rule (GEARS labels vs expected pattern)

GEARS K562 labels: ['mixed', 'mixed'] (drivers ['RPSA', 'SMC1A'])
GEARS RPE1 labels: ['mixed', 'mixed'] (drivers ['BET1', 'SF3B2'])

**OUTCOME: PARTIAL ALIGN**
