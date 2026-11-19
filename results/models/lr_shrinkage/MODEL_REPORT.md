# Tier B — shrinkage estimator

> **Type:** model · **Approach:** `lr_shrinkage` · **0/1 tests pass** ❌

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E14 · Shrinkage | — | n_variants=1042 | [report](e14_shrinkage/REPORT.md) |

## E14 · Shrinkage — DESCRIPTIVE •

1042 variants shrunk; mean |shrunk-raw| 0.1432; 48 variants flip a 0.5 grade. Rare variants pulled toward class mean; standard estimator choice, not an invention.

| Metric | Value |
|---|---|
| `n_variants` | 1042 |
| `mean_abs_shrinkage` | 0.1432 |
| `max_abs_shrinkage` | 0.9294 |
| `n_grade_flips_0.5` | 48 |
| `n_mech_classes` | 9 |

_Caveats: Descriptive estimator output for evidence cards; ablate, don't headline. No accuracy claim. Class priors estimated from the same data (empirical Bayes)._

