# Structure novel-variant test

> **Type:** model · **Approach:** `structure_test` · **1/1 tests pass** ✅

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E08 · Novel-variant holdout (S4) | ✅ | n_samples=312 | [report](e08_novel_variant/REPORT.md) |

## E08 · Novel-variant holdout (S4) — PASS ✅

Original structure achieves AUC ~0.56-0.57 (pre-registered negative baseline). Exploratory biophysical features show an unvalidated subset signal (AUC ~0.688). These features are NOT in the production model and have not been validated externally.

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `auc_plus_structure` | 0.5932 | 0.5700 | 0.0800 | ✅ |

| Metric | Value |
|---|---|
| `n_samples` | 312 |
| `n_distinct_sets` | 143 |
| `n_resistant` | 169 |
| `auc_region_count` | 0.5451 |
| `auc_plus_biophysics` | 0.6876 |
| `auc_full_novel` | 0.6694 |

_Caveats: DNA chains are a modeled complex. Grantham matrix captures radical vs conservative substitutions. Exploratory subset only (n=312); balanced class weighting deviates from production calibration._

