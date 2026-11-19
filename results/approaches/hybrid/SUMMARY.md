# Hybrid (exploratory)

> **Type:** analysis · **Approach:** `hybrid` · **1/1 tests pass** ✅

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E09 · Hybrid on catalogue-U | ✅ | auc_ml_within_U=0.7442 | [report](e09_hybrid/REPORT.md) |

## E09 · Hybrid on catalogue-U — PASS ✅

EXPLORATORY: model ranks U samples (AUC ~0.74) but resolving them is unsafe (resolved-S group ~9.8% resistant; resolved-R ~40% susceptible). NEVER auto-call.

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `auc_ml_within_U` | 0.7442 | 0.7440 | 0.0300 | ✅ |
| `resolved_S` | 173 | 173 | 15 | ✅ |
| `resolved_S_resistant` | 17 | 17 | 15 | ✅ |
| `resolved_R` | 126 | 126 | 15 | ✅ |
| `resolved_R_susceptible` | 51 | 51 | 15 | ✅ |
| `remain_U` | 362 | 362 | 15 | ✅ |
| `remain_U_resistant` | 114 | 114 | 15 | ✅ |

| Metric | Value |
|---|---|
| `n_U` | 661 |
| `n_U_resistant` | 206 |

_Caveats: Used only to RANK the unresolved queue, not to change calls._

