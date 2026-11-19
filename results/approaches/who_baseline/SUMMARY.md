# WHO v2 baseline

> **Type:** baseline · **Approach:** `who_baseline` · **2/2 tests pass** ✅

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E02 · WHO v2 as-is (FRS ≥ 0.1) | ✅ | sens_firm=0.8011 | [report](e02_who_as_is_frs0.1/REPORT.md) |
| E02 · WHO v2 as-is (FRS ≥ 0.75) | ✅ | sens_firm=0.5606 | [report](e02_who_as_is_frs0.75/REPORT.md) |

## E02 · WHO v2 as-is (FRS ≥ 0.1) — PASS ✅

WHO v2 as-is at FRS>=0.1, pooled. Anchor sens_U_as_S must match published.

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `firm` | 0.9607 | 0.9610 | 0.0050 | ✅ |
| `sens_firm` | 0.8011 | 0.8010 | 0.0050 | ✅ |
| `spec_firm` | 0.9812 | 0.9810 | 0.0050 | ✅ |
| `sens_U_as_S` | 0.6790 | 0.6790 | 0.0050 | ✅ |

| Metric | Value |
|---|---|
| `spec_U_as_S` | 0.9817 |
| `sens_U_as_R` | 0.8314 |
| `spec_U_as_R` | 0.9514 |
| `n_U` | 466 |
| `n_U_resistant` | 132 |

- **variant_call_counts:** `{'U': 467, 'S': 346, 'R': 229}`

_Caveats: U counted as S for the anchor sensitivity. NICD enriches resistance; read per-source._

## E02 · WHO v2 as-is (FRS ≥ 0.75) — PASS ✅

WHO v2 as-is at FRS>=0.75, pooled. Anchor sens_U_as_S must match published.

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `firm` | 0.9680 | 0.9680 | 0.0050 | ✅ |
| `sens_firm` | 0.5606 | 0.5610 | 0.0050 | ✅ |
| `spec_firm` | 0.9868 | 0.9870 | 0.0050 | ✅ |
| `sens_U_as_S` | 0.4965 | 0.4970 | 0.0050 | ✅ |

| Metric | Value |
|---|---|
| `spec_U_as_S` | 0.9872 |
| `sens_U_as_R` | 0.6109 |
| `spec_U_as_R` | 0.9616 |
| `n_U` | 380 |
| `n_U_resistant` | 99 |

- **variant_call_counts:** `{'U': 340, 'S': 311, 'R': 159}`

_Caveats: U counted as S for the anchor sensitivity. NICD enriches resistance; read per-source._

