# Tier A — catalogue

> **Type:** baseline · **Approach:** `catalogue` · **1/1 tests pass** ✅

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E03 · Catalogue CV (S1, patient-grouped) | ✅ | sens_firm=0.7924 | [report](e03_s1_catalogue_cv/REPORT.md) |

## E03 · Catalogue CV (S1, patient-grouped) — PASS ✅

Rebuilt catalogue under patient-grouped CV (stricter than the paper's shuffled CV).

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `firm` | 0.9443 | 0.9440 | 0.0100 | ✅ |
| `sens_firm` | 0.7924 | 0.7920 | 0.0150 | ✅ |
| `spec_firm` | 0.9835 | 0.9840 | 0.0040 | ✅ |
| `sens_U_as_S` | 0.6039 | 0.6040 | 0.0150 | ✅ |
| `spec_U_as_S` | 0.9842 | 0.9840 | 0.0040 | ✅ |
| `sens_U_as_R` | 0.8418 | 0.8420 | 0.0150 | ✅ |
| `spec_U_as_R` | 0.9428 | 0.9430 | 0.0100 | ✅ |
| `n_U` | 661 | 661 | 15 | ✅ |
| `n_U_resistant` | 206 | 206 | 15 | ✅ |

_Caveats: In-sample (built on all data) would be optimistic; this is out-of-fold._

