# Data preparation

> **Type:** data · **Approach:** `data_prep` · **1/1 tests pass** ✅

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E01 · Data preparation | ✅ | n_samples=11868 | [report](e01_prep/REPORT.md) |

## E01 · Data preparation — PASS ✅

Exact reproduction of the published cohort counts (11,868 samples / 866 resistant).

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `n_samples` | 11868 | 11868 | 0 | ✅ |
| `n_resistant` | 866 | 866 | 0 | ✅ |
| `cryptic_n` | 10704 | 10704 | 0 | ✅ |
| `cryptic_R` | 83 | 83 | 0 | ✅ |
| `nicd_n` | 1164 | 1164 | 0 | ✅ |
| `nicd_R` | 783 | 783 | 0 | ✅ |
| `unique_subjid` | 11611 | 11611 | 0 | ✅ |
| `mutation_rows` | 25554 | 25554 | 0 | ✅ |
| `minor_allele_rows` | 1115 | 1115 | 0 | ✅ |
| `null_rows` | 11 | 11 | 0 | ✅ |
| `null_rows_samples` | 4 | 4 | 0 | ✅ |

| Metric | Value |
|---|---|
| `median_FRSf_minor` | 0.3750 |

_Caveats: median_FRSf_minor (~0.375) is approximate, not an exact acceptance row. Repo DST files hold 22,347 R/S IDs; the QUALITY=HIGH + genome join yields 11,868._

