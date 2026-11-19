# Variant priority list

> **Type:** analysis · **Approach:** `priority_list` · **1/1 tests pass** ✅

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E12 · Variant priority list | ✅ | n_distinct_sets=332 | [report](e12_variant_priority/REPORT.md) |

## E12 · Variant priority list — PASS ✅

332 distinct sets cause 661 referrals; top-10 sets cover ~34%. Phenotyping a handful of variants would cut referrals sharply (descriptive).

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `n_distinct_sets` | 332 | 332 | 15 | ✅ |
| `top5_coverage` | 0.2769 | 0.2800 | 0.0300 | ✅ |
| `top10_coverage` | 0.3434 | 0.3400 | 0.0300 | ✅ |
| `top20_coverage` | 0.3979 | 0.4000 | 0.0300 | ✅ |
| `top50_coverage` | 0.4962 | 0.5000 | 0.0300 | ✅ |

| Metric | Value |
|---|---|
| `n_U` | 661 |
| `largest_set_n` | 102 |
| `largest_set_pctR` | 0.0882 |

- **largest_set:** `Rv0678@c-11a`

_Caveats: NICD enriches resistance; evidence cards split by source. Descriptive, not causal._

