# Value-of-information queue

> **Type:** analysis · **Approach:** `voi_queue` · **0/1 tests pass** ❌

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E15 · Retrospective comparison of variant sets | — | n_distinct_sets=332 | [report](e15_label_reveal_sim/REPORT.md) |

## E15 · Retrospective comparison of variant sets — DESCRIPTIVE •

The strategies compare equal numbers of variant sets. Those sets contain different numbers of samples.

| Metric | Value |
|---|---|
| `info_gain_pct_K1` | 0.2405 |
| `info_gain_pct_K2` | 0.2405 |
| `info_gain_pct_K5` | 0.3147 |
| `info_gain_pct_K10` | 0.3812 |
| `info_gain_pct_K20` | 0.4372 |
| `info_gain_pct_K50` | 0.5371 |
| `count_pct_K1` | 0.1921 |
| `count_pct_K2` | 0.2405 |
| `count_pct_K5` | 0.3147 |
| `count_pct_K10` | 0.3812 |
| `count_pct_K20` | 0.4372 |
| `count_pct_K50` | 0.5265 |
| `random_pct_K1` | 0.1936 |
| `random_pct_K2` | 0.2148 |
| `random_pct_K5` | 0.2209 |
| `random_pct_K10` | 0.2315 |
| `random_pct_K20` | 0.2466 |
| `random_pct_K50` | 0.2995 |
| `total_U` | 661 |
| `n_distinct_sets` | 332 |
| `info_gain_minus_count_K10` | 0.0000 |
| `info_gain_minus_count_K50` | 0.0106 |

_Caveats: The comparison uses unequal numbers of sample labels at the same K. Unrevealed samples counted as resolved include cases the proxy can already count using firm-case evidence alone. These totals cannot establish gains at an equal lab testing budget or production catalogue performance._

