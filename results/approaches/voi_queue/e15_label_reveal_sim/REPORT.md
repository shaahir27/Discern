# E15 · Retrospective comparison of variant sets

> **Approach:** `voi_queue` · **Test:** `e15_label_reveal_sim` · **DESCRIPTIVE •**

The strategies compare equal numbers of variant sets. Those sets contain different numbers of samples.

## Metrics

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

## Method
Retrospective comparison of variant sets. K counts selected sets. Each selected set reveals every member's saved lab label. The simulated total resolved includes those directly revealed samples and unrevealed samples counted by the proxy.

## Novelty claim
Descriptive archive simulation. The resolution proxy does not amend the catalogue.

## Caveats / negative results
The comparison uses unequal numbers of sample labels at the same K. Unrevealed samples counted as resolved include cases the proxy can already count using firm-case evidence alone. These totals cannot establish gains at an equal lab testing budget or production catalogue performance.

## Reproduce
```bash
python eval/e15_voi.py
```

## Provenance
- **Date (UTC):** 2026-11-19T04:28:47.172915+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `resolution_curves.csv`
- `info_gain_ranking.csv`
