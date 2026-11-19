# E15 · Label Reveal Sim

> **Approach:** `voi_queue` · **Test:** `e15_label_reveal_sim` · **DESCRIPTIVE •**

Info-gain vs count-only %U-resolved: K=10 0.381 vs 0.381 (delta +0.000); K=50 0.537 vs 0.526 (delta +0.011).

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
score(set)=n_U_carriers x expected beta-binomial entropy reduction after one label (prior from non-U firm evidence only). Retrospective reveal of top-K sets' labels; variant unlock via binomial-test proxy (background 0.1, alpha 0.05, min 3). Baselines: count-only, random (seed 1). K=1,2,5,10,20,50.

## Novelty claim
Value-of-information lab queue (novelty P-E15): prioritizes phenotyping by expected information per referral cleared.

## Caveats / negative results
Retrospective, not prospective; assumes catalogued variants transfer across lineages. Variant unlock is a binomial-test proxy for catomatic strict_unlock. If info-gain barely beats count-only, that is a valid negative result. Priors exclude U labels (no leakage).

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
