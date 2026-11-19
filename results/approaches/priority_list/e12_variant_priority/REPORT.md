# E12 · Variant priority list

> **Approach:** `priority_list` · **Test:** `e12_variant_priority` · **PASS ✅**

332 distinct sets cause 661 referrals; top-10 sets cover ~34%. Phenotyping a handful of variants would cut referrals sharply (descriptive).

## Metrics

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

## Method
Among catalogue-U samples, group by sorted core-variant set; rank by carrier count; top-K coverage. Per-variant Wilson 95% CI + by-source counts.

## Novelty claim
Referral concentration (basis for the VOI queue novelty).

## Caveats / negative results
NICD enriches resistance; evidence cards split by source. Descriptive, not causal.

## Reproduce
```bash
python eval/e12_priority.py
```

## Provenance
- **Date (UTC):** 2026-11-18T19:05:50.721667+00:00
- **Git commit:** uncommitted
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `top_sets.csv`
- `evidence_cards.csv`
