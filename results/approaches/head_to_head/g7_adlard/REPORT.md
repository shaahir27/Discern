# G7 · Adlard

> **Approach:** `head_to_head` · **Test:** `g7_adlard` · **DESCRIPTIVE •**

Shuffled sens_firm 0.7955 vs patient-grouped 0.7924 (delta +0.0031). Patient grouping is the stricter, leakage-controlled protocol. Ghosh comparison infeasible (data not local).

## Metrics

| Metric | Value |
|---|---|
| `shuffled_firm` | 0.9449 |
| `shuffled_sens_firm` | 0.7955 |
| `shuffled_spec_firm` | 0.9838 |
| `shuffled_n_U` | 654 |
| `grouped_firm` | 0.9443 |
| `grouped_sens_firm` | 0.7924 |
| `grouped_spec_firm` | 0.9835 |
| `grouped_n_U` | 661 |
| `sens_firm_shuffled_minus_grouped` | 0.0031 |

- **ghosh_reconstructable:** `False`

## Method
Adlard-style catalogue rebuilt under shuffled 5-fold (their protocol) vs patient-grouped 5-fold (ours, E03), identical 11,868 samples. Ghosh 632 isolates NOT reconstructable from local data (Tables S9/S10 absent) -> not compared.

## Novelty claim
Same-data head-to-head vs Adlard under both protocols (P-G7).

## Caveats / negative results
Shuffled CV can leak related germs across folds (optimistic). Grouped controls patient duplicates but not clones/outbreaks fully. Ghosh omitted honestly; NEVER cross-cohort AUC.

## Reproduce
```bash
python eval/e18_head_to_head.py
```

## Provenance
- **Date (UTC):** 2026-11-19T04:37:36.940617+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `shuffled_vs_grouped.csv`
