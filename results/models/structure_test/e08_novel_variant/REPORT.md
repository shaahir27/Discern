# E08 · Novel-variant holdout (S4)

> **Approach:** `structure_test` · **Test:** `e08_novel_variant` · **PASS ✅**

Original structure achieves AUC ~0.56-0.57 (pre-registered negative baseline). Exploratory biophysical features show an unvalidated subset signal (AUC ~0.688). These features are NOT in the production model and have not been validated externally.

## Metrics

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `auc_plus_structure` | 0.5932 | 0.5700 | 0.0800 | ✅ |

| Metric | Value |
|---|---|
| `n_samples` | 312 |
| `n_distinct_sets` | 143 |
| `n_resistant` | 169 |
| `auc_region_count` | 0.5451 |
| `auc_plus_biophysics` | 0.6876 |
| `auc_full_novel` | 0.6694 |

## Method
Subset: Rv0678 all-missense, no atpE/pepQ (n=312, 143 clusters). GroupKFold(5) by variant-set key. StandardScaler + LR(C=0.5, balanced). Baseline structure vs exploratory biophysical (Grantham distance, recognition helix alpha4, helix breakers).

## Novelty claim
Reproduces the pre-registered E08 negative baseline; documents an exploratory subset hypothesis.

## Caveats / negative results
DNA chains are a modeled complex. Grantham matrix captures radical vs conservative substitutions. Exploratory subset only (n=312); balanced class weighting deviates from production calibration.

## Reproduce
```bash
python eval/e08_novel_variant.py
```

## Provenance
- **Date (UTC):** 2026-11-19T17:18:58.314813+00:00
- **Git commit:** 22016ec
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `subset.csv`
