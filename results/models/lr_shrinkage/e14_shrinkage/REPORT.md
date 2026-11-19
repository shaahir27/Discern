# E14 · Shrinkage

> **Approach:** `lr_shrinkage` · **Test:** `e14_shrinkage` · **DESCRIPTIVE •**

1042 variants shrunk; mean |shrunk-raw| 0.1432; 48 variants flip a 0.5 grade. Rare variants pulled toward class mean; standard estimator choice, not an invention.

## Metrics

| Metric | Value |
|---|---|
| `n_variants` | 1042 |
| `mean_abs_shrinkage` | 0.1432 |
| `max_abs_shrinkage` | 0.9294 |
| `n_grade_flips_0.5` | 48 |
| `n_mech_classes` | 9 |

## Method
Empirical-Bayes beta-binomial shrinkage of per-variant R proportion toward mechanism-class mean (method-of-moments Beta prior). Wilson 95% CI on raw proportion.

## Novelty claim
M3 shrinkage rung (beta-binomial mechanism-class evidence grading).

## Caveats / negative results
Descriptive estimator output for evidence cards; ablate, don't headline. No accuracy claim. Class priors estimated from the same data (empirical Bayes).

## Reproduce
```bash
python eval/e14_conformal.py
```

## Provenance
- **Date (UTC):** 2026-11-19T04:25:38.551123+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `variant_shrinkage.csv`
- `class_priors.csv`
