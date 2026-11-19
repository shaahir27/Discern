# E03 · Catalogue CV (S1, patient-grouped)

> **Approach:** `catalogue` · **Test:** `e03_s1_catalogue_cv` · **PASS ✅**

Rebuilt catalogue under patient-grouped CV (stricter than the paper's shuffled CV).

## Metrics

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

## Method
catomatic BinaryBuilder per training fold (Binomial, background 0.1, two tails, strict_unlock, frs 0.1, no_mmpL5 wildcards); piezo predict on test; patient-grouped 5-fold.

## Novelty claim
Baseline reproduction of Adlard 2025 catalogue performance.

## Caveats / negative results
In-sample (built on all data) would be optimistic; this is out-of-fold.

## Reproduce
```bash
python eval/e03_catalogue_cv.py
```

## Provenance
- **Date (UTC):** 2026-11-18T19:05:18.291750+00:00
- **Git commit:** uncommitted
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `metrics.csv`
