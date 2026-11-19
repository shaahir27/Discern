# E09 · Hybrid on catalogue-U

> **Approach:** `hybrid` · **Test:** `e09_hybrid` · **PASS ✅**

EXPLORATORY: model ranks U samples (AUC ~0.74) but resolving them is unsafe (resolved-S group ~9.8% resistant; resolved-R ~40% susceptible). NEVER auto-call.

## Metrics

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `auc_ml_within_U` | 0.7442 | 0.7440 | 0.0300 | ✅ |
| `resolved_S` | 173 | 173 | 15 | ✅ |
| `resolved_S_resistant` | 17 | 17 | 15 | ✅ |
| `resolved_R` | 126 | 126 | 15 | ✅ |
| `resolved_R_susceptible` | 51 | 51 | 15 | ✅ |
| `remain_U` | 362 | 362 | 15 | ✅ |
| `remain_U_resistant` | 114 | 114 | 15 | ✅ |

| Metric | Value |
|---|---|
| `n_U` | 661 |
| `n_U_resistant` | 206 |

## Method
Among catalogue-U samples, apply per-fold model lo/hi to resolve some. ML AUC within U measures ranking value.

## Novelty claim
Queue ranking value; rejects unsafe auto-resolution (honest).

## Caveats / negative results
Used only to RANK the unresolved queue, not to change calls.

## Reproduce
```bash
python eval/e09_hybrid.py
```

## Provenance
- **Date (UTC):** 2026-11-18T19:05:39.182199+00:00
- **Git commit:** uncommitted
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `U_samples.csv`
