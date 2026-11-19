# E13 · PPV/NPV at prevalence

> **Approach:** `predictive_values` · **Test:** `e13_ppv_npv` · **PASS ✅**

At 2% prevalence ~half of catalogue R calls are false alarms (PPV 50%). R is a flag for confirmation; S is reliable (NPV ~99.6%).

## Metrics

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `model_ppv_0.01` | 0.1436 | 0.1400 | 0.0200 | ✅ |
| `model_npv_0.01` | 0.9982 | 0.9980 | 0.0030 | ✅ |
| `cat_ppv_0.01` | 0.3328 | 0.3300 | 0.0200 | ✅ |
| `cat_npv_0.01` | 0.9978 | 0.9980 | 0.0030 | ✅ |
| `model_ppv_0.02` | 0.2530 | 0.2500 | 0.0200 | ✅ |
| `model_npv_0.02` | 0.9964 | 0.9960 | 0.0030 | ✅ |
| `cat_ppv_0.02` | 0.5019 | 0.5000 | 0.0200 | ✅ |
| `cat_npv_0.02` | 0.9957 | 0.9960 | 0.0030 | ✅ |
| `model_ppv_0.05` | 0.4663 | 0.4700 | 0.0200 | ✅ |
| `model_npv_0.05` | 0.9907 | 0.9910 | 0.0030 | ✅ |
| `cat_ppv_0.05` | 0.7221 | 0.7200 | 0.0200 | ✅ |
| `cat_npv_0.05` | 0.9889 | 0.9890 | 0.0030 | ✅ |
| `model_ppv_0.1` | 0.6484 | 0.6500 | 0.0200 | ✅ |
| `model_npv_0.1` | 0.9805 | 0.9810 | 0.0030 | ✅ |
| `cat_ppv_0.1` | 0.8458 | 0.8500 | 0.0200 | ✅ |
| `cat_npv_0.1` | 0.9768 | 0.9770 | 0.0030 | ✅ |

## Method
Bayes' rule at prevalence 1/2/5/10% using CV sens/spec: model (0.83/0.95), catalogue (0.79/0.984).

## Novelty claim
Prevalence-aware predictive values (clinical framing).

## Caveats / negative results
Illustrative operating points; real lab prevalence must replace the placeholder.

## Reproduce
```bash
python eval/e13_ppv_npv.py
```

## Provenance
- **Date (UTC):** 2026-11-18T19:05:51.853935+00:00
- **Git commit:** uncommitted
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `ppv_npv.csv`
