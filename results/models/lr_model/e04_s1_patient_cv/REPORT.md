# E04 · Model CV (S1, patient-grouped)

> **Approach:** `lr_model` · **Test:** `e04_s1_patient_cv` · **PASS ✅**

Calibrated LR matches the catalogue at matched operating points; it does not beat it. Adds ranking + calibrated probability.

## Metrics

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `auc` | 0.9085 | 0.9090 | 0.0100 | ✅ |
| `auprc` | 0.6836 | 0.6840 | 0.0200 | ✅ |
| `binary_sens` | 0.8349 | 0.8350 | 0.0200 | ✅ |
| `binary_spec` | 0.9481 | 0.9480 | 0.0100 | ✅ |
| `ternary_firm` | 0.9469 | 0.9470 | 0.0100 | ✅ |
| `sens_at_spec_0.95` | 0.8303 | 0.8300 | 0.0300 | ✅ |
| `sens_at_spec_0.97` | 0.7818 | 0.7820 | 0.0300 | ✅ |
| `sens_at_spec_0.98` | 0.6524 | 0.6520 | 0.0300 | ✅ |
| `sens_at_spec_0.985` | 0.6039 | 0.6040 | 0.0300 | ✅ |

| Metric | Value |
|---|---|
| `ternary_sens_among_firm` | 0.7347 |
| `ternary_spec_among_firm` | 0.9889 |
| `n_features_all_data` | 120 |

## Method
L2 LogisticRegression C=0.3; nested patient-grouped Platt calibration; thresholds from inner OOF (NPV>=0.985, PPV>=0.80, Youden binary); F2 features; vocabulary>=3 training carriers.

## Novelty claim
Baseline reproduction; model adds calibration/ranking, not superior accuracy.

## Caveats / negative results
ternary sens/spec among firm (0.735/0.989) are descriptive (no tolerance). Model != better accuracy.

## Reproduce
```bash
python eval/e03_catalogue_cv.py && python eval/e04_lr_f2.py
```

## Provenance
- **Date (UTC):** 2026-11-18T19:05:21.788844+00:00
- **Git commit:** uncommitted
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `metrics.csv`
- `thresholds_per_fold.csv`
