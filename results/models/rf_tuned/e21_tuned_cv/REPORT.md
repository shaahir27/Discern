# E21 · Tuned CV

> **Approach:** `rf_tuned` · **Test:** `e21_tuned_cv` · **DESCRIPTIVE •**

## Metrics

| Metric | Value |
|---|---|
| `auc` | 0.9103 |
| `binary_sens` | 0.8314 |
| `binary_spec` | 0.9580 |
| `sens_at_spec_0.95` | 0.8406 |
| `sens_at_spec_0.97` | 0.8152 |
| `sens_at_spec_0.985` | 0.6940 |
| `auprc` | 0.7011 |
| `ternary_firm` | 0.9736 |

## Method
F2 with a vocabulary of variants carried by at least three training samples. Select by ROC AUC on three patient-grouped inner folds, then fit Platt scaling on raw out-of-fold scores from those folds and derive the golden thresholds from the calibrated scores. Evaluation uses the five saved S1 patient-grouped test folds.

## Caveats / negative results
Research use only. The selected settings and Platt calibration use the same inner folds. The outer and KM results are held-out evaluations.

## Reproduce
```bash
C:/Users/shaah/miniconda3/envs/discern/python.exe eval/e21_tuned_comparison.py
```

## Provenance
- **Date (UTC):** 2026-11-19T05:28:58.089421+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `metrics.csv`
- `best_params_per_fold.csv`
