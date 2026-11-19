# E21 · KM External

> **Approach:** `rf_tuned` · **Test:** `e21_km_external` · **DESCRIPTIVE •**

## Metrics

| Metric | Value |
|---|---|
| `auc` | 0.8675 |
| `binary_sens` | 0.8286 |
| `binary_spec` | 0.7860 |
| `sens_at_spec_0.95` | 0.7714 |
| `n_isolates` | 792 |
| `n_resistant` | 35 |

- **zero_training_overlap:** `True`

## Method
F2 with a vocabulary of variants carried by at least three training samples. Select by ROC AUC on three patient-grouped inner folds, then fit Platt scaling on raw out-of-fold scores from those folds and derive the golden thresholds from the calibrated scores. Refit using the entire golden cohort and evaluate on KMVALIDATION2025.

## Caveats / negative results
Research use only. The selected settings and Platt calibration use the same inner folds. The outer and KM results are held-out evaluations. KM is a sensitivity analysis with undocumented assay and critical concentration, medium-quality phenotypes, and few resistant isolates. Sensitivity at fixed specificity describes the held-out ROC curve.

## Reproduce
```bash
C:/Users/shaah/miniconda3/envs/discern/python.exe eval/e21_tuned_comparison.py
```

## Provenance
- **Date (UTC):** 2026-11-19T05:29:23.006385+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `metrics.csv`
- `best_params.csv`
