# E05b Model CV

> **Approach:** `svm_model` · **Test:** `e05b_model_cv` · **DESCRIPTIVE •**

Exploratory model-class candidate. AUC 0.8953, sens@spec0.985 0.6732 vs golden LR 0.9085/0.6039.

## Metrics

| Metric | Value |
|---|---|
| `auc` | 0.8953 |
| `auprc` | 0.6829 |
| `binary_sens` | 0.7979 |
| `binary_spec` | 0.9765 |
| `ternary_firm` | 0.0709 |
| `sens_at_spec_0.95` | 0.8014 |
| `sens_at_spec_0.97` | 0.8014 |
| `sens_at_spec_0.98` | 0.7552 |
| `sens_at_spec_0.985` | 0.6732 |

## Method
linear SVM C=0.3 (decision_function -> Platt). Identical F2 features, S1 patient-grouped 5-fold, nested patient-grouped Platt calibration + thresholds (NPV>=0.985, PPV>=0.8, Youden).

## Novelty claim
Model-class ablation: does a non-linear model beat calibrated LR on identical features/splits?

## Caveats / negative results
No pre-registered acceptance numbers; hyperparameters fixed a priori (D007), not tuned. Descriptive comparison only.

## Reproduce
```bash
python eval/e05b_model_comparison.py
```

## Provenance
- **Date (UTC):** 2026-11-19T04:05:25.310623+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `metrics.csv`
- `thresholds_per_fold.csv`
