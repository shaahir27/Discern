# E05b Comparison

> **Approach:** `head_to_head` · **Test:** `e05b_comparison` · **DESCRIPTIVE •**

Best by AUC: rf_model (0.9094). No candidate beats golden LR by >0.005 AUC; KEEP golden LR (F2) as production.

## Metrics

| Metric | Value |
|---|---|
| `lr_model_auc` | 0.9085 |
| `lr_model_auprc` | 0.6836 |
| `lr_model_binary_sens` | 0.8349 |
| `lr_model_binary_spec` | 0.9481 |
| `lr_model_sens_at_spec_0.97` | 0.7818 |
| `lr_model_sens_at_spec_0.985` | 0.6039 |
| `lr_l1_auc` | 0.9053 |
| `lr_l1_auprc` | 0.6753 |
| `lr_l1_binary_sens` | 0.8395 |
| `lr_l1_binary_spec` | 0.9455 |
| `lr_l1_sens_at_spec_0.97` | 0.7898 |
| `lr_l1_sens_at_spec_0.985` | 0.5993 |
| `rf_model_auc` | 0.9094 |
| `rf_model_auprc` | 0.6818 |
| `rf_model_binary_sens` | 0.8337 |
| `rf_model_binary_spec` | 0.9548 |
| `rf_model_sens_at_spec_0.97` | 0.8048 |
| `rf_model_sens_at_spec_0.985` | 0.6594 |
| `hgb_model_auc` | 0.9052 |
| `hgb_model_auprc` | 0.6787 |
| `hgb_model_binary_sens` | 0.8372 |
| `hgb_model_binary_spec` | 0.9528 |
| `hgb_model_sens_at_spec_0.97` | 0.8141 |
| `hgb_model_sens_at_spec_0.985` | 0.6259 |
| `svm_model_auc` | 0.8953 |
| `svm_model_auprc` | 0.6829 |
| `svm_model_binary_sens` | 0.7979 |
| `svm_model_binary_spec` | 0.9765 |
| `svm_model_sens_at_spec_0.97` | 0.8014 |
| `svm_model_sens_at_spec_0.985` | 0.6732 |

## Method
Same-data model-class comparison on identical 11,868 samples, F2 features, S1 patient-grouped nested CV. Ranked by AUC; clinical operating point sens@spec0.985 also reported.

## Novelty claim
Justifies the model class choice (calibrated LR) with a same-data ablation.

## Caveats / negative results
Exploratory; hyperparameters fixed a priori (D007), not tuned. match != beat. AUC differences <0.005 are within noise.

## Reproduce
```bash
python eval/e05b_model_comparison.py
```

## Provenance
- **Date (UTC):** 2026-11-19T04:05:25.414443+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `comparison.csv`
