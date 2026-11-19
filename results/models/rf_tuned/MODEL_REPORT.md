# Tier B — tuned random forest

> **Type:** model · **Approach:** `rf_tuned` · **0/2 tests pass** ❌

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E21 · KM External | — | auc=0.8675 | [report](e21_km_external/REPORT.md) |
| E21 · Tuned CV | — | auc=0.9103 | [report](e21_tuned_cv/REPORT.md) |

## E21 · KM External — DESCRIPTIVE •

| Metric | Value |
|---|---|
| `auc` | 0.8675 |
| `binary_sens` | 0.8286 |
| `binary_spec` | 0.7860 |
| `sens_at_spec_0.95` | 0.7714 |
| `n_isolates` | 792 |
| `n_resistant` | 35 |

- **zero_training_overlap:** `True`

_Caveats: Research use only. The selected settings and Platt calibration use the same inner folds. The outer and KM results are held-out evaluations. KM is a sensitivity analysis with undocumented assay and critical concentration, medium-quality phenotypes, and few resistant isolates. Sensitivity at fixed specificity describes the held-out ROC curve._

## E21 · Tuned CV — DESCRIPTIVE •

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

_Caveats: Research use only. The selected settings and Platt calibration use the same inner folds. The outer and KM results are held-out evaluations._

