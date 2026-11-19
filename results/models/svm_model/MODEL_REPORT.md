# Tier B — linear SVM

> **Type:** model · **Approach:** `svm_model` · **0/1 tests pass** ❌

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E05b Model CV | — | auc=0.8953 | [report](e05b_model_cv/REPORT.md) |

## E05b Model CV — DESCRIPTIVE •

Exploratory model-class candidate. AUC 0.8953, sens@spec0.985 0.6732 vs golden LR 0.9085/0.6039.

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

_Caveats: No pre-registered acceptance numbers; hyperparameters fixed a priori (D007), not tuned. Descriptive comparison only._

