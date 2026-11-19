# Tier B — histogram gradient boosting

> **Type:** model · **Approach:** `hgb_model` · **0/1 tests pass** ❌

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E05b Model CV | — | auc=0.9052 | [report](e05b_model_cv/REPORT.md) |

## E05b Model CV — DESCRIPTIVE •

Exploratory model-class candidate. AUC 0.9052, sens@spec0.985 0.6259 vs golden LR 0.9085/0.6039.

| Metric | Value |
|---|---|
| `auc` | 0.9052 |
| `auprc` | 0.6787 |
| `binary_sens` | 0.8372 |
| `binary_spec` | 0.9528 |
| `ternary_firm` | 0.9549 |
| `sens_at_spec_0.95` | 0.8395 |
| `sens_at_spec_0.97` | 0.8141 |
| `sens_at_spec_0.98` | 0.7194 |
| `sens_at_spec_0.985` | 0.6259 |

_Caveats: No pre-registered acceptance numbers; hyperparameters fixed a priori (D007), not tuned. Descriptive comparison only._

