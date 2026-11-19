# Tier B — random forest

> **Type:** model · **Approach:** `rf_model` · **0/1 tests pass** ❌

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E05b Model CV | — | auc=0.9094 | [report](e05b_model_cv/REPORT.md) |

## E05b Model CV — DESCRIPTIVE •

Exploratory model-class candidate. AUC 0.9094, sens@spec0.985 0.6594 vs golden LR 0.9085/0.6039.

| Metric | Value |
|---|---|
| `auc` | 0.9094 |
| `auprc` | 0.6818 |
| `binary_sens` | 0.8337 |
| `binary_spec` | 0.9548 |
| `ternary_firm` | 0.9643 |
| `sens_at_spec_0.95` | 0.8337 |
| `sens_at_spec_0.97` | 0.8048 |
| `sens_at_spec_0.98` | 0.7425 |
| `sens_at_spec_0.985` | 0.6594 |

_Caveats: No pre-registered acceptance numbers; hyperparameters fixed a priori (D007), not tuned. Descriptive comparison only._

