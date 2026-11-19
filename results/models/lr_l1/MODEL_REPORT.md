# Tier B — L1 logistic regression

> **Type:** model · **Approach:** `lr_l1` · **0/1 tests pass** ❌

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E05b Model CV | — | auc=0.9053 | [report](e05b_model_cv/REPORT.md) |

## E05b Model CV — DESCRIPTIVE •

Exploratory model-class candidate. AUC 0.9053, sens@spec0.985 0.5993 vs golden LR 0.9085/0.6039.

| Metric | Value |
|---|---|
| `auc` | 0.9053 |
| `auprc` | 0.6753 |
| `binary_sens` | 0.8395 |
| `binary_spec` | 0.9455 |
| `ternary_firm` | 0.9414 |
| `sens_at_spec_0.95` | 0.8268 |
| `sens_at_spec_0.97` | 0.7898 |
| `sens_at_spec_0.98` | 0.6651 |
| `sens_at_spec_0.985` | 0.5993 |

_Caveats: No pre-registered acceptance numbers; hyperparameters fixed a priori (D007), not tuned. Descriptive comparison only._

