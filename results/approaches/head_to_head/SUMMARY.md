# Head-to-head comparisons

> **Type:** comparison · **Approach:** `head_to_head` · **0/4 tests pass** ❌

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E05b Comparison | — | lr_model_auc=0.9085 | [report](e05b_comparison/REPORT.md) |
| E21 · Comparison | — |  | [report](e21_comparison/REPORT.md) |
| G7 · Adlard | — | shuffled_firm=0.9449 | [report](g7_adlard/REPORT.md) |
| G7 · WHO V2 | — | who_frs0.1_firm=0.9607 | [report](g7_who_v2/REPORT.md) |

## E05b Comparison — DESCRIPTIVE •

Best by AUC: rf_model (0.9094). No candidate beats golden LR by >0.005 AUC; KEEP golden LR (F2) as production.

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

_Caveats: Exploratory; hyperparameters fixed a priori (D007), not tuned. match != beat. AUC differences <0.005 are within noise._

## E21 · Comparison — DESCRIPTIVE •

Keep golden LR (F2). No tuned candidate met both D013 criteria.

- **S1_leader:** `rf_tuned`
- **KM_leader:** `['hgb_tuned']`

_Caveats: Research use only. The selected settings and Platt calibration use the same inner folds. The outer and KM results are held-out evaluations. KM is a sensitivity analysis with undocumented assay and critical concentration, medium-quality phenotypes, and few resistant isolates. Sensitivity at fixed specificity describes the held-out ROC curve._

## G7 · Adlard — DESCRIPTIVE •

Shuffled sens_firm 0.7955 vs patient-grouped 0.7924 (delta +0.0031). Patient grouping is the stricter, leakage-controlled protocol. Ghosh comparison infeasible (data not local).

| Metric | Value |
|---|---|
| `shuffled_firm` | 0.9449 |
| `shuffled_sens_firm` | 0.7955 |
| `shuffled_spec_firm` | 0.9838 |
| `shuffled_n_U` | 654 |
| `grouped_firm` | 0.9443 |
| `grouped_sens_firm` | 0.7924 |
| `grouped_spec_firm` | 0.9835 |
| `grouped_n_U` | 661 |
| `sens_firm_shuffled_minus_grouped` | 0.0031 |

- **ghosh_reconstructable:** `False`

_Caveats: Shuffled CV can leak related germs across folds (optimistic). Grouped controls patient duplicates but not clones/outbreaks fully. Ghosh omitted honestly; NEVER cross-cohort AUC._

## G7 · WHO V2 — DESCRIPTIVE •

WHO v2 ternary sens_firm at FRS0.25 = 0.7680; our model ternary sens among firm = 0.7347. Our ternary does NOT beat WHO at matched FRS (honest). Binary U-as-S sensitivities reported at each FRS convention.

| Metric | Value |
|---|---|
| `who_frs0.1_firm` | 0.9607 |
| `who_frs0.1_sens_firm` | 0.8011 |
| `who_frs0.1_spec_firm` | 0.9812 |
| `who_frs0.1_sens_U_as_S` | 0.6790 |
| `who_frs0.1_sens_U_as_R` | 0.8314 |
| `who_frs0.25_firm` | 0.9611 |
| `who_frs0.25_sens_firm` | 0.7680 |
| `who_frs0.25_spec_firm` | 0.9828 |
| `who_frs0.25_sens_U_as_S` | 0.6420 |
| `who_frs0.25_sens_U_as_R` | 0.8060 |
| `who_frs0.75_firm` | 0.9680 |
| `who_frs0.75_sens_firm` | 0.5606 |
| `who_frs0.75_spec_firm` | 0.9868 |
| `who_frs0.75_sens_U_as_S` | 0.4965 |
| `who_frs0.75_sens_U_as_R` | 0.6109 |
| `our_catalogue_grouped_sens_firm` | 0.7924 |
| `our_catalogue_grouped_spec_firm` | 0.9835 |
| `our_model_ternary_sens_among_firm` | 0.7347 |
| `our_model_ternary_spec_among_firm` | 0.9889 |

_Caveats: Same-data comparison only; NEVER compare AUC across cohorts. Ternary sens among firm is descriptive. WHO FRS convention matters (0.25 vs 0.1 vs 0.75). match != beat._

