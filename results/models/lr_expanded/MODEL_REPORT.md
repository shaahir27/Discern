# Tier B — expanded-training (CRyPTIC v3.4)

> **Type:** model · **Approach:** `lr_expanded` · **0/4 tests pass** ❌

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E18 · Cross Source | — | CRyPTIC_to_NICD_auc=0.9020 | [report](e18_cross_source/REPORT.md) |
| E18 · Expanded CV | — | auc=0.8381 | [report](e18_expanded_cv/REPORT.md) |
| E18 · KM External | — | n_isolates=792 | [report](e18_km_external/REPORT.md) |
| E18b High Only | — | auc=0.6668 | [report](e18b_high_only/REPORT.md) |

## E18 · Cross Source — DESCRIPTIVE •

Cross-source generalization with richer CRyPTIC resistant data (255 R vs 83). Compare NICD->CRyPTIC sensitivity to golden 0.675 (the weak spot).

| Metric | Value |
|---|---|
| `CRyPTIC_to_NICD_auc` | 0.9020 |
| `CRyPTIC_to_NICD_sens` | 0.8557 |
| `CRyPTIC_to_NICD_spec` | 0.8707 |
| `CRyPTIC_to_NICD_golden_auc` | 0.9105 |
| `NICD_to_CRyPTIC_auc` | 0.6533 |
| `NICD_to_CRyPTIC_sens` | 0.3320 |
| `NICD_to_CRyPTIC_spec` | 0.9667 |
| `NICD_to_CRyPTIC_golden_auc` | 0.8277 |

_Caveats: Descriptive; different cohort. CRyPTIC test R count larger than golden (83). No tuning._

## E18 · Expanded CV — DESCRIPTIVE •

Expanded training AUC 0.8381 vs golden 0.9085 (delta -0.0704); sens@spec0.985 0.6195. +182 new resistant isolates, CRyPTIC R 83->255.

| Metric | Value |
|---|---|
| `n_samples` | 23441 |
| `n_resistant` | 1046 |
| `auc` | 0.8381 |
| `auprc` | 0.5665 |
| `binary_sens` | 0.7247 |
| `binary_spec` | 0.9650 |
| `sens_at_spec_0.95` | 0.7266 |
| `sens_at_spec_0.97` | 0.7132 |
| `sens_at_spec_0.985` | 0.6195 |
| `n_features` | 137 |
| `golden_auc` | 0.9085 |
| `delta_auc_vs_golden` | -0.0704 |

_Caveats: Different cohort than golden (not the published 11,868 reproduction). MEDIUM labels add noise. Descriptive; golden baseline + production untouched; no tuning. KM held out._

## E18 · KM External — DESCRIPTIVE •

Expanded model external KM AUC 0.8636 vs golden production 0.8682 (delta -0.0046). Tests whether richer training improves generalization to unseen isolates.

| Metric | Value |
|---|---|
| `n_isolates` | 792 |
| `n_resistant` | 35 |
| `expanded_model_auc` | 0.8636 |
| `expanded_model_sens` | 0.8571 |
| `expanded_model_spec` | 0.7728 |
| `expanded_sens_at_spec_0.95` | 0.6571 |
| `golden_production_km_auc` | 0.8682 |
| `delta_auc_vs_golden_production` | -0.0046 |

_Caveats: KM assay undocumented (MEDIUM); 35 R -> wide CI; sensitivity analysis. No tuning; KM excluded from expanded training. match != beat._

## E18b High Only — DESCRIPTIVE •

HIGH-only expansion (16923/201 R): S1 AUC 0.6668 (golden 0.9085, HIGH+MEDIUM 0.8381); KM external AUC 0.8633 (golden production 0.8682). Cross-source: N/A (NICD is MEDIUM, excluded from HIGH-only).

| Metric | Value |
|---|---|
| `n_samples` | 16923 |
| `n_resistant` | 201 |
| `auc` | 0.6668 |
| `binary_sens` | 0.4030 |
| `binary_spec` | 0.9637 |
| `sens_at_spec_0.95` | 0.4080 |
| `sens_at_spec_0.985` | 0.2687 |
| `golden_auc` | 0.9085 |
| `high_medium_auc` | 0.8381 |
| `km_external_auc` | 0.8633 |
| `golden_production_km_auc` | 0.8682 |

- **cross_source_note:** `NICD is MEDIUM quality; absent from HIGH-only arm -> cross-source N/A`

_Caveats: Descriptive; different cohort than golden reproduction. KM assay undocumented; 35 R wide CI. Golden + production untouched; no tuning. Tests if CLEAN richer data helps generalization._

