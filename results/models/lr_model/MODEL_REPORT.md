# Tier B — model (LR F2)

> **Type:** model · **Approach:** `lr_model` · **5/5 tests pass** ✅

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E04 · Model CV (S1, patient-grouped) | ✅ | auc=0.9085 | [report](e04_s1_patient_cv/REPORT.md) |
| E05 · Feature ablation | ✅ | F0_frs0.1_auc=0.7910 | [report](e05_ablation/REPORT.md) |
| E06 · Cross-source (S3) | ✅ | CRyPTIC_to_NICD_auc=0.9105 | [report](e06_cross_source/REPORT.md) |
| E07 · Leave-one-site-out (S2) | ✅ | pooled_sens=0.7875 | [report](e07_loso/REPORT.md) |
| E11 · Imbalance + learning curve | ✅ | lc_auc_0.1=0.8963 | [report](e11_imbalance_learning_curve/REPORT.md) |

## E04 · Model CV (S1, patient-grouped) — PASS ✅

Calibrated LR matches the catalogue at matched operating points; it does not beat it. Adds ranking + calibrated probability.

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `auc` | 0.9085 | 0.9090 | 0.0100 | ✅ |
| `auprc` | 0.6836 | 0.6840 | 0.0200 | ✅ |
| `binary_sens` | 0.8349 | 0.8350 | 0.0200 | ✅ |
| `binary_spec` | 0.9481 | 0.9480 | 0.0100 | ✅ |
| `ternary_firm` | 0.9469 | 0.9470 | 0.0100 | ✅ |
| `sens_at_spec_0.95` | 0.8303 | 0.8300 | 0.0300 | ✅ |
| `sens_at_spec_0.97` | 0.7818 | 0.7820 | 0.0300 | ✅ |
| `sens_at_spec_0.98` | 0.6524 | 0.6520 | 0.0300 | ✅ |
| `sens_at_spec_0.985` | 0.6039 | 0.6040 | 0.0300 | ✅ |

| Metric | Value |
|---|---|
| `ternary_sens_among_firm` | 0.7347 |
| `ternary_spec_among_firm` | 0.9889 |
| `n_features_all_data` | 120 |

_Caveats: ternary sens/spec among firm (0.735/0.989) are descriptive (no tolerance). Model != better accuracy._

## E05 · Feature ablation — PASS ✅

Type/region features ~+0.11 AUC; minor alleles ~+0.11 AUC; mmpL5 interaction ~+0.004 AUC.

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `F0_frs0.1_auc` | 0.7910 | 0.7910 | 0.0150 | ✅ |
| `F1_frs0.1_auc` | 0.9045 | 0.9050 | 0.0150 | ✅ |
| `F2_frs0.1_auc` | 0.9085 | 0.9090 | 0.0150 | ✅ |
| `F2_frs0.75_auc` | 0.7986 | 0.7990 | 0.0150 | ✅ |
| `F0_frs0.1_auprc` | 0.5026 | 0.5030 | 0.0150 | ✅ |
| `F1_frs0.1_auprc` | 0.6721 | 0.6720 | 0.0150 | ✅ |
| `F2_frs0.1_auprc` | 0.6836 | 0.6840 | 0.0150 | ✅ |
| `F2_frs0.75_auprc` | 0.5479 | 0.5480 | 0.0150 | ✅ |

_Caveats: Ablation, not tuning. frs 0.75 drops minor-allele signal._

## E06 · Cross-source (S3) — PASS ✅

Accuracy depends on lab/source: NICD->CRyPTIC drops (few CRyPTIC resistant training germs, lineage/prevalence shift).

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `CRyPTIC_to_NICD_auc` | 0.9105 | 0.9100 | 0.0300 | ✅ |
| `CRyPTIC_to_NICD_sens` | 0.8442 | 0.8440 | 0.0300 | ✅ |
| `CRyPTIC_to_NICD_spec` | 0.9633 | 0.9630 | 0.0300 | ✅ |
| `NICD_to_CRyPTIC_auc` | 0.8277 | 0.8300 | 0.0300 | ✅ |
| `NICD_to_CRyPTIC_sens` | 0.6747 | 0.6750 | 0.0300 | ✅ |
| `NICD_to_CRyPTIC_spec` | 0.9532 | 0.9530 | 0.0300 | ✅ |

_Caveats: NICD test has 783 R; CRyPTIC test only 83 R (wide CI)._

## E07 · Leave-one-site-out (S2) — PASS ✅

Pooled sens/spec across held-out sites. Per-site specificity stable; per-site sensitivity unreliable (most sites <25 R).

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `pooled_sens` | 0.7875 | 0.7880 | 0.0200 | ✅ |
| `pooled_spec` | 0.9537 | 0.9540 | 0.0100 | ✅ |

| Metric | Value |
|---|---|
| `per_site_spec_min_n300` | 0.9441 |
| `per_site_spec_max_n300` | 0.9785 |

_Caveats: Do NOT report pooled LOSO AUC (mixes 23 calibrations). Site 5 (inferred India, unverified) has ~5 R -> false-alarm audit only._

## E11 · Imbalance + learning curve — PASS ✅

Balancing does not help and breaks calibration; AUC is flat across resistant fractions (signal saturates early).

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `lc_auc_0.1` | 0.8963 | 0.8980 | 0.0100 | ✅ |
| `lc_auc_0.25` | 0.9050 | 0.9050 | 0.0100 | ✅ |
| `lc_auc_0.5` | 0.9050 | 0.9080 | 0.0100 | ✅ |
| `lc_auc_0.75` | 0.9075 | 0.9060 | 0.0100 | ✅ |
| `lc_auc_1.0` | 0.9085 | 0.9090 | 0.0100 | ✅ |
| `balanced_auc` | 0.9051 | 0.9040 | 0.0100 | ✅ |
| `balanced_auprc` | 0.6587 | 0.6530 | 0.0200 | ✅ |
| `balanced_sens_at_spec98` | 0.6062 | 0.5920 | 0.0300 | ✅ |

_Caveats: Balanced model is diagnostic only; production model uses no class_weight (rule 4)._

