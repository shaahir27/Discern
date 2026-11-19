# Evaluation results

Research use only. For evaluation and education. Clinical decisions require laboratory confirmation.

Metrics are reproduced from the indexed result JSON.

## hgb_model/e05b_model_cv

Descriptive · [Source JSON](../results/models/hgb_model/e05b_model_cv/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `auc` | 0.9052400403369869 | — | — |
| `auprc` | 0.678652962897792 | — | — |
| `binary_sens` | 0.8371824480369515 | — | — |
| `binary_spec` | 0.9528267587711325 | — | — |
| `ternary_firm` | 0.9549207954162454 | — | — |
| `sens_at_spec_0.95` | 0.8394919168591224 | — | — |
| `sens_at_spec_0.97` | 0.8140877598152425 | — | — |
| `sens_at_spec_0.98` | 0.7193995381062356 | — | — |
| `sens_at_spec_0.985` | 0.625866050808314 | — | — |

[metrics.csv](../results/models/hgb_model/e05b_model_cv/metrics.csv) · [thresholds_per_fold.csv](../results/models/hgb_model/e05b_model_cv/thresholds_per_fold.csv)

## hgb_tuned/e21_km_external

Descriptive · [Source JSON](../results/models/hgb_tuned/e21_km_external/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `auc` | 0.8892243819588602 | — | — |
| `binary_sens` | 0.8571428571428571 | — | — |
| `binary_spec` | 0.7793923381770145 | — | — |
| `sens_at_spec_0.95` | 0.7714285714285715 | — | — |
| `n_isolates` | 792 | — | — |
| `n_resistant` | 35 | — | — |
| `zero_training_overlap` | True | — | — |

[metrics.csv](../results/models/hgb_tuned/e21_km_external/metrics.csv) · [best_params.csv](../results/models/hgb_tuned/e21_km_external/best_params.csv)

## hgb_tuned/e21_tuned_cv

Descriptive · [Source JSON](../results/models/hgb_tuned/e21_tuned_cv/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `auc` | 0.9057513897326247 | — | — |
| `binary_sens` | 0.8325635103926097 | — | — |
| `binary_spec` | 0.9536447918560261 | — | — |
| `sens_at_spec_0.95` | 0.8383371824480369 | — | — |
| `sens_at_spec_0.97` | 0.8140877598152425 | — | — |
| `sens_at_spec_0.985` | 0.628175519630485 | — | — |
| `auprc` | 0.6753630222453036 | — | — |
| `ternary_firm` | 0.9590495449949444 | — | — |

[metrics.csv](../results/models/hgb_tuned/e21_tuned_cv/metrics.csv) · [best_params_per_fold.csv](../results/models/hgb_tuned/e21_tuned_cv/best_params_per_fold.csv)

## lr_expanded/e18_cross_source

Descriptive · [Source JSON](../results/models/lr_expanded/e18_cross_source/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `CRyPTIC_to_NICD_auc` | 0.9020478034216548 | — | — |
| `CRyPTIC_to_NICD_sens` | 0.855683269476373 | — | — |
| `CRyPTIC_to_NICD_spec` | 0.8707124010554089 | — | — |
| `CRyPTIC_to_NICD_golden_auc` | 0.9104980172497594 | — | — |
| `NICD_to_CRyPTIC_auc` | 0.6532820593313453 | — | — |
| `NICD_to_CRyPTIC_sens` | 0.332 | — | — |
| `NICD_to_CRyPTIC_spec` | 0.9666718987076859 | — | — |
| `NICD_to_CRyPTIC_golden_auc` | 0.8277106164985712 | — | — |

[cross_source.csv](../results/models/lr_expanded/e18_cross_source/cross_source.csv)

## lr_expanded/e18_expanded_cv

Descriptive · [Source JSON](../results/models/lr_expanded/e18_expanded_cv/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `n_samples` | 23441 | — | — |
| `n_resistant` | 1046 | — | — |
| `auc` | 0.838130054125541 | — | — |
| `auprc` | 0.5664602456715522 | — | — |
| `binary_sens` | 0.7246653919694073 | — | — |
| `binary_spec` | 0.9650368385800402 | — | — |
| `sens_at_spec_0.95` | 0.7265774378585086 | — | — |
| `sens_at_spec_0.97` | 0.7131931166347992 | — | — |
| `sens_at_spec_0.985` | 0.6195028680688337 | — | — |
| `n_features` | 137 | — | — |
| `golden_auc` | 0.9085389366535497 | — | — |
| `delta_auc_vs_golden` | -0.0704088825280087 | — | — |

[metrics.csv](../results/models/lr_expanded/e18_expanded_cv/metrics.csv) · [thresholds_per_fold.csv](../results/models/lr_expanded/e18_expanded_cv/thresholds_per_fold.csv)

## lr_expanded/e18_km_external

Descriptive · [Source JSON](../results/models/lr_expanded/e18_km_external/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `n_isolates` | 792 | — | — |
| `n_resistant` | 35 | — | — |
| `expanded_model_auc` | 0.8636157765616155 | — | — |
| `expanded_model_sens` | 0.8571428571428571 | — | — |
| `expanded_model_spec` | 0.7727873183619551 | — | — |
| `expanded_sens_at_spec_0.95` | 0.6571428571428571 | — | — |
| `golden_production_km_auc` | 0.8681826759765995 | — | — |
| `delta_auc_vs_golden_production` | -0.004566899414984005 | — | — |

[metrics.csv](../results/models/lr_expanded/e18_km_external/metrics.csv)

## lr_expanded/e18b_high_only

Descriptive · [Source JSON](../results/models/lr_expanded/e18b_high_only/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `n_samples` | 16923 | — | — |
| `n_resistant` | 201 | — | — |
| `auc` | 0.6668331289373013 | — | — |
| `binary_sens` | 0.40298507462686567 | — | — |
| `binary_spec` | 0.9637005142925488 | — | — |
| `sens_at_spec_0.95` | 0.4079601990049751 | — | — |
| `sens_at_spec_0.985` | 0.26865671641791045 | — | — |
| `golden_auc` | 0.908539 | — | — |
| `high_medium_auc` | 0.8381 | — | — |
| `cross_source_note` | NICD is MEDIUM quality; absent from HIGH-only arm -> cross-source N/A | — | — |
| `km_external_auc` | 0.8633327042838272 | — | — |
| `golden_production_km_auc` | 0.8681826759765995 | — | — |

[metrics.csv](../results/models/lr_expanded/e18b_high_only/metrics.csv)

## lr_l1/e05b_model_cv

Descriptive · [Source JSON](../results/models/lr_l1/e05b_model_cv/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `auc` | 0.9053351836512614 | — | — |
| `auprc` | 0.6752893399877362 | — | — |
| `binary_sens` | 0.8394919168591224 | — | — |
| `binary_spec` | 0.9454644610070896 | — | — |
| `ternary_firm` | 0.9413549039433772 | — | — |
| `sens_at_spec_0.95` | 0.8267898383371824 | — | — |
| `sens_at_spec_0.97` | 0.789838337182448 | — | — |
| `sens_at_spec_0.98` | 0.6651270207852193 | — | — |
| `sens_at_spec_0.985` | 0.5993071593533488 | — | — |

[metrics.csv](../results/models/lr_l1/e05b_model_cv/metrics.csv) · [thresholds_per_fold.csv](../results/models/lr_l1/e05b_model_cv/thresholds_per_fold.csv)

## lr_l1_tuned/e21_km_external

Descriptive · [Source JSON](../results/models/lr_l1_tuned/e21_km_external/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `auc` | 0.8690130213247783 | — | — |
| `binary_sens` | 0.8285714285714286 | — | — |
| `binary_spec` | 0.7807133421400264 | — | — |
| `sens_at_spec_0.95` | 0.6571428571428571 | — | — |
| `n_isolates` | 792 | — | — |
| `n_resistant` | 35 | — | — |
| `zero_training_overlap` | True | — | — |

[metrics.csv](../results/models/lr_l1_tuned/e21_km_external/metrics.csv) · [best_params.csv](../results/models/lr_l1_tuned/e21_km_external/best_params.csv)

## lr_l1_tuned/e21_tuned_cv

Descriptive · [Source JSON](../results/models/lr_l1_tuned/e21_tuned_cv/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `auc` | 0.9085739397371798 | — | — |
| `binary_sens` | 0.8394919168591224 | — | — |
| `binary_spec` | 0.9475549900018179 | — | — |
| `sens_at_spec_0.95` | 0.836027713625866 | — | — |
| `sens_at_spec_0.97` | 0.7863741339491916 | — | — |
| `sens_at_spec_0.985` | 0.6293302540415704 | — | — |
| `auprc` | 0.6959023331605724 | — | — |
| `ternary_firm` | 0.9486012807549713 | — | — |

[metrics.csv](../results/models/lr_l1_tuned/e21_tuned_cv/metrics.csv) · [best_params_per_fold.csv](../results/models/lr_l1_tuned/e21_tuned_cv/best_params_per_fold.csv)

## lr_l2_tuned/e21_km_external

Descriptive · [Source JSON](../results/models/lr_l2_tuned/e21_km_external/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `auc` | 0.8680694470654842 | — | — |
| `binary_sens` | 0.8571428571428571 | — | — |
| `binary_spec` | 0.7727873183619551 | — | — |
| `sens_at_spec_0.95` | 0.6571428571428571 | — | — |
| `n_isolates` | 792 | — | — |
| `n_resistant` | 35 | — | — |
| `zero_training_overlap` | True | — | — |

[metrics.csv](../results/models/lr_l2_tuned/e21_km_external/metrics.csv) · [best_params.csv](../results/models/lr_l2_tuned/e21_km_external/best_params.csv)

## lr_l2_tuned/e21_tuned_cv

Descriptive · [Source JSON](../results/models/lr_l2_tuned/e21_tuned_cv/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `auc` | 0.9058801716924867 | — | — |
| `binary_sens` | 0.8383371824480369 | — | — |
| `binary_spec` | 0.9503726595164516 | — | — |
| `sens_at_spec_0.95` | 0.8383371824480369 | — | — |
| `sens_at_spec_0.97` | 0.7886836027713626 | — | — |
| `sens_at_spec_0.985` | 0.6304849884526559 | — | — |
| `auprc` | 0.6971430019924434 | — | — |
| `ternary_firm` | 0.9512133468149646 | — | — |

[metrics.csv](../results/models/lr_l2_tuned/e21_tuned_cv/metrics.csv) · [best_params_per_fold.csv](../results/models/lr_l2_tuned/e21_tuned_cv/best_params_per_fold.csv)

## lr_mechanism/e16_lof_ablation

Descriptive · [Source JSON](../results/models/lr_mechanism/e16_lof_ablation/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `lr_f2_auc` | 0.9085389366535497 | — | — |
| `lr_f2mech_auc` | 0.9087876842043836 | — | — |
| `lr_delta_auc` | 0.00024874755083392763 | — | — |
| `lr_f2_sens_at_spec_0.985` | 0.6039260969976905 | — | — |
| `lr_f2mech_sens_at_spec_0.985` | 0.6166281755196305 | — | — |
| `cat_s1_sens_firm_off` | 0.7924242424242425 | — | — |
| `cat_s1_sens_firm_on` | 0.7917933130699089 | — | — |
| `cat_s1_spec_firm_off` | 0.983502417749123 | — | — |
| `cat_s1_spec_firm_on` | 0.9871526456033498 | — | — |
| `cat_s1_n_U_off` | 661 | — | — |
| `cat_s1_n_U_on` | 702 | — | — |
| `s1_n_flips` | 41 | — | — |
| `s1_n_flips_resistant` | 2 | — | — |
| `s1_n_flips_susceptible` | 39 | — | — |
| `cs_CRyPTIC_to_NICD_sens_firm_off` | 0.6078431372549019 | — | — |
| `cs_CRyPTIC_to_NICD_sens_firm_on` | 0.6078431372549019 | — | — |
| `cs_CRyPTIC_to_NICD_spec_firm_off` | 0.991869918699187 | — | — |
| `cs_CRyPTIC_to_NICD_spec_firm_on` | 0.991869918699187 | — | — |
| `cs_CRyPTIC_to_NICD_n_flips` | 0 | — | — |
| `cs_CRyPTIC_to_NICD_n_flips_resistant` | 0 | — | — |
| `cs_NICD_to_CRyPTIC_sens_firm_off` | 0.5517241379310345 | — | — |
| `cs_NICD_to_CRyPTIC_sens_firm_on` | 0.5517241379310345 | — | — |
| `cs_NICD_to_CRyPTIC_spec_firm_off` | 0.9780872613946241 | — | — |
| `cs_NICD_to_CRyPTIC_spec_firm_on` | 0.9818164043406002 | — | — |
| `cs_NICD_to_CRyPTIC_n_flips` | 39 | — | — |
| `cs_NICD_to_CRyPTIC_n_flips_resistant` | 0 | — | — |

[lr_compare.csv](../results/models/lr_mechanism/e16_lof_ablation/lr_compare.csv) · [s1_veto_metrics.csv](../results/models/lr_mechanism/e16_lof_ablation/s1_veto_metrics.csv) · [cross_source_veto.csv](../results/models/lr_mechanism/e16_lof_ablation/cross_source_veto.csv)

## lr_mechanism/e16_lof_descriptive

Descriptive · [Source JSON](../results/models/lr_mechanism/e16_lof_descriptive/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `snobre_frac_susceptible_frameshift` | 0.250386398763524 | — | — |
| `snobre_reported` | 0.248 | — | — |
| `n_frameshift_carriers` | 647 | — | — |
| `n_frameshift_susceptible` | 162 | — | — |
| `n_truncating_variants` | 142 | — | — |
| `n_veto_fire_samples` | 43 | — | — |
| `n_variants_DBD` | 65 | — | — |
| `n_variants_mid` | 32 | — | — |
| `n_variants_late` | 43 | — | — |
| `n_variants_full_del` | 2 | — | — |

[lof_effects.csv](../results/models/lr_mechanism/e16_lof_descriptive/lof_effects.csv) · [region_summary.csv](../results/models/lr_mechanism/e16_lof_descriptive/region_summary.csv)

## lr_model/e04_s1_patient_cv

Passed · [Source JSON](../results/models/lr_model/e04_s1_patient_cv/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `auc` | 0.9085389366535497 | 0.909 | 0.01 |
| `auprc` | 0.6836352010125266 | 0.684 | 0.02 |
| `binary_sens` | 0.8348729792147807 | 0.835 | 0.02 |
| `binary_spec` | 0.9481003453917469 | 0.948 | 0.01 |
| `ternary_firm` | 0.9469160768452983 | 0.947 | 0.01 |
| `ternary_sens_among_firm` | 0.7346938775510204 | — | — |
| `ternary_spec_among_firm` | 0.9889201877934273 | — | — |
| `sens_at_spec_0.95` | 0.8302540415704388 | 0.83 | 0.03 |
| `sens_at_spec_0.97` | 0.7817551963048499 | 0.782 | 0.03 |
| `sens_at_spec_0.98` | 0.6524249422632794 | 0.652 | 0.03 |
| `sens_at_spec_0.985` | 0.6039260969976905 | 0.604 | 0.03 |
| `n_features_all_data` | 120 | — | — |

[metrics.csv](../results/models/lr_model/e04_s1_patient_cv/metrics.csv) · [thresholds_per_fold.csv](../results/models/lr_model/e04_s1_patient_cv/thresholds_per_fold.csv)

## lr_model/e05_ablation

Passed · [Source JSON](../results/models/lr_model/e05_ablation/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `F0_frs0.1_auc` | 0.7910392525734351 | 0.791 | 0.015 |
| `F1_frs0.1_auc` | 0.9045152613444627 | 0.905 | 0.015 |
| `F2_frs0.1_auc` | 0.9085389366535497 | 0.909 | 0.015 |
| `F2_frs0.75_auc` | 0.7986128283205279 | 0.799 | 0.015 |
| `F0_frs0.1_auprc` | 0.5025603185843016 | 0.503 | 0.015 |
| `F1_frs0.1_auprc` | 0.672147644479894 | 0.672 | 0.015 |
| `F2_frs0.1_auprc` | 0.6836352010125266 | 0.684 | 0.015 |
| `F2_frs0.75_auprc` | 0.5478810303005124 | 0.548 | 0.015 |

[ablation.csv](../results/models/lr_model/e05_ablation/ablation.csv)

## lr_model/e06_cross_source

Passed · [Source JSON](../results/models/lr_model/e06_cross_source/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `CRyPTIC_to_NICD_auc` | 0.9104980172497594 | 0.91 | 0.03 |
| `CRyPTIC_to_NICD_sens` | 0.8441890166028098 | 0.844 | 0.03 |
| `CRyPTIC_to_NICD_spec` | 0.963254593175853 | 0.963 | 0.03 |
| `NICD_to_CRyPTIC_auc` | 0.8277106164985712 | 0.83 | 0.03 |
| `NICD_to_CRyPTIC_sens` | 0.6746987951807228 | 0.675 | 0.03 |
| `NICD_to_CRyPTIC_spec` | 0.953205912814236 | 0.953 | 0.03 |

[cross_source.csv](../results/models/lr_model/e06_cross_source/cross_source.csv)

## lr_model/e07_loso

Passed · [Source JSON](../results/models/lr_model/e07_loso/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `pooled_sens` | 0.7875288683602771 | 0.788 | 0.02 |
| `pooled_spec` | 0.9537356844210143 | 0.954 | 0.01 |
| `per_site_spec_min_n300` | 0.9440993788819876 | — | — |
| `per_site_spec_max_n300` | 0.9785217940619078 | — | — |

[per_site.csv](../results/models/lr_model/e07_loso/per_site.csv)

## lr_model/e11_imbalance_learning_curve

Passed · [Source JSON](../results/models/lr_model/e11_imbalance_learning_curve/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `lc_auc_0.1` | 0.8963497818788353 | 0.898 | 0.01 |
| `lc_auc_0.25` | 0.9049963831896195 | 0.905 | 0.01 |
| `lc_auc_0.5` | 0.9050051995585098 | 0.908 | 0.01 |
| `lc_auc_0.75` | 0.9074998646057635 | 0.906 | 0.01 |
| `lc_auc_1.0` | 0.9085389366535497 | 0.909 | 0.01 |
| `balanced_auc` | 0.9050684884923296 | 0.904 | 0.01 |
| `balanced_auprc` | 0.6586597345941371 | 0.653 | 0.02 |
| `balanced_sens_at_spec98` | 0.6062355658198614 | 0.592 | 0.03 |

[learning_curve.csv](../results/models/lr_model/e11_imbalance_learning_curve/learning_curve.csv)

## lr_shrinkage/e14_shrinkage

Descriptive · [Source JSON](../results/models/lr_shrinkage/e14_shrinkage/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `n_variants` | 1042 | — | — |
| `mean_abs_shrinkage` | 0.14319372374084754 | — | — |
| `max_abs_shrinkage` | 0.9294249132357397 | — | — |
| `n_grade_flips_0.5` | 48 | — | — |
| `n_mech_classes` | 9 | — | — |

[variant_shrinkage.csv](../results/models/lr_shrinkage/e14_shrinkage/variant_shrinkage.csv) · [class_priors.csv](../results/models/lr_shrinkage/e14_shrinkage/class_priors.csv)

## rf_model/e05b_model_cv

Descriptive · [Source JSON](../results/models/rf_model/e05b_model_cv/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `auc` | 0.9093575470006923 | — | — |
| `auprc` | 0.6818038144317877 | — | — |
| `binary_sens` | 0.8337182448036952 | — | — |
| `binary_spec` | 0.9548263952008725 | — | — |
| `ternary_firm` | 0.964273677114931 | — | — |
| `sens_at_spec_0.95` | 0.8337182448036952 | — | — |
| `sens_at_spec_0.97` | 0.8048498845265589 | — | — |
| `sens_at_spec_0.98` | 0.7424942263279446 | — | — |
| `sens_at_spec_0.985` | 0.6593533487297921 | — | — |

[metrics.csv](../results/models/rf_model/e05b_model_cv/metrics.csv) · [thresholds_per_fold.csv](../results/models/rf_model/e05b_model_cv/thresholds_per_fold.csv)

## rf_tuned/e21_km_external

Descriptive · [Source JSON](../results/models/rf_tuned/e21_km_external/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `auc` | 0.8674655595395359 | — | — |
| `binary_sens` | 0.8285714285714286 | — | — |
| `binary_spec` | 0.7859973579920739 | — | — |
| `sens_at_spec_0.95` | 0.7714285714285715 | — | — |
| `n_isolates` | 792 | — | — |
| `n_resistant` | 35 | — | — |
| `zero_training_overlap` | True | — | — |

[metrics.csv](../results/models/rf_tuned/e21_km_external/metrics.csv) · [best_params.csv](../results/models/rf_tuned/e21_km_external/best_params.csv)

## rf_tuned/e21_tuned_cv

Descriptive · [Source JSON](../results/models/rf_tuned/e21_tuned_cv/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `auc` | 0.9103450327947932 | — | — |
| `binary_sens` | 0.8314087759815243 | — | — |
| `binary_spec` | 0.958007634975459 | — | — |
| `sens_at_spec_0.95` | 0.8406466512702079 | — | — |
| `sens_at_spec_0.97` | 0.815242494226328 | — | — |
| `sens_at_spec_0.985` | 0.6939953810623557 | — | — |
| `auprc` | 0.7010731255343863 | — | — |
| `ternary_firm` | 0.9736265588136165 | — | — |

[metrics.csv](../results/models/rf_tuned/e21_tuned_cv/metrics.csv) · [best_params_per_fold.csv](../results/models/rf_tuned/e21_tuned_cv/best_params_per_fold.csv)

## structure_test/e08_novel_variant

Passed · [Source JSON](../results/models/structure_test/e08_novel_variant/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `n_samples` | 312 | — | — |
| `n_distinct_sets` | 143 | — | — |
| `n_resistant` | 169 | — | — |
| `auc_region_count` | 0.5451442049075186 | — | — |
| `auc_plus_structure` | 0.5932469896966939 | 0.57 | 0.08 |

## svm_model/e05b_model_cv

Descriptive · [Source JSON](../results/models/svm_model/e05b_model_cv/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `auc` | 0.8953331181019786 | — | — |
| `auprc` | 0.6828598661710699 | — | — |
| `binary_sens` | 0.7979214780600462 | — | — |
| `binary_spec` | 0.9765497182330486 | — | — |
| `ternary_firm` | 0.07086282440175261 | — | — |
| `sens_at_spec_0.95` | 0.8013856812933026 | — | — |
| `sens_at_spec_0.97` | 0.8013856812933026 | — | — |
| `sens_at_spec_0.98` | 0.7551963048498845 | — | — |
| `sens_at_spec_0.985` | 0.6732101616628176 | — | — |

[metrics.csv](../results/models/svm_model/e05b_model_cv/metrics.csv) · [thresholds_per_fold.csv](../results/models/svm_model/e05b_model_cv/thresholds_per_fold.csv)

## svm_tuned/e21_km_external

Descriptive · [Source JSON](../results/models/svm_tuned/e21_km_external/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `auc` | 0.871994715984148 | — | — |
| `binary_sens` | 0.8 | — | — |
| `binary_spec` | 0.9326287978863936 | — | — |
| `sens_at_spec_0.95` | 0.7714285714285715 | — | — |
| `n_isolates` | 792 | — | — |
| `n_resistant` | 35 | — | — |
| `zero_training_overlap` | True | — | — |

[metrics.csv](../results/models/svm_tuned/e21_km_external/metrics.csv) · [best_params.csv](../results/models/svm_tuned/e21_km_external/best_params.csv)

## svm_tuned/e21_tuned_cv

Descriptive · [Source JSON](../results/models/svm_tuned/e21_tuned_cv/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `auc` | 0.8984782002684374 | — | — |
| `binary_sens` | 0.815242494226328 | — | — |
| `binary_spec` | 0.9670059989092892 | — | — |
| `sens_at_spec_0.95` | 0.8210161662817552 | — | — |
| `sens_at_spec_0.97` | 0.8013856812933026 | — | — |
| `sens_at_spec_0.985` | 0.6859122401847575 | — | — |
| `auprc` | 0.6652478869507001 | — | — |
| `ternary_firm` | 0.7944893832153691 | — | — |

[metrics.csv](../results/models/svm_tuned/e21_tuned_cv/metrics.csv) · [best_params_per_fold.csv](../results/models/svm_tuned/e21_tuned_cv/best_params_per_fold.csv)

## catalogue/e03_s1_catalogue_cv

Passed · [Source JSON](../results/approaches/catalogue/e03_s1_catalogue_cv/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `firm` | 0.944304010785305 | 0.944 | 0.01 |
| `sens_firm` | 0.7924242424242425 | 0.792 | 0.015 |
| `spec_firm` | 0.983502417749123 | 0.984 | 0.004 |
| `sens_U_as_S` | 0.6039260969976905 | 0.604 | 0.015 |
| `spec_U_as_S` | 0.984184693692056 | 0.984 | 0.004 |
| `sens_U_as_R` | 0.8418013856812933 | 0.842 | 0.015 |
| `spec_U_as_R` | 0.9428285766224322 | 0.943 | 0.01 |
| `n_U` | 661 | 661 | 15 |
| `n_U_resistant` | 206 | 206 | 15 |

[metrics.csv](../results/approaches/catalogue/e03_s1_catalogue_cv/metrics.csv)

## conformal/e14_certificate_audit

Descriptive · [Source JSON](../results/approaches/conformal/e14_certificate_audit/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `alpha` | 0.015 | — | — |
| `beta` | 0.2 | — | — |
| `delta` | 0.05 | — | — |
| `pi_train` | 0.07296932928884395 | — | — |
| `baseline_emp_s_error` | 0.014595808383233533 | — | — |
| `baseline_emp_r_error` | 0.21454545454545454 | — | — |
| `baseline_n_S` | 10688 | — | — |
| `baseline_n_R` | 550 | — | — |
| `baseline_n_U` | 630 | — | — |
| `baseline_n_U_resistant` | 278 | — | — |
| `certified_t_R` | 0.8766 | — | — |
| `certified_bound_s` | 1.0 | — | — |
| `certified_bound_r` | 0.19630691917422763 | — | — |
| `certified_s_certified` | False | — | — |
| `certified_r_certified` | True | — | — |
| `certified_emp_r_error` | 0.14342629482071714 | — | — |
| `certified_n_S` | 0 | — | — |
| `certified_n_R` | 251 | — | — |
| `certified_n_U` | 11617 | — | — |
| `certified_n_U_resistant` | 651 | — | — |
| `fold_n_tS_certified` | 0 | — | — |
| `fold_emp_r_min` | 0.16363636363636364 | — | — |
| `fold_emp_r_max` | 0.21311475409836064 | — | — |
| `s_best_achievable_bound` | 0.01551592944072069 | — | — |
| `s_best_achievable_t` | 0.022 | — | — |
| `out_of_site_n_sites` | 11 | — | — |
| `out_of_site_violations_s` | 0 | — | — |
| `out_of_site_violations_r` | 5 | — | — |
| `cs_CRyPTIC_to_NICD_bound_s` | 0.00972879828658099 | — | — |
| `cs_CRyPTIC_to_NICD_bound_r` | 1.0 | — | — |
| `cs_CRyPTIC_to_NICD_emp_s_error` | 0.671280276816609 | — | — |
| `cs_CRyPTIC_to_NICD_violates_s` | True | — | — |
| `cs_CRyPTIC_to_NICD_violates_r` | False | — | — |
| `cs_NICD_to_CRyPTIC_bound_s` | 1.0 | — | — |
| `cs_NICD_to_CRyPTIC_bound_r` | 0.17913304251203174 | — | — |
| `cs_NICD_to_CRyPTIC_emp_r_error` | 0.9772288691624855 | — | — |
| `cs_NICD_to_CRyPTIC_violates_s` | False | — | — |
| `cs_NICD_to_CRyPTIC_violates_r` | True | — | — |
| `weighted_emp_r_error_pi002` | 0.39239733124688225 | — | — |
| `weighted_emp_r_error_pi005` | 0.20026682227918696 | — | — |

[policy_compare.csv](../results/approaches/conformal/e14_certificate_audit/policy_compare.csv) · [per_fold.csv](../results/approaches/conformal/e14_certificate_audit/per_fold.csv) · [out_of_site.csv](../results/approaches/conformal/e14_certificate_audit/out_of_site.csv) · [cross_source.csv](../results/approaches/conformal/e14_certificate_audit/cross_source.csv)

## data_prep/e01_prep

Passed · [Source JSON](../results/approaches/data_prep/e01_prep/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `n_samples` | 11868 | 11868 | 0 |
| `n_resistant` | 866 | 866 | 0 |
| `cryptic_n` | 10704 | 10704 | 0 |
| `cryptic_R` | 83 | 83 | 0 |
| `nicd_n` | 1164 | 1164 | 0 |
| `nicd_R` | 783 | 783 | 0 |
| `unique_subjid` | 11611 | 11611 | 0 |
| `mutation_rows` | 25554 | 25554 | 0 |
| `minor_allele_rows` | 1115 | 1115 | 0 |
| `median_FRSf_minor` | 0.375 | — | — |
| `null_rows` | 11 | 11 | 0 |
| `null_rows_samples` | 4 | 4 | 0 |

[counts.csv](../results/approaches/data_prep/e01_prep/counts.csv)

## external_choi/e_choi_inventory

Descriptive · [Source JSON](../results/approaches/external_choi/e_choi_inventory/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `n_accessions` | 49 | — | — |
| `country` | South Korea | — | — |
| `breakpoint_caveat` | paper breakpoint >=0.5 != WHO CC; reclassify at 1.0 / 0.25 first | — | — |
| `genotype_pipeline_blocked` | True | — | — |

## external_km/e_km_cohort

Descriptive · [Source JSON](../results/approaches/external_km/e_km_cohort/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `n_isolates` | 792 | — | — |
| `n_resistant` | 35 | — | — |
| `n_susceptible` | 757 | — | — |
| `prevalence` | 0.04419191919191919 | — | — |
| `quality` | MEDIUM | — | — |
| `assay_cc_undocumented` | True | — | — |
| `zero_training_overlap` | True | — | — |
| `n_overlap` | 0 | — | — |
| `genotype_pipeline_blocked` | True | — | — |
| `n_lineage2` | 559 | — | — |
| `n_lineage4` | 198 | — | — |
| `n_mixed` | 20 | — | — |
| `n_lineage1` | 9 | — | — |
| `n_lineage3` | 6 | — | — |
| `nR_lineage2` | 22 | — | — |
| `nR_lineage4` | 9 | — | — |
| `nR_lineage1` | 4 | — | — |

[lineage_stratification.csv](../results/approaches/external_km/e_km_cohort/lineage_stratification.csv)

## external_km/e_km_validation

Descriptive · [Source JSON](../results/approaches/external_km/e_km_validation/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `n_isolates` | 792 | — | — |
| `n_resistant` | 35 | — | — |
| `zero_training_overlap` | True | — | — |
| `n_with_variants` | 790 | — | — |
| `cat_firm` | 0.7752525252525253 | — | — |
| `cat_sens_firm` | 0.7916666666666666 | — | — |
| `cat_spec_firm` | 0.988135593220339 | — | — |
| `cat_sens_U_as_S` | 0.5428571428571428 | — | — |
| `cat_spec_U_as_S` | 0.9907529722589168 | — | — |
| `cat_sens_U_as_R` | 0.8571428571428571 | — | — |
| `cat_spec_U_as_R` | 0.7701453104359313 | — | — |
| `cat_n_U` | 178 | — | — |
| `cat_n_U_resistant` | 11 | — | — |
| `model_auc` | 0.8681826759765995 | — | — |
| `model_binary_sens` | 0.8571428571428571 | — | — |
| `model_binary_spec` | 0.7727873183619551 | — | — |
| `model_sens_at_spec_0.95` | 0.6857142857142857 | — | — |
| `model_sens_at_spec_0.97` | 0.45714285714285713 | — | — |
| `lin_lineage1_sens_firm` | 1.0 | — | — |
| `lin_lineage1_spec_firm` | 1.0 | — | — |
| `lin_lineage2_sens_firm` | 0.6428571428571429 | — | — |
| `lin_lineage2_spec_firm` | 0.9896907216494846 | — | — |
| `lin_lineage4_sens_firm` | 1.0 | — | — |
| `lin_lineage4_spec_firm` | 0.9887005649717514 | — | — |

[metrics.csv](../results/approaches/external_km/e_km_validation/metrics.csv) · [lineage_stratification.csv](../results/approaches/external_km/e_km_validation/lineage_stratification.csv)

## external_sobol/e_sobol_inventory

Descriptive · [Source JSON](../results/approaches/external_sobol/e_sobol_inventory/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `n_accessions` | 744 | — | — |
| `country` | South Africa | — | — |
| `claimed_bdst_pdst_n` | 423 | — | — |
| `claimed_pdst_n_verified` | False | — | — |
| `assay` | MGIT 1.0 (per paper; assay-matched to NICD) | — | — |
| `paper_89_is_genotypic` | True | — | — |
| `genotype_pipeline_blocked` | True | — | — |

## head_to_head/e05b_comparison

Descriptive · [Source JSON](../results/approaches/head_to_head/e05b_comparison/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `lr_model_auc` | 0.9085389366535497 | — | — |
| `lr_model_auprc` | 0.6836352010125266 | — | — |
| `lr_model_binary_sens` | 0.8348729792147807 | — | — |
| `lr_model_binary_spec` | 0.9481003453917469 | — | — |
| `lr_model_sens_at_spec_0.97` | 0.7817551963048499 | — | — |
| `lr_model_sens_at_spec_0.985` | 0.6039260969976905 | — | — |
| `lr_l1_auc` | 0.9053351836512614 | — | — |
| `lr_l1_auprc` | 0.6752893399877362 | — | — |
| `lr_l1_binary_sens` | 0.8394919168591224 | — | — |
| `lr_l1_binary_spec` | 0.9454644610070896 | — | — |
| `lr_l1_sens_at_spec_0.97` | 0.789838337182448 | — | — |
| `lr_l1_sens_at_spec_0.985` | 0.5993071593533488 | — | — |
| `rf_model_auc` | 0.9093575470006923 | — | — |
| `rf_model_auprc` | 0.6818038144317877 | — | — |
| `rf_model_binary_sens` | 0.8337182448036952 | — | — |
| `rf_model_binary_spec` | 0.9548263952008725 | — | — |
| `rf_model_sens_at_spec_0.97` | 0.8048498845265589 | — | — |
| `rf_model_sens_at_spec_0.985` | 0.6593533487297921 | — | — |
| `hgb_model_auc` | 0.9052400403369869 | — | — |
| `hgb_model_auprc` | 0.678652962897792 | — | — |
| `hgb_model_binary_sens` | 0.8371824480369515 | — | — |
| `hgb_model_binary_spec` | 0.9528267587711325 | — | — |
| `hgb_model_sens_at_spec_0.97` | 0.8140877598152425 | — | — |
| `hgb_model_sens_at_spec_0.985` | 0.625866050808314 | — | — |
| `svm_model_auc` | 0.8953331181019786 | — | — |
| `svm_model_auprc` | 0.6828598661710699 | — | — |
| `svm_model_binary_sens` | 0.7979214780600462 | — | — |
| `svm_model_binary_spec` | 0.9765497182330486 | — | — |
| `svm_model_sens_at_spec_0.97` | 0.8013856812933026 | — | — |
| `svm_model_sens_at_spec_0.985` | 0.6732101616628176 | — | — |

[comparison.csv](../results/approaches/head_to_head/e05b_comparison/comparison.csv)

## head_to_head/e21_comparison

Descriptive · [Source JSON](../results/approaches/head_to_head/e21_comparison/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `S1_leader` | rf_tuned | — | — |

[comparison.csv](../results/approaches/head_to_head/e21_comparison/comparison.csv)

## head_to_head/g7_adlard

Descriptive · [Source JSON](../results/approaches/head_to_head/g7_adlard/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `shuffled_firm` | 0.9448938321536906 | — | — |
| `shuffled_sens_firm` | 0.7955223880597015 | — | — |
| `shuffled_spec_firm` | 0.9837822458270106 | — | — |
| `shuffled_n_U` | 654 | — | — |
| `grouped_firm` | 0.944304010785305 | — | — |
| `grouped_sens_firm` | 0.7924242424242425 | — | — |
| `grouped_spec_firm` | 0.983502417749123 | — | — |
| `grouped_n_U` | 661 | — | — |
| `sens_firm_shuffled_minus_grouped` | 0.003098145635459071 | — | — |
| `ghosh_reconstructable` | False | — | — |

[shuffled_vs_grouped.csv](../results/approaches/head_to_head/g7_adlard/shuffled_vs_grouped.csv)

## head_to_head/g7_who_v2

Descriptive · [Source JSON](../results/approaches/head_to_head/g7_who_v2/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `who_frs0.1_firm` | 0.9607347489046174 | — | — |
| `who_frs0.1_sens_firm` | 0.8010899182561307 | — | — |
| `who_frs0.1_spec_firm` | 0.9811586051743532 | — | — |
| `who_frs0.1_sens_U_as_S` | 0.6789838337182448 | — | — |
| `who_frs0.1_sens_U_as_R` | 0.8314087759815243 | — | — |
| `who_frs0.25_firm` | 0.961071789686552 | — | — |
| `who_frs0.25_sens_firm` | 0.7679558011049724 | — | — |
| `who_frs0.25_spec_firm` | 0.9827747612806591 | — | — |
| `who_frs0.25_sens_U_as_S` | 0.6420323325635104 | — | — |
| `who_frs0.25_sens_U_as_R` | 0.8060046189376443 | — | — |
| `who_frs0.75_firm` | 0.9679811257162116 | — | — |
| `who_frs0.75_sens_firm` | 0.560625814863103 | — | — |
| `who_frs0.75_spec_firm` | 0.9868482417684917 | — | — |
| `who_frs0.75_sens_U_as_S` | 0.49653579676674364 | — | — |
| `who_frs0.75_sens_U_as_R` | 0.6108545034642032 | — | — |
| `our_catalogue_grouped_sens_firm` | 0.7924242424242425 | — | — |
| `our_catalogue_grouped_spec_firm` | 0.983502417749123 | — | — |
| `our_model_ternary_sens_among_firm` | 0.7346938775510204 | — | — |
| `our_model_ternary_spec_among_firm` | 0.9889201877934273 | — | — |

[who_vs_ours.csv](../results/approaches/head_to_head/g7_who_v2/who_vs_ours.csv)

## hybrid/e09_hybrid

Passed · [Source JSON](../results/approaches/hybrid/e09_hybrid/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `n_U` | 661 | — | — |
| `n_U_resistant` | 206 | — | — |
| `auc_ml_within_U` | 0.7442227675237384 | 0.744 | 0.03 |
| `resolved_S` | 173 | 173 | 15 |
| `resolved_S_resistant` | 17 | 17 | 15 |
| `resolved_R` | 126 | 126 | 15 |
| `resolved_R_susceptible` | 51 | 51 | 15 |
| `remain_U` | 362 | 362 | 15 |
| `remain_U_resistant` | 114 | 114 | 15 |

## mic_audit/e17_borderline_mic

Descriptive · [Source JSON](../results/approaches/mic_audit/e17_borderline_mic/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `n_borderline` | 454 | — | — |
| `n_borderline_resistant` | 191 | — | — |
| `borderline_share_of_positives` | 0.22055427251732102 | — | — |
| `n_cry_borderline` | 236 | — | — |
| `n_nicd_borderline` | 218 | — | — |
| `n_nicd_mic_available` | 612 | — | — |
| `full_auc` | 0.9085389366535497 | — | — |
| `full_binary_sens` | 0.8348729792147807 | — | — |
| `full_binary_spec` | 0.9481003453917469 | — | — |
| `full_sens_at_spec_0.985` | 0.6039260969976905 | — | — |
| `full_n` | 11868 | — | — |
| `full_n_resistant` | 866 | — | — |
| `excl_borderline_auc` | 0.9044757598645298 | — | — |
| `excl_borderline_binary_sens` | 0.8237037037037037 | — | — |
| `excl_borderline_binary_spec` | 0.9555824564670826 | — | — |
| `excl_borderline_sens_at_spec_0.985` | 0.6577777777777778 | — | — |
| `excl_borderline_n` | 11414 | — | — |
| `excl_borderline_n_resistant` | 675 | — | — |
| `delta_auc_excl_minus_full` | -0.004063176789019907 | — | — |

[borderline_by_mic.csv](../results/approaches/mic_audit/e17_borderline_mic/borderline_by_mic.csv) · [full_vs_excluded.csv](../results/approaches/mic_audit/e17_borderline_mic/full_vs_excluded.csv)

## predictive_values/e13_ppv_npv

Passed · [Source JSON](../results/approaches/predictive_values/e13_ppv_npv/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `model_ppv_0.01` | 0.1435986159169549 | 0.14 | 0.02 |
| `model_npv_0.01` | 0.9981957121630227 | 0.998 | 0.003 |
| `cat_ppv_0.01` | 0.33277169334456597 | 0.33 | 0.02 |
| `cat_npv_0.01` | 0.9978489336856985 | 0.998 | 0.003 |
| `model_ppv_0.02` | 0.2530487804878047 | 0.25 | 0.02 |
| `model_npv_0.02` | 0.996361301369863 | 0.996 | 0.003 |
| `cat_ppv_0.02` | 0.5019059720457432 | 0.5 | 0.02 |
| `cat_npv_0.02` | 0.9956634865568084 | 0.996 | 0.003 |
| `model_ppv_0.05` | 0.4662921348314604 | 0.47 | 0.02 |
| `model_npv_0.05` | 0.990669593852909 | 0.991 | 0.003 |
| `cat_ppv_0.05` | 0.7221206581352833 | 0.72 | 0.02 |
| `cat_npv_0.05` | 0.9888924151063155 | 0.989 | 0.003 |
| `model_ppv_0.1` | 0.6484374999999998 | 0.65 | 0.02 |
| `model_npv_0.1` | 0.9805045871559633 | 0.981 | 0.003 |
| `cat_ppv_0.1` | 0.8458244111349036 | 0.85 | 0.02 |
| `cat_npv_0.1` | 0.9768365320979484 | 0.977 | 0.003 |

[ppv_npv.csv](../results/approaches/predictive_values/e13_ppv_npv/ppv_npv.csv)

## prevalence_correction/e10_prevalence_shift

Passed · [Source JSON](../results/approaches/prevalence_correction/e10_prevalence_shift/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `CRyPTIC_to_NICD_ece_raw` | 0.5684131836680266 | 0.568 | 0.03 |
| `CRyPTIC_to_NICD_ece_adj` | 0.10120066429985355 | 0.101 | 0.03 |
| `CRyPTIC_to_NICD_brier_raw` | 0.5009009962128409 | 0.501 | 0.02 |
| `CRyPTIC_to_NICD_brier_adj` | 0.12376224807239711 | 0.124 | 0.02 |
| `NICD_to_CRyPTIC_ece_raw` | 0.28816773432562826 | 0.288 | 0.01 |
| `NICD_to_CRyPTIC_ece_adj` | 0.012789839603910978 | 0.013 | 0.01 |
| `NICD_to_CRyPTIC_brier_raw` | 0.10203671807064739 | 0.102 | 0.005 |
| `NICD_to_CRyPTIC_brier_adj` | 0.016917730721914482 | 0.017 | 0.005 |

[prevalence_shift.csv](../results/approaches/prevalence_correction/e10_prevalence_shift/prevalence_shift.csv)

## priority_list/e12_variant_priority

Passed · [Source JSON](../results/approaches/priority_list/e12_variant_priority/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `n_U` | 661 | — | — |
| `n_distinct_sets` | 332 | 332 | 15 |
| `top5_coverage` | 0.27685325264750377 | 0.28 | 0.03 |
| `top10_coverage` | 0.3434190620272315 | 0.34 | 0.03 |
| `top20_coverage` | 0.3978819969742814 | 0.4 | 0.03 |
| `top50_coverage` | 0.4962178517397882 | 0.5 | 0.03 |
| `largest_set` | Rv0678@c-11a | — | — |
| `largest_set_n` | 102 | — | — |
| `largest_set_pctR` | 0.08823529411764706 | — | — |

[top_sets.csv](../results/approaches/priority_list/e12_variant_priority/top_sets.csv) · [evidence_cards.csv](../results/approaches/priority_list/e12_variant_priority/evidence_cards.csv)

## research_archive/g8_evaluation_context

Descriptive · [Source JSON](../results/approaches/research_archive/g8_evaluation_context/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `interval_level` | 0.95 | — | — |
| `bootstrap_replicates` | 500 | — | — |
| `n_resistant` | 866 | — | — |
| `n_resistant_without_retained_core_variant` | 137 | — | — |
| `fraction_without_retained_core_variant` | 0.1581986143187067 | — | — |

[intervals.csv](../results/approaches/research_archive/g8_evaluation_context/intervals.csv)

## ship_decision/e19_km_veto

Descriptive · [Source JSON](../results/approaches/ship_decision/e19_km_veto/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `cat_sens_firm_off` | 0.7916666666666666 | — | — |
| `cat_spec_firm_off` | 0.988135593220339 | — | — |
| `cat_n_U_off` | 178 | — | — |
| `cat_sens_firm_on` | 0.7916666666666666 | — | — |
| `cat_spec_firm_on` | 0.988135593220339 | — | — |
| `cat_n_U_on` | 178 | — | — |
| `n_flips` | 0 | — | — |
| `n_flips_resistant` | 0 | — | — |
| `n_flips_susceptible` | 0 | — | — |
| `n_true_R_to_S` | 0 | — | — |
| `all_flips_R_to_U` | True | — | — |
| `n_isolates` | 792 | — | — |
| `n_resistant` | 35 | — | — |
| `zero_training_overlap` | True | — | — |
| `model_auc` | 0.8681826759765995 | — | — |

[lineage_stratification.csv](../results/approaches/ship_decision/e19_km_veto/lineage_stratification.csv)

## ship_decision/e19_ladder_summary

Descriptive · [Source JSON](../results/approaches/ship_decision/e19_ladder_summary/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `model` | models/lr_model | — | — |
| `epistatic_veto` | True | — | — |

[ladder.csv](../results/approaches/ship_decision/e19_ladder_summary/ladder.csv)

## voi_queue/e15_label_reveal_sim

Descriptive · [Source JSON](../results/approaches/voi_queue/e15_label_reveal_sim/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `info_gain_pct_K1` | 0.2405446293494705 | — | — |
| `info_gain_pct_K2` | 0.2405446293494705 | — | — |
| `info_gain_pct_K5` | 0.3146747352496218 | — | — |
| `info_gain_pct_K10` | 0.3812405446293495 | — | — |
| `info_gain_pct_K20` | 0.4372163388804841 | — | — |
| `info_gain_pct_K50` | 0.5370650529500757 | — | — |
| `count_pct_K1` | 0.19213313161875945 | — | — |
| `count_pct_K2` | 0.2405446293494705 | — | — |
| `count_pct_K5` | 0.3146747352496218 | — | — |
| `count_pct_K10` | 0.3812405446293495 | — | — |
| `count_pct_K20` | 0.4372163388804841 | — | — |
| `count_pct_K50` | 0.5264750378214826 | — | — |
| `random_pct_K1` | 0.19364599092284418 | — | — |
| `random_pct_K2` | 0.21482602118003025 | — | — |
| `random_pct_K5` | 0.22087745839636913 | — | — |
| `random_pct_K10` | 0.2314674735249622 | — | — |
| `random_pct_K20` | 0.24659606656580937 | — | — |
| `random_pct_K50` | 0.29954614220877457 | — | — |
| `total_U` | 661 | — | — |
| `n_distinct_sets` | 332 | — | — |
| `info_gain_minus_count_K10` | 0.0 | — | — |
| `info_gain_minus_count_K50` | 0.010590015128593033 | — | — |

[resolution_curves.csv](../results/approaches/voi_queue/e15_label_reveal_sim/resolution_curves.csv) · [info_gain_ranking.csv](../results/approaches/voi_queue/e15_label_reveal_sim/info_gain_ranking.csv)

## who_baseline/e02_who_as_is_frs0.1

Passed · [Source JSON](../results/approaches/who_baseline/e02_who_as_is_frs0.1/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `firm` | 0.9607347489046174 | 0.961 | 0.005 |
| `sens_firm` | 0.8010899182561307 | 0.801 | 0.005 |
| `spec_firm` | 0.9811586051743532 | 0.981 | 0.005 |
| `sens_U_as_S` | 0.6789838337182448 | 0.679 | 0.005 |
| `spec_U_as_S` | 0.981730594437375 | — | — |
| `sens_U_as_R` | 0.8314087759815243 | — | — |
| `spec_U_as_R` | 0.9513724777313216 | — | — |
| `n_U` | 466 | — | — |
| `n_U_resistant` | 132 | — | — |

[metrics.csv](../results/approaches/who_baseline/e02_who_as_is_frs0.1/metrics.csv)

## who_baseline/e02_who_as_is_frs0.75

Passed · [Source JSON](../results/approaches/who_baseline/e02_who_as_is_frs0.75/results.json)

| Metric | Observed | Expected | Tolerance |
|---|---|---|---|
| `firm` | 0.9679811257162116 | 0.968 | 0.005 |
| `sens_firm` | 0.560625814863103 | 0.561 | 0.005 |
| `spec_firm` | 0.9868482417684917 | 0.987 | 0.005 |
| `sens_U_as_S` | 0.49653579676674364 | 0.497 | 0.005 |
| `spec_U_as_S` | 0.987184148336666 | — | — |
| `sens_U_as_R` | 0.6108545034642032 | — | — |
| `spec_U_as_R` | 0.9616433375749863 | — | — |
| `n_U` | 380 | — | — |
| `n_U_resistant` | 99 | — | — |

[metrics.csv](../results/approaches/who_baseline/e02_who_as_is_frs0.75/metrics.csv)
