# Production ship decision

> **Type:** analysis · **Approach:** `ship_decision` · **0/2 tests pass** ❌

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E19 · KM Veto | — | cat_sens_firm_off=0.7917 | [report](e19_km_veto/REPORT.md) |
| E19 · Ladder Summary | — |  | [report](e19_ladder_summary/REPORT.md) |

## E19 · KM Veto — DESCRIPTIVE •

Adopt the audited Tier A veto modifier.

| Metric | Value |
|---|---|
| `cat_sens_firm_off` | 0.7917 |
| `cat_spec_firm_off` | 0.9881 |
| `cat_n_U_off` | 178 |
| `cat_sens_firm_on` | 0.7917 |
| `cat_spec_firm_on` | 0.9881 |
| `cat_n_U_on` | 178 |
| `n_flips` | 0 |
| `n_flips_resistant` | 0 |
| `n_flips_susceptible` | 0 |
| `n_true_R_to_S` | 0 |
| `n_isolates` | 792 |
| `n_resistant` | 35 |
| `model_auc` | 0.8682 |

- **all_flips_R_to_U:** `True`
- **zero_training_overlap:** `True`
- **per_lineage:** `{'lineage1': {'n': 9, 'n_resistant': 4, 'cat_sens_firm_off': 1.0, 'cat_spec_firm_off': 1.0, 'cat_n_U_off': 1, 'cat_sens_firm_on': 1.0, 'cat_spec_firm_on': 1.0, 'cat_n_U_on': 1, 'n_flips': 0, 'n_flips_resistant': 0, 'n_flips_susceptible': 0, 'n_true_R_to_S': 0, 'all_flips_R_to_U': True}, 'lineage2': {'n': 559, 'n_resistant': 22, 'cat_sens_firm_off': 0.6428571428571429, 'cat_spec_firm_off': 0.9896907216494846, 'cat_n_U_off': 157, 'cat_sens_firm_on': 0.6428571428571429, 'cat_spec_firm_on': 0.9896907216494846, 'cat_n_U_on': 157, 'n_flips': 0, 'n_flips_resistant': 0, 'n_flips_susceptible': 0, 'n_true_R_to_S': 0, 'all_flips_R_to_U': True}, 'lineage3': {'n': 6, 'n_resistant': 0, 'cat_sens_firm_off': None, 'cat_spec_firm_off': 1.0, 'cat_n_U_off': 0, 'cat_sens_firm_on': None, 'cat_spec_firm_on': 1.0, 'cat_n_U_on': 0, 'n_flips': 0, 'n_flips_resistant': 0, 'n_flips_susceptible': 0, 'n_true_R_to_S': 0, 'all_flips_R_to_U': True}, 'lineage4': {'n': 198, 'n_resistant': 9, 'cat_sens_firm_off': 1.0, 'cat_spec_firm_off': 0.9887005649717514, 'cat_n_U_off': 14, 'cat_sens_firm_on': 1.0, 'cat_spec_firm_on': 0.9887005649717514, 'cat_n_U_on': 14, 'n_flips': 0, 'n_flips_resistant': 0, 'n_flips_susceptible': 0, 'n_true_R_to_S': 0, 'all_flips_R_to_U': True}, 'mixed': {'n': 20, 'n_resistant': 0, 'cat_sens_firm_off': None, 'cat_spec_firm_off': 0.9285714285714286, 'cat_n_U_off': 6, 'cat_sens_firm_on': None, 'cat_spec_firm_on': 0.9285714285714286, 'cat_n_U_on': 6, 'n_flips': 0, 'n_flips_resistant': 0, 'n_flips_susceptible': 0, 'n_true_R_to_S': 0, 'all_flips_R_to_U': True}}`

_Caveats: Research use only. KM sensitivity analysis with 35 resistant isolates, undocumented assay and critical concentration, and medium-quality phenotypes._

## E19 · Ladder Summary — DESCRIPTIVE •

Keep golden L2 LR (F2). Adopt the audited Tier A veto modifier.

- **model:** `models/lr_model`
- **epistatic_veto:** `True`

_Caveats: Research use only. KM sensitivity analysis with 35 resistant isolates, undocumented assay and critical concentration, and medium-quality phenotypes. Report the certificate bounds with their confidence level and lab-shift results._

