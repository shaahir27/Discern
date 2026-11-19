# E19 · KM Veto

> **Approach:** `ship_decision` · **Test:** `e19_km_veto` · **DESCRIPTIVE •**

Adopt the audited Tier A veto modifier.

## Metrics

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

## Method
Apply the frozen golden catalogue to the fixed KM cohort. Audit the existing veto on each call and use the saved production KM AUC as the LR anchor.

## Caveats / negative results
Research use only. KM sensitivity analysis with 35 resistant isolates, undocumented assay and critical concentration, and medium-quality phenotypes.

## Reproduce
```bash
C:/Users/shaah/miniconda3/envs/discern/python.exe eval/e19_ship_decision.py
```

## Provenance
- **Date (UTC):** 2026-11-19T06:00:12.448251+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `call_flips.csv`
- `lineage_stratification.csv`
