# External — KM validation

> **Type:** validation · **Approach:** `external_km` · **0/2 tests pass** ❌

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E KM Cohort | — | n_isolates=792 | [report](e_km_cohort/REPORT.md) |
| E KM Validation | — | n_isolates=792 | [report](e_km_validation/REPORT.md) |

## E KM Cohort — DESCRIPTIVE •

Verified external cohort: 792 isolates, 35 resistant (prevalence 0.0442), ZERO training overlap. All MEDIUM quality; METHOD_3/METHOD_CC undocumented. Lineage mix L2=559/L4=198. Per P-EXT this is a SENSITIVITY ANALYSIS, not headline proof. Genotype->phenotype sens/spec is BLOCKED: external variant calls are not local.

| Metric | Value |
|---|---|
| `n_isolates` | 792 |
| `n_resistant` | 35 |
| `n_susceptible` | 757 |
| `prevalence` | 0.0442 |
| `n_overlap` | 0 |
| `n_lineage2` | 559 |
| `n_lineage4` | 198 |
| `n_mixed` | 20 |
| `n_lineage1` | 9 |
| `n_lineage3` | 6 |
| `nR_lineage2` | 22 |
| `nR_lineage4` | 9 |
| `nR_lineage1` | 4 |

- **quality:** `MEDIUM`
- **assay_cc_undocumented:** `True`
- **zero_training_overlap:** `True`
- **genotype_pipeline_blocked:** `True`

_Caveats: Assay/CC UNDOCUMENTED (MEDIUM) -> cannot assume WHO critical concentration; stratify by lineage. NO external sens/spec computed (no external genotypes). NO Indian external resistant isolates exist -> no Indian-sensitivity claim ever. To finish: acquire CRyPTIC MUTATIONS (BDQ genes) for these 792 isolates, then run catalogue+model._

## E KM Validation — DESCRIPTIVE •

EXTERNAL test on 792 non-overlapping isolates: catalogue sens_firm 0.7917 / spec_firm 0.9881; model AUC 0.8682, binary sens 0.8571 / spec 0.7728. Stratified by lineage in tables.

| Metric | Value |
|---|---|
| `n_isolates` | 792 |
| `n_resistant` | 35 |
| `n_with_variants` | 790 |
| `cat_firm` | 0.7753 |
| `cat_sens_firm` | 0.7917 |
| `cat_spec_firm` | 0.9881 |
| `cat_sens_U_as_S` | 0.5429 |
| `cat_spec_U_as_S` | 0.9908 |
| `cat_sens_U_as_R` | 0.8571 |
| `cat_spec_U_as_R` | 0.7701 |
| `cat_n_U` | 178 |
| `cat_n_U_resistant` | 11 |
| `model_auc` | 0.8682 |
| `model_binary_sens` | 0.8571 |
| `model_binary_spec` | 0.7728 |
| `model_sens_at_spec_0.95` | 0.6857 |
| `model_sens_at_spec_0.97` | 0.4571 |
| `lin_lineage1_sens_firm` | 1.0000 |
| `lin_lineage1_spec_firm` | 1.0000 |
| `lin_lineage2_sens_firm` | 0.6429 |
| `lin_lineage2_spec_firm` | 0.9897 |
| `lin_lineage4_sens_firm` | 1.0000 |
| `lin_lineage4_spec_firm` | 0.9887 |

- **zero_training_overlap:** `True`

_Caveats: SENSITIVITY ANALYSIS, not headline proof: KM assay/CC UNDOCUMENTED (MEDIUM quality). Only 35 resistant -> wide sensitivity CI. No Indian external resistant isolates exist. Production model has NOT seen these isolates (zero overlap verified). match != beat._

