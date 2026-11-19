# E KM Validation

> **Approach:** `external_km` · **Test:** `e_km_validation` · **DESCRIPTIVE •**

EXTERNAL test on 792 non-overlapping isolates: catalogue sens_firm 0.7917 / spec_firm 0.9881; model AUC 0.8682, binary sens 0.8571 / spec 0.7728. Stratified by lineage in tables.

## Metrics

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

## Method
Production catalogue (built on all 11,868) + production calibrated LR applied to KMVALIDATION2025 DNA variants from CRyPTIC v3.4.0 MUTATIONS.parquet (Zenodo 15680920, CC-BY 4.0). 792 isolates, ZERO training overlap. Compared to lab phenotype (35 R/757 S).

## Novelty claim
First external evaluation of a BDQ triage tool on non-overlapping isolates (KM cohort).

## Caveats / negative results
SENSITIVITY ANALYSIS, not headline proof: KM assay/CC UNDOCUMENTED (MEDIUM quality). Only 35 resistant -> wide sensitivity CI. No Indian external resistant isolates exist. Production model has NOT seen these isolates (zero overlap verified). match != beat.

## Reproduce
```bash
python eval/ext_km_validate.py
```

## Provenance
- **Date (UTC):** 2026-11-19T17:13:55.852724+00:00
- **Git commit:** 22016ec
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `metrics.csv`
- `lineage_stratification.csv`
