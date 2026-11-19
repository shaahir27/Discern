# E KM Cohort

> **Approach:** `external_km` · **Test:** `e_km_cohort` · **DESCRIPTIVE •**

Verified external cohort: 792 isolates, 35 resistant (prevalence 0.0442), ZERO training overlap. All MEDIUM quality; METHOD_3/METHOD_CC undocumented. Lineage mix L2=559/L4=198. Per P-EXT this is a SENSITIVITY ANALYSIS, not headline proof. Genotype->phenotype sens/spec is BLOCKED: external variant calls are not local.

## Metrics

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

## Method
KMVALIDATION2025 BDQ DST from CRyPTIC v3.4.0 (DST_MEASUREMENTS); lineage from GENOMES; training-ID overlap asserted against the 11,868-sample cohort.

## Novelty claim
First external evaluation cohort for a BDQ triage tool (characterized; accuracy pending external genotypes).

## Caveats / negative results
Assay/CC UNDOCUMENTED (MEDIUM) -> cannot assume WHO critical concentration; stratify by lineage. NO external sens/spec computed (no external genotypes). NO Indian external resistant isolates exist -> no Indian-sensitivity claim ever. To finish: acquire CRyPTIC MUTATIONS (BDQ genes) for these 792 isolates, then run catalogue+model.

## Reproduce
```bash
python eval/ext_validate.py
```

## Provenance
- **Date (UTC):** 2026-11-19T04:32:30.703816+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `lineage_stratification.csv`
