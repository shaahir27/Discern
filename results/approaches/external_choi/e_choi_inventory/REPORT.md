# E Choi Inventory

> **Approach:** `external_choi` · **Test:** `e_choi_inventory` · **DESCRIPTIVE •**

49 South Korean run accessions inventoried. Paper breakpoint >=0.5 mg/L differs from WHO CC (MGIT 1.0 / plate 0.25); phenotypes must be reclassified first.

## Metrics

| Metric | Value |
|---|---|
| `n_accessions` | 49 |

- **country:** `South Korea`
- **breakpoint_caveat:** `paper breakpoint >=0.5 != WHO CC; reclassify at 1.0 / 0.25 first`
- **genotype_pipeline_blocked:** `True`

## Method
Choi 2026 Korea run-accession inventory (data/choi_runs.tsv). Breakpoint differs from WHO critical concentration; reclassification required before comparison.

## Novelty claim
External validation cohort (Choi Korea), pending genotype acquisition + reclassification.

## Caveats / negative results
BLOCKED: needs ENA PRJNA1301050 raw reads -> variant calling + MIC values for reclassification. Small n (49) -> wide intervals. No accuracy fabricated.

## Reproduce
```bash
python eval/ext_validate.py
```

## Provenance
- **Date (UTC):** 2026-11-19T04:32:30.902968+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `accessions.csv`
