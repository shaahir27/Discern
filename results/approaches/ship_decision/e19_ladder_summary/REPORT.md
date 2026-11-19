# E19 · Ladder Summary

> **Approach:** `ship_decision` · **Test:** `e19_ladder_summary` · **DESCRIPTIVE •**

Keep golden L2 LR (F2). Adopt the audited Tier A veto modifier.

## Metrics

- **model:** `models/lr_model`
- **epistatic_veto:** `True`

## Method
Apply D011 to the saved E16, shrinkage, conformal, and E21 records and the KM veto audit.

## Caveats / negative results
Research use only. KM sensitivity analysis with 35 resistant isolates, undocumented assay and critical concentration, and medium-quality phenotypes. Report the certificate bounds with their confidence level and lab-shift results.

## Reproduce
```bash
C:/Users/shaah/miniconda3/envs/discern/python.exe eval/e19_ship_decision.py
```

## Provenance
- **Date (UTC):** 2026-11-19T06:00:12.636531+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `ladder.csv`
