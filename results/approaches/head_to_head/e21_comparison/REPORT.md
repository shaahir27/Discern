# E21 · Comparison

> **Approach:** `head_to_head` · **Test:** `e21_comparison` · **DESCRIPTIVE •**

Keep golden LR (F2). No tuned candidate met both D013 criteria.

## Metrics

- **S1_leader:** `rf_tuned`
- **KM_leader:** `['hgb_tuned']`

## Method
Rank S1 and KM separately, using the saved golden LR references. D013 requires an S1 AUC gain greater than 0.005 and KM AUC at least as high as golden. The comparison shows the settings selected on all golden training for the KM evaluation.

## Caveats / negative results
Research use only. The selected settings and Platt calibration use the same inner folds. The outer and KM results are held-out evaluations. KM is a sensitivity analysis with undocumented assay and critical concentration, medium-quality phenotypes, and few resistant isolates. Sensitivity at fixed specificity describes the held-out ROC curve.

## Reproduce
```bash
C:/Users/shaah/miniconda3/envs/discern/python.exe eval/e21_tuned_comparison.py
```

## Provenance
- **Date (UTC):** 2026-11-19T05:32:05.318019+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `comparison.csv`
