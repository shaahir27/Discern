# E18 · Cross Source

> **Approach:** `lr_expanded` · **Test:** `e18_cross_source` · **DESCRIPTIVE •**

Cross-source generalization with richer CRyPTIC resistant data (255 R vs 83). Compare NICD->CRyPTIC sensitivity to golden 0.675 (the weak spot).

## Metrics

| Metric | Value |
|---|---|
| `CRyPTIC_to_NICD_auc` | 0.9020 |
| `CRyPTIC_to_NICD_sens` | 0.8557 |
| `CRyPTIC_to_NICD_spec` | 0.8707 |
| `CRyPTIC_to_NICD_golden_auc` | 0.9105 |
| `NICD_to_CRyPTIC_auc` | 0.6533 |
| `NICD_to_CRyPTIC_sens` | 0.3320 |
| `NICD_to_CRyPTIC_spec` | 0.9667 |
| `NICD_to_CRyPTIC_golden_auc` | 0.8277 |

## Method
Expanded cohort, train one source test the other (CRyPTIC<->NICD), nested calibration.

## Novelty claim
Tests whether more diverse resistant data fixes cross-source generalization.

## Caveats / negative results
Descriptive; different cohort. CRyPTIC test R count larger than golden (83). No tuning.

## Reproduce
```bash
python eval/e18_expanded_training.py
```

## Provenance
- **Date (UTC):** 2026-11-19T05:05:30.882435+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `cross_source.csv`
