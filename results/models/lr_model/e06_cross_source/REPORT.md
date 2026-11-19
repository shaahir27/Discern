# E06 · Cross-source (S3)

> **Approach:** `lr_model` · **Test:** `e06_cross_source` · **PASS ✅**

Accuracy depends on lab/source: NICD->CRyPTIC drops (few CRyPTIC resistant training germs, lineage/prevalence shift).

## Metrics

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `CRyPTIC_to_NICD_auc` | 0.9105 | 0.9100 | 0.0300 | ✅ |
| `CRyPTIC_to_NICD_sens` | 0.8442 | 0.8440 | 0.0300 | ✅ |
| `CRyPTIC_to_NICD_spec` | 0.9633 | 0.9630 | 0.0300 | ✅ |
| `NICD_to_CRyPTIC_auc` | 0.8277 | 0.8300 | 0.0300 | ✅ |
| `NICD_to_CRyPTIC_sens` | 0.6747 | 0.6750 | 0.0300 | ✅ |
| `NICD_to_CRyPTIC_spec` | 0.9532 | 0.9530 | 0.0300 | ✅ |

## Method
Train on one SRC, test on the other; nested calibration inside training source; binary at per-fold Youden cut.

## Novelty claim
Cross-lab generalization limit (honest).

## Caveats / negative results
NICD test has 783 R; CRyPTIC test only 83 R (wide CI).

## Reproduce
```bash
python eval/e06_cross_source.py
```

## Provenance
- **Date (UTC):** 2026-11-18T19:05:30.083030+00:00
- **Git commit:** uncommitted
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `cross_source.csv`
