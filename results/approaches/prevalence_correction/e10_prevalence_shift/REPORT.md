# E10 · Calibration with known target prevalence

> **Approach:** `prevalence_correction` · **Test:** `e10_prevalence_shift` · **PASS ✅**

Adjusting to the observed prevalence improved ECE and Brier score in both archived comparisons.

## Metrics

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `CRyPTIC_to_NICD_ece_raw` | 0.5684 | 0.5680 | 0.0300 | ✅ |
| `CRyPTIC_to_NICD_ece_adj` | 0.1012 | 0.1010 | 0.0300 | ✅ |
| `CRyPTIC_to_NICD_brier_raw` | 0.5009 | 0.5010 | 0.0200 | ✅ |
| `CRyPTIC_to_NICD_brier_adj` | 0.1238 | 0.1240 | 0.0200 | ✅ |
| `NICD_to_CRyPTIC_ece_raw` | 0.2882 | 0.2880 | 0.0100 | ✅ |
| `NICD_to_CRyPTIC_ece_adj` | 0.0128 | 0.0130 | 0.0100 | ✅ |
| `NICD_to_CRyPTIC_brier_raw` | 0.1020 | 0.1020 | 0.0050 | ✅ |
| `NICD_to_CRyPTIC_brier_adj` | 0.0169 | 0.0170 | 0.0050 | ✅ |

## Method
Cross-source out-of-fold probabilities adjusted to known target prevalence. The target prevalence came from recorded lab labels. Calibration was measured using quantile-binned ECE and Brier score.

## Novelty claim
Retrospective calibration comparison with known target prevalence.

## Caveats / negative results
This comparison used known target prevalence. It does not validate estimating prevalence from an unlabeled batch or transferring calibration to a new lab.

## Reproduce
```bash
python eval/e10_prevalence_shift.py
```

## Provenance
- **Date (UTC):** 2026-11-18T19:05:41.528065+00:00
- **Git commit:** uncommitted
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `prevalence_shift.csv`
