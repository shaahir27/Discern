# E10 · Prevalence-shift correction

> **Approach:** `prevalence_correction` · **Test:** `e10_prevalence_shift` · **PASS ✅**

Calibrated probabilities do NOT transfer across prevalence unless corrected. NICD->CRyPTIC ECE 0.288->0.013; reverse improves less (only 83 CRyPTIC R).

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
Cross-source OOF probabilities; odds-ratio prevalence correction (section 10); ECE (10 quantile bins) + Brier raw vs adjusted.

## Novelty claim
Prevalence correction validated across labs (math is Saerens 2002; application here).

## Caveats / negative results
At low prevalence the PPV-80% R threshold is rarely reached — correct behaviour.

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
