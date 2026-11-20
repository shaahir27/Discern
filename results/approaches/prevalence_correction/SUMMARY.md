# Prevalence correction

> **Type:** analysis · **Approach:** `prevalence_correction` · **1/1 tests pass** ✅

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E10 · Calibration with known target prevalence | ✅ | CRyPTIC_to_NICD_ece_raw=0.5684 | [report](e10_prevalence_shift/REPORT.md) |

## E10 · Calibration with known target prevalence — PASS ✅

Adjusting to the observed prevalence improved ECE and Brier score in both archived comparisons.

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

_Caveats: This comparison used known target prevalence. It does not validate estimating prevalence from an unlabeled batch or transferring calibration to a new lab._

