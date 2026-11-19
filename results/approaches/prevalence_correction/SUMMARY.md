# Prevalence correction

> **Type:** analysis · **Approach:** `prevalence_correction` · **1/1 tests pass** ✅

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E10 · Prevalence-shift correction | ✅ | CRyPTIC_to_NICD_ece_raw=0.5684 | [report](e10_prevalence_shift/REPORT.md) |

## E10 · Prevalence-shift correction — PASS ✅

Calibrated probabilities do NOT transfer across prevalence unless corrected. NICD->CRyPTIC ECE 0.288->0.013; reverse improves less (only 83 CRyPTIC R).

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

_Caveats: At low prevalence the PPV-80% R threshold is rarely reached — correct behaviour._

