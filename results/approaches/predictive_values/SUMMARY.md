# Predictive values (PPV/NPV)

> **Type:** analysis · **Approach:** `predictive_values` · **1/1 tests pass** ✅

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E13 · PPV/NPV at prevalence | ✅ | model_ppv_0.01=0.1436 | [report](e13_ppv_npv/REPORT.md) |

## E13 · PPV/NPV at prevalence — PASS ✅

At 2% prevalence ~half of catalogue R calls are false alarms (PPV 50%). R is a flag for confirmation; S is reliable (NPV ~99.6%).

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `model_ppv_0.01` | 0.1436 | 0.1400 | 0.0200 | ✅ |
| `model_npv_0.01` | 0.9982 | 0.9980 | 0.0030 | ✅ |
| `cat_ppv_0.01` | 0.3328 | 0.3300 | 0.0200 | ✅ |
| `cat_npv_0.01` | 0.9978 | 0.9980 | 0.0030 | ✅ |
| `model_ppv_0.02` | 0.2530 | 0.2500 | 0.0200 | ✅ |
| `model_npv_0.02` | 0.9964 | 0.9960 | 0.0030 | ✅ |
| `cat_ppv_0.02` | 0.5019 | 0.5000 | 0.0200 | ✅ |
| `cat_npv_0.02` | 0.9957 | 0.9960 | 0.0030 | ✅ |
| `model_ppv_0.05` | 0.4663 | 0.4700 | 0.0200 | ✅ |
| `model_npv_0.05` | 0.9907 | 0.9910 | 0.0030 | ✅ |
| `cat_ppv_0.05` | 0.7221 | 0.7200 | 0.0200 | ✅ |
| `cat_npv_0.05` | 0.9889 | 0.9890 | 0.0030 | ✅ |
| `model_ppv_0.1` | 0.6484 | 0.6500 | 0.0200 | ✅ |
| `model_npv_0.1` | 0.9805 | 0.9810 | 0.0030 | ✅ |
| `cat_ppv_0.1` | 0.8458 | 0.8500 | 0.0200 | ✅ |
| `cat_npv_0.1` | 0.9768 | 0.9770 | 0.0030 | ✅ |

_Caveats: Illustrative operating points; real lab prevalence must replace the placeholder._

