# E11 · Imbalance + learning curve

> **Approach:** `lr_model` · **Test:** `e11_imbalance_learning_curve` · **PASS ✅**

Balancing does not help and breaks calibration; AUC is flat across resistant fractions (signal saturates early).

## Metrics

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `lc_auc_0.1` | 0.8963 | 0.8980 | 0.0100 | ✅ |
| `lc_auc_0.25` | 0.9050 | 0.9050 | 0.0100 | ✅ |
| `lc_auc_0.5` | 0.9050 | 0.9080 | 0.0100 | ✅ |
| `lc_auc_0.75` | 0.9075 | 0.9060 | 0.0100 | ✅ |
| `lc_auc_1.0` | 0.9085 | 0.9090 | 0.0100 | ✅ |
| `balanced_auc` | 0.9051 | 0.9040 | 0.0100 | ✅ |
| `balanced_auprc` | 0.6587 | 0.6530 | 0.0200 | ✅ |
| `balanced_sens_at_spec98` | 0.6062 | 0.5920 | 0.0300 | ✅ |

## Method
1:1 class_weight balanced LR (NOT production); learning curve retains fraction of resistant training rows (all S kept), nested calibration, S1.

## Novelty claim
Justifies no class_weight / no resampling (honest).

## Caveats / negative results
Balanced model is diagnostic only; production model uses no class_weight (rule 4).

## Reproduce
```bash
python eval/e11_imbalance.py
```

## Provenance
- **Date (UTC):** 2026-11-18T19:05:46.838900+00:00
- **Git commit:** uncommitted
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `learning_curve.csv`
