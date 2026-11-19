# E18 · KM External

> **Approach:** `lr_expanded` · **Test:** `e18_km_external` · **DESCRIPTIVE •**

Expanded model external KM AUC 0.8636 vs golden production 0.8682 (delta -0.0046). Tests whether richer training improves generalization to unseen isolates.

## Metrics

| Metric | Value |
|---|---|
| `n_isolates` | 792 |
| `n_resistant` | 35 |
| `expanded_model_auc` | 0.8636 |
| `expanded_model_sens` | 0.8571 |
| `expanded_model_spec` | 0.7728 |
| `expanded_sens_at_spec_0.95` | 0.6571 |
| `golden_production_km_auc` | 0.8682 |
| `delta_auc_vs_golden_production` | -0.0046 |

## Method
Expanded-trained LR applied to held-out KM 792 isolates (zero overlap with expanded training). Compare to golden-production KM AUC 0.868.

## Novelty claim
External generalization benefit of richer training data.

## Caveats / negative results
KM assay undocumented (MEDIUM); 35 R -> wide CI; sensitivity analysis. No tuning; KM excluded from expanded training. match != beat.

## Reproduce
```bash
python eval/e18_expanded_training.py
```

## Provenance
- **Date (UTC):** 2026-11-19T05:05:31.803098+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `metrics.csv`
