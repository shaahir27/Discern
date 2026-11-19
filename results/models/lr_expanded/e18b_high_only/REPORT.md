# E18b High Only

> **Approach:** `lr_expanded` · **Test:** `e18b_high_only` · **DESCRIPTIVE •**

HIGH-only expansion (16923/201 R): S1 AUC 0.6668 (golden 0.9085, HIGH+MEDIUM 0.8381); KM external AUC 0.8633 (golden production 0.8682). Cross-source: N/A (NICD is MEDIUM, excluded from HIGH-only).

## Metrics

| Metric | Value |
|---|---|
| `n_samples` | 16923 |
| `n_resistant` | 201 |
| `auc` | 0.6668 |
| `binary_sens` | 0.4030 |
| `binary_spec` | 0.9637 |
| `sens_at_spec_0.95` | 0.4080 |
| `sens_at_spec_0.985` | 0.2687 |
| `golden_auc` | 0.9085 |
| `high_medium_auc` | 0.8381 |
| `km_external_auc` | 0.8633 |
| `golden_production_km_auc` | 0.8682 |

- **cross_source_note:** `NICD is MEDIUM quality; absent from HIGH-only arm -> cross-source N/A`

## Method
HIGH-quality-only expanded CRyPTIC v3.4 BDQ cohort (exclude KM). LR F2, S1 patient-grouped nested CV + cross-source + KM external. Isolates data quantity from MEDIUM label noise.

## Novelty claim
Data-richness ablation: quality-controlled expansion vs golden.

## Caveats / negative results
Descriptive; different cohort than golden reproduction. KM assay undocumented; 35 R wide CI. Golden + production untouched; no tuning. Tests if CLEAN richer data helps generalization.

## Reproduce
```bash
python eval/e18b_high_only.py
```

## Provenance
- **Date (UTC):** 2026-11-19T05:07:32.696575+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `metrics.csv`
