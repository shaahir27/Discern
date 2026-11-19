# E18 · Expanded CV

> **Approach:** `lr_expanded` · **Test:** `e18_expanded_cv` · **DESCRIPTIVE •**

Expanded training AUC 0.8381 vs golden 0.9085 (delta -0.0704); sens@spec0.985 0.6195. +182 new resistant isolates, CRyPTIC R 83->255.

## Metrics

| Metric | Value |
|---|---|
| `n_samples` | 23441 |
| `n_resistant` | 1046 |
| `auc` | 0.8381 |
| `auprc` | 0.5665 |
| `binary_sens` | 0.7247 |
| `binary_spec` | 0.9650 |
| `sens_at_spec_0.95` | 0.7266 |
| `sens_at_spec_0.97` | 0.7132 |
| `sens_at_spec_0.985` | 0.6195 |
| `n_features` | 137 |
| `golden_auc` | 0.9085 |
| `delta_auc_vs_golden` | -0.0704 |

## Method
Expanded CRyPTIC v3.4 BDQ cohort (HIGH+MEDIUM, exclude KM): 23,441 isolates / 1,046 R. LR F2, S1 patient-grouped 5-fold nested Platt + thresholds (same protocol as golden E04).

## Novelty claim
Data-richness lever: more resistant isolates to improve generalization.

## Caveats / negative results
Different cohort than golden (not the published 11,868 reproduction). MEDIUM labels add noise. Descriptive; golden baseline + production untouched; no tuning. KM held out.

## Reproduce
```bash
python eval/e18_expanded_training.py
```

## Provenance
- **Date (UTC):** 2026-11-19T05:05:30.127477+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `metrics.csv`
- `thresholds_per_fold.csv`
