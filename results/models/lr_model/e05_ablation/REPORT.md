# E05 · Feature ablation

> **Approach:** `lr_model` · **Test:** `e05_ablation` · **PASS ✅**

Type/region features ~+0.11 AUC; minor alleles ~+0.11 AUC; mmpL5 interaction ~+0.004 AUC.

## Metrics

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `F0_frs0.1_auc` | 0.7910 | 0.7910 | 0.0150 | ✅ |
| `F1_frs0.1_auc` | 0.9045 | 0.9050 | 0.0150 | ✅ |
| `F2_frs0.1_auc` | 0.9085 | 0.9090 | 0.0150 | ✅ |
| `F2_frs0.75_auc` | 0.7986 | 0.7990 | 0.0150 | ✅ |
| `F0_frs0.1_auprc` | 0.5026 | 0.5030 | 0.0150 | ✅ |
| `F1_frs0.1_auprc` | 0.6721 | 0.6720 | 0.0150 | ✅ |
| `F2_frs0.1_auprc` | 0.6836 | 0.6840 | 0.0150 | ✅ |
| `F2_frs0.75_auprc` | 0.5479 | 0.5480 | 0.0150 | ✅ |

## Method
S1 nested CV; F0 vocab only, F1 +type/region, F2 +mmpL5/interaction; F2 frs_min 0.75 ignores minor alleles.

## Novelty claim
Quantifies feature value (minor alleles published by Fowler lab).

## Caveats / negative results
Ablation, not tuning. frs 0.75 drops minor-allele signal.

## Reproduce
```bash
python eval/e05_ablation.py
```

## Provenance
- **Date (UTC):** 2026-11-18T19:05:27.640111+00:00
- **Git commit:** uncommitted
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `ablation.csv`
