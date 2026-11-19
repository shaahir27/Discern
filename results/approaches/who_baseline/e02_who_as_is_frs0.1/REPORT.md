# E02 · WHO v2 as-is (FRS ≥ 0.1)

> **Approach:** `who_baseline` · **Test:** `e02_who_as_is_frs0.1` · **PASS ✅**

WHO v2 as-is at FRS>=0.1, pooled. Anchor sens_U_as_S must match published.

## Metrics

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `firm` | 0.9607 | 0.9610 | 0.0050 | ✅ |
| `sens_firm` | 0.8011 | 0.8010 | 0.0050 | ✅ |
| `spec_firm` | 0.9812 | 0.9810 | 0.0050 | ✅ |
| `sens_U_as_S` | 0.6790 | 0.6790 | 0.0050 | ✅ |

| Metric | Value |
|---|---|
| `spec_U_as_S` | 0.9817 |
| `sens_U_as_R` | 0.8314 |
| `spec_U_as_R` | 0.9514 |
| `n_U` | 466 |
| `n_U_resistant` | 132 |

- **variant_call_counts:** `{'U': 467, 'S': 346, 'R': 229}`

## Method
WHO v2 GARC catalogue via piezo; non-syn/non-null variants FRSf>=0.1; exception->S; aggregate R>U>S.

## Novelty claim
Baseline reproduction (published WHO v2 performance).

## Caveats / negative results
U counted as S for the anchor sensitivity. NICD enriches resistance; read per-source.

## Reproduce
```bash
python eval/e02_who.py
```

## Provenance
- **Date (UTC):** 2026-11-18T19:04:42.674378+00:00
- **Git commit:** uncommitted
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `metrics.csv`
