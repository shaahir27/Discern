# E07 · Leave-one-site-out (S2)

> **Approach:** `lr_model` · **Test:** `e07_loso` · **PASS ✅**

Pooled sens/spec across held-out sites. Per-site specificity stable; per-site sensitivity unreliable (most sites <25 R).

## Metrics

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `pooled_sens` | 0.7875 | 0.7880 | 0.0200 | ✅ |
| `pooled_spec` | 0.9537 | 0.9540 | 0.0100 | ✅ |

| Metric | Value |
|---|---|
| `per_site_spec_min_n300` | 0.9441 |
| `per_site_spec_max_n300` | 0.9785 |

## Method
Leave-one-site-out (23 sites); nested calibration per training set; pooled binary sens/spec; per-site spec for n>=300; per-site sens only if nR>=25.

## Novelty claim
Cross-site generalization (honest).

## Caveats / negative results
Do NOT report pooled LOSO AUC (mixes 23 calibrations). Site 5 (inferred India, unverified) has ~5 R -> false-alarm audit only.

## Reproduce
```bash
python eval/e07_loso.py
```

## Provenance
- **Date (UTC):** 2026-11-18T19:05:34.579889+00:00
- **Git commit:** uncommitted
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `per_site.csv`
