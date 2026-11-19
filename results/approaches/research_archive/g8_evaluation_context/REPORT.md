# G8 · Evaluation Context

> **Approach:** `research_archive` · **Test:** `g8_evaluation_context` · **DESCRIPTIVE •**

## Metrics

| Metric | Value |
|---|---|
| `interval_level` | 0.9500 |
| `bootstrap_replicates` | 500 |
| `n_resistant` | 866 |
| `n_resistant_without_retained_core_variant` | 137 |
| `fraction_without_retained_core_variant` | 0.1582 |

## Method
Patient resampling of the archived S1 predictions and Wilson intervals from saved KM firm-call totals.

## Caveats / negative results
Research use only. For evaluation and education. Clinical decisions require laboratory confirmation. Resampling fixed predictions measures sample uncertainty. Refitting variability is outside these intervals.

## Reproduce
```bash
python eval/g8_context.py
```

## Provenance
- **Date (UTC):** 2026-11-19T06:22:19.415498+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `intervals.csv`
