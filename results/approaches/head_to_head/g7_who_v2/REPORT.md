# G7 · WHO V2

> **Approach:** `head_to_head` · **Test:** `g7_who_v2` · **DESCRIPTIVE •**

WHO v2 ternary sens_firm at FRS0.25 = 0.7680; our model ternary sens among firm = 0.7347. Our ternary does NOT beat WHO at matched FRS (honest). Binary U-as-S sensitivities reported at each FRS convention.

## Metrics

| Metric | Value |
|---|---|
| `who_frs0.1_firm` | 0.9607 |
| `who_frs0.1_sens_firm` | 0.8011 |
| `who_frs0.1_spec_firm` | 0.9812 |
| `who_frs0.1_sens_U_as_S` | 0.6790 |
| `who_frs0.1_sens_U_as_R` | 0.8314 |
| `who_frs0.25_firm` | 0.9611 |
| `who_frs0.25_sens_firm` | 0.7680 |
| `who_frs0.25_spec_firm` | 0.9828 |
| `who_frs0.25_sens_U_as_S` | 0.6420 |
| `who_frs0.25_sens_U_as_R` | 0.8060 |
| `who_frs0.75_firm` | 0.9680 |
| `who_frs0.75_sens_firm` | 0.5606 |
| `who_frs0.75_spec_firm` | 0.9868 |
| `who_frs0.75_sens_U_as_S` | 0.4965 |
| `who_frs0.75_sens_U_as_R` | 0.6109 |
| `our_catalogue_grouped_sens_firm` | 0.7924 |
| `our_catalogue_grouped_spec_firm` | 0.9835 |
| `our_model_ternary_sens_among_firm` | 0.7347 |
| `our_model_ternary_spec_among_firm` | 0.9889 |

## Method
WHO v2 GARC catalogue as-is via piezo at FRS 0.1/0.25/0.75 on identical 11,868; ternary R>U>S. Compared to our patient-grouped catalogue (E03) and model ternary (E04).

## Novelty claim
Same-data head-to-head vs WHO v2 (P-G7), honest matched-FC comparison.

## Caveats / negative results
Same-data comparison only; NEVER compare AUC across cohorts. Ternary sens among firm is descriptive. WHO FRS convention matters (0.25 vs 0.1 vs 0.75). match != beat.

## Reproduce
```bash
python eval/e18_head_to_head.py
```

## Provenance
- **Date (UTC):** 2026-11-19T04:35:58.406100+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `who_vs_ours.csv`
