# E17 · Borderline MIC

> **Approach:** `mic_audit` · **Test:** `e17_borderline_mic` · **DESCRIPTIVE •**

454 borderline samples (191 resistant, 0.221 of all positives). Excluding them shifts AUC 0.9085->0.9045 (delta -0.0041). Label noise near the cutoff is modest; sensitivity analysis only.

## Metrics

| Metric | Value |
|---|---|
| `n_borderline` | 454 |
| `n_borderline_resistant` | 191 |
| `borderline_share_of_positives` | 0.2206 |
| `n_cry_borderline` | 236 |
| `n_nicd_borderline` | 218 |
| `n_nicd_mic_available` | 612 |
| `full_auc` | 0.9085 |
| `full_binary_sens` | 0.8349 |
| `full_binary_spec` | 0.9481 |
| `full_sens_at_spec_0.985` | 0.6039 |
| `full_n` | 11868 |
| `full_n_resistant` | 866 |
| `excl_borderline_auc` | 0.9045 |
| `excl_borderline_binary_sens` | 0.8237 |
| `excl_borderline_binary_spec` | 0.9556 |
| `excl_borderline_sens_at_spec_0.985` | 0.6578 |
| `excl_borderline_n` | 11414 |
| `excl_borderline_n_resistant` | 675 |
| `delta_auc_excl_minus_full` | -0.0041 |

## Method
Borderline = exact MIC within +/-1 two-fold dilution of source CC (CRyPTIC 0.25, MGIT 1.0). E04 OOF metrics recomputed excluding borderline samples, side-by-side. Production keeps ALL rows.

## Novelty claim
Label-noise audit (P-E17): quantifies cutoff-adjacent uncertainty honestly.

## Caveats / negative results
SENSITIVITY ANALYSIS, not a model change; production keeps all rows. Censored MICs (<=,>) excluded from borderline. NICD MGIT MICs largely unavailable. Margin-weighted fit NOT promoted from this audit (pre-registered).

## Reproduce
```bash
python eval/e17_mic_audit.py
```

## Provenance
- **Date (UTC):** 2026-11-19T04:35:00.638613+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `borderline_by_mic.csv`
- `full_vs_excluded.csv`
