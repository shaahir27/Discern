# MIC label-noise audit

> **Type:** analysis · **Approach:** `mic_audit` · **0/1 tests pass** ❌

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E17 · Borderline MIC | — | n_borderline=454 | [report](e17_borderline_mic/REPORT.md) |

## E17 · Borderline MIC — DESCRIPTIVE •

454 borderline samples (191 resistant, 0.221 of all positives). Excluding them shifts AUC 0.9085->0.9045 (delta -0.0041). Label noise near the cutoff is modest; sensitivity analysis only.

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

_Caveats: SENSITIVITY ANALYSIS, not a model change; production keeps all rows. Censored MICs (<=,>) excluded from borderline. NICD MGIT MICs largely unavailable. Margin-weighted fit NOT promoted from this audit (pre-registered)._

