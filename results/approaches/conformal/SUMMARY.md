# Conformal certificates

> **Type:** method · **Approach:** `conformal` · **0/1 tests pass** ❌

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E14 · Certificate Audit | — | alpha=0.0150 | [report](e14_certificate_audit/REPORT.md) |

## E14 · Certificate Audit — DESCRIPTIVE •

R-zone certified: bound 0.1963 <= beta 0.2 (achieved 0.14342629482071714). S-zone NOT certifiable at alpha=0.015/delta=0.05: best achievable bound 0.0155 (at tS=0.022); the strict NPV>=98.5% certificate fails given the score distribution. Certified policy abstains from S (tS=None), U-zone 630->11617 (widens, as pre-registered). Out-of-site S violations 0/11; cross-source transfer in tables. The bound, not the point target, is the guarantee.

| Metric | Value |
|---|---|
| `alpha` | 0.0150 |
| `beta` | 0.2000 |
| `delta` | 0.0500 |
| `pi_train` | 0.0730 |
| `baseline_emp_s_error` | 0.0146 |
| `baseline_emp_r_error` | 0.2145 |
| `baseline_n_S` | 10688 |
| `baseline_n_R` | 550 |
| `baseline_n_U` | 630 |
| `baseline_n_U_resistant` | 278 |
| `certified_t_R` | 0.8766 |
| `certified_bound_s` | 1.0000 |
| `certified_bound_r` | 0.1963 |
| `certified_emp_r_error` | 0.1434 |
| `certified_n_S` | 0 |
| `certified_n_R` | 251 |
| `certified_n_U` | 11617 |
| `certified_n_U_resistant` | 651 |
| `fold_n_tS_certified` | 0 |
| `fold_emp_r_min` | 0.1636 |
| `fold_emp_r_max` | 0.2131 |
| `s_best_achievable_bound` | 0.0155 |
| `s_best_achievable_t` | 0.0220 |
| `out_of_site_n_sites` | 11 |
| `out_of_site_violations_s` | 0 |
| `out_of_site_violations_r` | 5 |
| `cs_CRyPTIC_to_NICD_bound_s` | 0.0097 |
| `cs_CRyPTIC_to_NICD_bound_r` | 1.0000 |
| `cs_CRyPTIC_to_NICD_emp_s_error` | 0.6713 |
| `cs_NICD_to_CRyPTIC_bound_s` | 1.0000 |
| `cs_NICD_to_CRyPTIC_bound_r` | 0.1791 |
| `cs_NICD_to_CRyPTIC_emp_r_error` | 0.9772 |
| `weighted_emp_r_error_pi002` | 0.3924 |
| `weighted_emp_r_error_pi005` | 0.2003 |

- **certified_t_S:** `None`
- **certified_s_certified:** `False`
- **certified_r_certified:** `True`
- **certified_emp_s_error:** `None`
- **fold_tS_min:** `None`
- **fold_tS_max:** `None`
- **fold_emp_s_min:** `None`
- **fold_emp_s_max:** `None`
- **cs_CRyPTIC_to_NICD_emp_r_error:** `None`
- **cs_CRyPTIC_to_NICD_violates_s:** `True`
- **cs_CRyPTIC_to_NICD_violates_r:** `False`
- **cs_NICD_to_CRyPTIC_emp_s_error:** `None`
- **cs_NICD_to_CRyPTIC_violates_s:** `False`
- **cs_NICD_to_CRyPTIC_violates_r:** `True`
- **weighted_emp_s_error_pi002:** `None`
- **weighted_emp_s_error_pi005:** `None`

_Caveats: Patient grouping => effective n ~ patients; marginal split-conformal validity only, not i.i.d. R-zone bound is wide (~170 positives/fold). Prevalence-weighted certificate is approximate at an estimated user prevalence (Barber 2022 coverage-gap). Exchangeability breaks across labs — violations are a RESULT, motivating the per-lab prevalence slider._

