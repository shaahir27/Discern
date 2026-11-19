# Tier B — mechanism-enhanced

> **Type:** model · **Approach:** `lr_mechanism` · **0/2 tests pass** ❌

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E16 · Lof Ablation | — | lr_f2_auc=0.9085 | [report](e16_lof_ablation/REPORT.md) |
| E16 · Lof Descriptive | — | snobre_frac_susceptible_frameshift=0.2504 | [report](e16_lof_descriptive/REPORT.md) |

## E16 · Lof Ablation — DESCRIPTIVE •

Mechanism features delta AUC +0.0002 (sens@spec0.985 +0.0127). S1 veto downgrades 41 R->U: 2 true resistant (harm), 39 true susceptible (benefit). Ship decision DEFERRED to G5 (must not degrade external validation); golden kept unless a clear patient-grouped improvement survives external check.

| Metric | Value |
|---|---|
| `lr_f2_auc` | 0.9085 |
| `lr_f2mech_auc` | 0.9088 |
| `lr_delta_auc` | 0.0002 |
| `lr_f2_sens_at_spec_0.985` | 0.6039 |
| `lr_f2mech_sens_at_spec_0.985` | 0.6166 |
| `cat_s1_sens_firm_off` | 0.7924 |
| `cat_s1_sens_firm_on` | 0.7918 |
| `cat_s1_spec_firm_off` | 0.9835 |
| `cat_s1_spec_firm_on` | 0.9872 |
| `cat_s1_n_U_off` | 661 |
| `cat_s1_n_U_on` | 702 |
| `s1_n_flips` | 41 |
| `s1_n_flips_resistant` | 2 |
| `s1_n_flips_susceptible` | 39 |
| `cs_CRyPTIC_to_NICD_sens_firm_off` | 0.6078 |
| `cs_CRyPTIC_to_NICD_sens_firm_on` | 0.6078 |
| `cs_CRyPTIC_to_NICD_spec_firm_off` | 0.9919 |
| `cs_CRyPTIC_to_NICD_spec_firm_on` | 0.9919 |
| `cs_CRyPTIC_to_NICD_n_flips` | 0 |
| `cs_CRyPTIC_to_NICD_n_flips_resistant` | 0 |
| `cs_NICD_to_CRyPTIC_sens_firm_off` | 0.5517 |
| `cs_NICD_to_CRyPTIC_sens_firm_on` | 0.5517 |
| `cs_NICD_to_CRyPTIC_spec_firm_off` | 0.9781 |
| `cs_NICD_to_CRyPTIC_spec_firm_on` | 0.9818 |
| `cs_NICD_to_CRyPTIC_n_flips` | 39 |
| `cs_NICD_to_CRyPTIC_n_flips_resistant` | 0 |

_Caveats: Veto is conservative (R->U, send-for-lab, not R->S). 'No atpE variant' is a proxy for 'not atpE-driven R'. Homopolymer indels excluded from high-confidence LoF. Mechanism features may move metrics negligibly (few truncation carriers) — a valid outcome. External validation (G5) pending; no ship decision yet._

## E16 · Lof Descriptive — DESCRIPTIVE •

Snobre's susceptible-frameshift fraction REPLICATES in this lineage mix: 0.2504 (162/647) vs reported 0.248. DBD-break truncations dominate; 43 late variants are function-preserving candidates.

| Metric | Value |
|---|---|
| `snobre_frac_susceptible_frameshift` | 0.2504 |
| `snobre_reported` | 0.2480 |
| `n_frameshift_carriers` | 647 |
| `n_frameshift_susceptible` | 162 |
| `n_truncating_variants` | 142 |
| `n_veto_fire_samples` | 43 |
| `n_variants_DBD` | 65 |
| `n_variants_mid` | 32 |
| `n_variants_late` | 43 |
| `n_variants_full_del` | 2 |

_Caveats: ARF conservation approximated by late-stop preservation (full phylogenetic ARF alignment out of scope). Homopolymer indels flagged low-confidence (e.g. 193_del_g). Whole-gene deletion = full_del. A sample carrying multiple regions is counted in each._

