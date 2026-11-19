# E16 · Lof Ablation

> **Approach:** `lr_mechanism` · **Test:** `e16_lof_ablation` · **DESCRIPTIVE •**

Mechanism features delta AUC +0.0002 (sens@spec0.985 +0.0127). S1 veto downgrades 41 R->U: 2 true resistant (harm), 39 true susceptible (benefit). Ship decision DEFERRED to G5 (must not degrade external validation); golden kept unless a clear patient-grouped improvement survives external check.

## Metrics

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

## Method
S1 nested LR F2 vs F2+mechanism (rv_lof_any/highconf/dbd/late). Catalogue epistatic veto R->U when high-confidence Rv0678 LoF + mmpL5/S5 LoF + no atpE variant. S1 + cross-source. Every flip enumerated with true label.

## Novelty claim
Audited epistatic downgrade rule + mechanism-feature ablation (novelty P-E16).

## Caveats / negative results
Veto is conservative (R->U, send-for-lab, not R->S). 'No atpE variant' is a proxy for 'not atpE-driven R'. Homopolymer indels excluded from high-confidence LoF. Mechanism features may move metrics negligibly (few truncation carriers) — a valid outcome. External validation (G5) pending; no ship decision yet.

## Reproduce
```bash
python eval/e16_mechanism.py
```

## Provenance
- **Date (UTC):** 2026-11-19T04:19:24.547556+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `lr_compare.csv`
- `s1_veto_metrics.csv`
- `cross_source_veto.csv`
- `call_flips.csv`
