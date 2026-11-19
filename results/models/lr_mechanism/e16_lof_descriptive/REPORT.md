# E16 · Lof Descriptive

> **Approach:** `lr_mechanism` · **Test:** `e16_lof_descriptive` · **DESCRIPTIVE •**

Snobre's susceptible-frameshift fraction REPLICATES in this lineage mix: 0.2504 (162/647) vs reported 0.248. DBD-break truncations dominate; 43 late variants are function-preserving candidates.

## Metrics

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

## Method
Rv0678 CDS (NC_000962.3:778,990-779,487) translation; first-stop region DBD<=99 / mid 100-160 / late>160; high-confidence LoF excludes single-nt homopolymer indels; late-stop preserves dimerization domain (Snobre function-preserving proxy).

## Novelty claim
First predictive integration of truncation-position reasoning for BDQ (to our knowledge); descriptive replication of Snobre 2024.

## Caveats / negative results
ARF conservation approximated by late-stop preservation (full phylogenetic ARF alignment out of scope). Homopolymer indels flagged low-confidence (e.g. 193_del_g). Whole-gene deletion = full_del. A sample carrying multiple regions is counted in each.

## Reproduce
```bash
python eval/e16_mechanism.py
```

## Provenance
- **Date (UTC):** 2026-11-19T04:18:59.720443+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `lof_effects.csv`
- `region_summary.csv`
