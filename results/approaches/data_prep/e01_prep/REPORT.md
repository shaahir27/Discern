# E01 · Data preparation

> **Approach:** `data_prep` · **Test:** `e01_prep` · **PASS ✅**

Exact reproduction of the published cohort counts (11,868 samples / 866 resistant).

## Metrics

| Metric | Value | Expected | ±Tol | Status |
|---|---|---|---|---|
| `n_samples` | 11868 | 11868 | 0 | ✅ |
| `n_resistant` | 866 | 866 | 0 | ✅ |
| `cryptic_n` | 10704 | 10704 | 0 | ✅ |
| `cryptic_R` | 83 | 83 | 0 | ✅ |
| `nicd_n` | 1164 | 1164 | 0 | ✅ |
| `nicd_R` | 783 | 783 | 0 | ✅ |
| `unique_subjid` | 11611 | 11611 | 0 | ✅ |
| `mutation_rows` | 25554 | 25554 | 0 | ✅ |
| `minor_allele_rows` | 1115 | 1115 | 0 | ✅ |
| `null_rows` | 11 | 11 | 0 | ✅ |
| `null_rows_samples` | 4 | 4 | 0 | ✅ |

| Metric | Value |
|---|---|
| `median_FRSf_minor` | 0.3750 |

## Method
AGENTS.md §5: CRyPTIC QUALITY=HIGH + all NICD; dedup prefer R then MIC; inner-join GENOMES; mutations GENE@MUT with MINOR_MUTATION when minor, FRSf; drop synonymous/null.

## Novelty claim
Baseline reproduction (not a novelty claim).

## Caveats / negative results
median_FRSf_minor (~0.375) is approximate, not an exact acceptance row. Repo DST files hold 22,347 R/S IDs; the QUALITY=HIGH + genome join yields 11,868.

## Reproduce
```bash
python eval/make_manifest.py && python eval/e01_prep.py
```

## Provenance
- **Date (UTC):** 2026-11-18T19:04:30.131974+00:00
- **Git commit:** uncommitted
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `counts.csv`
