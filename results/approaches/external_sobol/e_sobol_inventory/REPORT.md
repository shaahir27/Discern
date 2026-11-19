# E Sobol Inventory

> **Approach:** `external_sobol` · **Test:** `e_sobol_inventory` · **DESCRIPTIVE •**

744 South African run accessions inventoried. The 423 BDQ-pDST subset and exact pDST-R count are NOT verifiable locally (paper's '89' is genotypic, not pDST). Strongest external anchor (MGIT 1.0, assay-matched to NICD) once genotypes acquired.

## Metrics

| Metric | Value |
|---|---|
| `n_accessions` | 744 |
| `claimed_bdst_pdst_n` | 423 |

- **country:** `South Africa`
- **claimed_pdst_n_verified:** `False`
- **assay:** `MGIT 1.0 (per paper; assay-matched to NICD)`
- **exact_pdst_R_count:** `None`
- **paper_89_is_genotypic:** `True`
- **genotype_pipeline_blocked:** `True`

## Method
Sobol 2026 KZN run-accession inventory (data/sobol_runs.tsv). pDST subset count and exact resistant count require the paper supplement; raw reads require ENA PRJNA1257273.

## Novelty claim
External validation anchor (Sobol KZN), pending genotype acquisition.

## Caveats / negative results
BLOCKED: needs ENA raw reads -> variant calling, plus supplement pDST counts. No accuracy fabricated. MGIT 1.0 matches NICD training assay (favourable for transfer).

## Reproduce
```bash
python eval/ext_validate.py
```

## Provenance
- **Date (UTC):** 2026-11-19T04:32:30.807804+00:00
- **Git commit:** 39dad3a
- **Data manifest:** `40a4ac904529…`
- **Versions:** scikit-learn 1.9.1 · pandas 3.1.1 · numpy 2.5.5 · scipy 1.18.1 · catomatic 0.2.3 · piezo 0.9.2 · pyarrow 26.0.1

## Files
- `results.json` (source of truth)
- `accessions.csv`
