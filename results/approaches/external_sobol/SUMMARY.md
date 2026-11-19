# External — Sobol validation

> **Type:** validation · **Approach:** `external_sobol` · **0/1 tests pass** ❌

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E Sobol Inventory | — | n_accessions=744 | [report](e_sobol_inventory/REPORT.md) |

## E Sobol Inventory — DESCRIPTIVE •

744 South African run accessions inventoried. The 423 BDQ-pDST subset and exact pDST-R count are NOT verifiable locally (paper's '89' is genotypic, not pDST). Strongest external anchor (MGIT 1.0, assay-matched to NICD) once genotypes acquired.

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

_Caveats: BLOCKED: needs ENA raw reads -> variant calling, plus supplement pDST counts. No accuracy fabricated. MGIT 1.0 matches NICD training assay (favourable for transfer)._

