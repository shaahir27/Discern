# External — Choi validation

> **Type:** validation · **Approach:** `external_choi` · **0/1 tests pass** ❌

## Tests in this approach

| Test | Status | Headline | Detail |
|---|---|---|---|
| E Choi Inventory | — | n_accessions=49 | [report](e_choi_inventory/REPORT.md) |

## E Choi Inventory — DESCRIPTIVE •

49 South Korean run accessions inventoried. Paper breakpoint >=0.5 mg/L differs from WHO CC (MGIT 1.0 / plate 0.25); phenotypes must be reclassified first.

| Metric | Value |
|---|---|
| `n_accessions` | 49 |

- **country:** `South Korea`
- **breakpoint_caveat:** `paper breakpoint >=0.5 != WHO CC; reclassify at 1.0 / 0.25 first`
- **genotype_pipeline_blocked:** `True`

_Caveats: BLOCKED: needs ENA PRJNA1301050 raw reads -> variant calling + MIC values for reclassification. Small n (49) -> wide intervals. No accuracy fabricated._

