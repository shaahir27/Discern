# Runbook: external validation of Discern on Pakistani WGS data (for an AI coding agent)

Goal: score the **frozen** Discern model on two public Pakistani datasets (lineage-3-heavy) and report honest,
per-dataset numbers. **Never train on these data.** Read "Ground rules" first.

Tags used below: **[VERIFIED]** = checked by the person who wrote this runbook (against the repo, the papers, or by running code).
**[CHECK]** = could not be verified; confirm on your machine before relying on it. If a [CHECK] fails, stop and report.

---------------------------------------------------------------------------------------------------

## 0. Ground rules (stop and report if any would be violated)

1. Do not modify `models/lr_model/bdq_lr.joblib`, `models/lr_model/catalogue_bdq.csv`, or anything under `models/`.
2. Run `make verify` before you start and after every phase that touches `src/`. It must print 14/14. [VERIFIED: it passed 14/14 on the archive]
3. External data are **test data only**. No fitting, no threshold tuning, no vocabulary changes.
4. Write the pre-registration (Phase 7) **before** you score anything (Phase 8). Commit it. Do not edit it afterwards.
5. Keep every dataset separate. Never pool Saeed and Rashid. Never pool them with CRyPTIC/NICD/KM.
6. Do not call `write_result()` (it registers into `results/index.json`). Write plain CSV/JSON into `results/external_pending/`.
7. Log every command and tool version to `results/external_pending/LOG.md`. Save a manifest of every file downloaded with its SHA-256.
8. Make no India claim anywhere. These are Pakistani isolates.

## 1. What is verified about the two datasets

| | Saeed 2022 (Aga Khan Univ., Karachi) | Rashid 2025 (National TB Reference Lab, Islamabad) |
|---|---|---|
| Reads | SRA BioProject **PRJNA717333** [VERIFIED in paper] | "BioProject no. 1195559" [VERIFIED wording]; probably **PRJNA1195559** [CHECK] |
| Licence | CC BY 4.0 [VERIFIED] | CC BY 4.0 [VERIFIED] |
| Isolates | 48: 10 BDQ-R, 8 intermediate (MIC 0.25), 20 BDQ-S MDR, 10 pan-susceptible [VERIFIED] | 100: 50 BDQ-R, 50 BDQ-S controls chosen from RIF-resistant [VERIFIED] |
| Phenotype | Broth microdilution MIC. S <= 0.12, I = 0.25, R >= 0.5 ug/ml [VERIFIED] | MGIT 960, **no MIC**, BDQ concentration not stated [VERIFIED] |
| Lineage (authors) | 42 of 48 lineage 3; 3 L4; 2 L2; 1 L1 [VERIFIED] | L3 73, L4 15, L2 7, L1 5 [VERIFIED] |
| Platform | MiSeq 2x250; coverage 9.9-230x; 10 isolates re-sequenced on MiniSeq [VERIFIED] | HiSeq2000 [VERIFIED] |
| Authors' calling | BWA + GATK (not built for low-frequency variants) [VERIFIED] | BWA-MEM + bcftools; heteroresistance window 5-95% [VERIFIED] |
| Per-isolate labels | Additional File 1 (XLSX) = phenotypes; File 3 (XLSX) = genotypes/coverage [VERIFIED to exist; contents NOT opened] | Tables 3 and 4 in the paper list variants + BDQ status for variant-carrying isolates only [VERIFIED]. IDs **BDQ_FR_004 and BDQ_FR_030 are R in Table 3 but S in Table 4** [VERIFIED]. Full R/S list for all 100: [CHECK] SRA/BioSample metadata or ask the authors |
| Design caveats | Convenience sample; BDQ-naive patients | Hand-picked 50/50 case-control (no PPV!); 4 patients had prior BDQ; no MIC |

Label rules for scoring (frozen now):
* Saeed: R = MIC >= 0.5. Intermediate (0.25) = **S** for the primary analysis (same cutoff as CRyPTIC, where R means MIC > 0.25) and reported separately as a sensitivity analysis.
* Rashid: the MGIT R/S label exactly as published.

## 2. Laptop setup

Needs: 150 GB free disk (estimate), 16 GB RAM (estimate), internet. Linux or macOS. **Windows: use WSL2 Ubuntu and keep all files inside the WSL filesystem.**

```bash
# 2.1 conda (miniforge)
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-*.sh -b -p "$HOME/miniforge3" && source "$HOME/miniforge3/etc/profile.d/conda.sh"

# 2.2 TB-Profiler + SRA tools  [VERIFIED install route: bioconda]
mamba create -n tbp -c conda-forge -c bioconda tb-profiler sra-tools jq -y
conda activate tbp
tb-profiler update_tbdb            # [VERIFIED in the TB-Profiler README]
tb-profiler --version | tee -a results/external_pending/LOG.md
```
Apple Silicon: if the solve fails, retry with `CONDA_SUBDIR=osx-64` (Rosetta) or use Docker. [CHECK]

```bash
# 2.3 Discern python env: follow the repo's own README/requirements, then
make verify        # must be 14/14
```

## 3. Fix indel naming FIRST (the single biggest known risk)

Problem [VERIFIED by running code on the training data]: the same DNA change has several valid names. Training uses the
left-most name (`192_ins_g`); HGVS/TB-Profiler uses the right-most (`c.198dupG`). Discern's `hgvs_to_garc` does not convert
between them, so known variants look "novel". Of 164 distinct plain Rv0678 indel names in training, **128 sit in repeats and
have at least one other valid spelling**. Confirmed pairs: `144dupC` = `141_ins_c` (168 carriers), `198dupG` = `192_ins_g` (99),
`139dupG` = `138_ins_g` (58). Training itself spells 2 events two ways (3 rows in total).

Files provided with this runbook (already tested against the training data): `src/discern/indel_canon.py`, `tests/test_indel_canon.py`.

Steps:
1. `indel_canon.py` is in `src/discern/indel_canon.py` and `test_indel_canon.py` is in `tests/test_indel_canon.py`.
2. `pytest tests/test_indel_canon.py` passes all checks.
3. Build the index from training names: all `Rv0678@...` entries in `models/lr_model/evidence_store.csv` (column `variant`), `IndelIndex(RV0678_CDS, names, counts)`.
4. Use `hgvs_indel_to_garc()` then `IndelIndex.to_training_name()` for every Rv0678 indel in the converter (Phase 6).
5. Other genes (atpE, pepQ, mmpL5, mmpS5) have many indel names too (pepQ 291, mmpL5 80, atpE 69 distinct in the raw table), but `indel_canon.py` embeds only the Rv0678 CDS.
   To extend: download H37Rv (GenBank AL123456.3) from NCBI, extract the CDS by locus tag (atpE Rv1305, pepQ Rv2535c, mmpL5 Rv0676c, mmpS5 Rv0677c), reverse-complement minus-strand genes [CHECK strand],
   then **only enable a gene if every plain `_del_` name in training for that gene passes `del_matches_reference`** (this validates both the CDS and the position convention). Otherwise leave that gene's indels unchanged and report how many were skipped.
6. Measure and report (do not assume) how many external indel calls change catalogue call after canonicalization. The catalogue's wildcard rules may already catch frameshifts, so the effect may be small.

## 4. Phase 0 - "practice exam" on CRyPTIC isolates with known answers

Purpose: prove that TB-Profiler -> converter reproduces the variant names Discern was trained on, **before** looking at external data.

1. Pick 60 isolates from the Fowler data (`data/GENOMES.csv.gz`, `data/MUTATIONS_BDQ.csv.gz`): 30 carrying a common variant (include `141_ins_c`, `192_ins_g`, `138_ins_g`, `c-11a`, `pepQ@G197R`, `Rv0678@E55D`), 30 with no variants. Fix the random seed.
2. Each has a per-sample VCF at: `<EBI CRyPTIC FTP root> + FTP_PATH + FTP_FILENAME_VCF` (columns exist in GENOMES [VERIFIED]; example `vcfs/jeffk-20230406/` + `00/03/25/08/32508/site.01.iso.1.subject.01-DR0682.lab_id.DR0682.seq_reps.1.vcf.gz`).
   The FTP root is believed to be `https://ftp.ebi.ac.uk/pub/databases/cryptic/release_june2022/` [CHECK: download one file; if 404, search the EBI CRyPTIC FTP listing for `vcfs/`].
3. Run: `tb-profiler profile --vcf <file>.vcf.gz -p <UNIQUEID> -d results/external_pending/phase0` [`--vcf` support VERIFIED in the LSHTM course; whether it yields allele frequencies for Clockwork VCFs is CHECK].
   If TB-Profiler cannot use the CRyPTIC VCFs, look for run accessions in the CRyPTIC reuse table / Zenodo 15680920 `GENOMES.parquet` (columns CHECK) and re-call from reads for 20 isolates instead. If neither works, STOP and report: the pipeline cannot be validated and external results must be labelled "unvalidated conversion".
4. Convert (Phase 6 code) and compare with the Clockwork calls for the same isolates, 5 genes, FRS >= 0.1, non-synonymous only.
5. Report: per-variant recall, per-isolate exact-set agreement, extra variants (precision), and the relation between TB-Profiler `freq` and training `FRSf` (training sets FRSf = 1.0 for major alleles, = FRS for minor alleles [VERIFIED in `prep.py`]).
6. **Freeze** the FRSf rule from this comparison (e.g. "FRSf = freq" or "FRSf = 1.0 if freq >= x"). Do not change it after Phase 8.
7. Suggested gate (the owner may change it before running): recall of training-variant occurrences >= 95% and no systematic indel loss after Phase 3. If it fails, STOP and report.

## 5. Get the external data

```bash
mkdir -p data_external/{saeed,rashid} results/external_pending
# list runs and file sizes (ENA Portal API; documented, not run by the author) [CHECK]
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA717333&result=read_run&fields=run_accession,sample_accession,sample_alias,library_layout,instrument_model,base_count,fastq_bytes,fastq_ftp&format=tsv" \
  > data_external/saeed/runs.tsv
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1195559&result=read_run&fields=run_accession,sample_accession,sample_alias,library_layout,instrument_model,base_count,fastq_bytes,fastq_ftp&format=tsv" \
  > data_external/rashid/runs.tsv
```
* Sum `fastq_bytes` and report total size per dataset **before downloading** (author's rough guess, unverified: Saeed 5-20 GB, Rashid 30-100 GB). If free disk < 1.3x total, process in batches and delete FASTQ after each isolate's JSON is validated.
* If `PRJNA1195559` returns nothing, resolve "1195559" through NCBI BioProject search and report the real accession.
* Saeed has 48 isolates but 10 were re-sequenced: expect more than 48 runs. Per isolate, merge runs only if they are the same library type; otherwise keep the run with higher coverage. Record the choice per isolate.
* Download with `wget -c` from the `fastq_ftp` column (prefix `ftp://` or `https://`), or `prefetch` + `fasterq-dump` from sra-tools. Verify each file is non-empty and paired files have equal read counts.

Labels (a **blocking** step):
* Saeed: parse Additional File 1 (XLSX, download from the article page) for BDQ MIC per isolate. Map isolate names (S1..S23, SS1..SS8 pattern in the paper) to `sample_alias`. [CHECK mapping]
* Rashid: look in BioSample attributes for phenotype; else use Tables 3-4 (incomplete, contradictory for 004/030); else email the corresponding author (drshaukat@cust.edu.pk) for the 100-row table.
* Write `data_external/<set>/labels.csv` with columns `isolate_id, run_ids, y, mic, label_source`.
* **STOP and report if labels cannot be assigned to at least 90% of isolates of a dataset.** Never guess a label.

## 6. Run TB-Profiler and build Discern's mutation table

```bash
# one isolate first; inspect; then batch
tb-profiler profile -1 R1.fastq.gz -2 R2.fastq.gz -p SAMPLE -t 4 -d results/external_pending/tbp --txt --csv
tb-profiler profile -h   # [CHECK] read the --af / --depth help and record the values used; Discern needs variants down to FRS 0.1
jq '.main_lineage, .sub_lineage, .qc' results/external_pending/tbp/results/SAMPLE.results.json   # [CHECK] key names
jq '.dr_variants[0], .other_variants[0]' results/external_pending/tbp/results/SAMPLE.results.json   # [CHECK] field names
```
Output goes to `<dir>/results/<prefix>.results.json` [VERIFIED from README/collate code]. Time the first isolate and extrapolate; run at most 2 isolates in parallel with `-t 4`.

Inspect one JSON that has an Rv0678 frameshift. Discern's web importer reads only the `change` field [VERIFIED in `importers.py`]. If frameshifts appear there as protein strings
(`p.Asp47fs`) the importer cannot use them; use the nucleotide field instead (likely `nucleotide_change`, [CHECK]).

Converter `eval/external_convert.py` (write it; ~80 lines). Per isolate, per variant in genes Rv0678, atpE, pepQ, mmpL5, mmpS5 (accept aliases: mmpR5 = Rv0678, Rv0678 etc., as in `BDQ_GENE_ALIAS`):
1. SNV / promoter: `hgvs_to_garc(change)` (existing).
2. Coding indel in Rv0678: `hgvs_indel_to_garc(cds, nucleotide_change)` then `IndelIndex.to_training_name()`; keep the `how` flag in a log column.
3. Anything unparseable: keep it, flag it, count it, **never silently drop**.
4. `FRSf` = frozen rule from Phase 0.
5. Output columns: `UNIQUEID, GENE, suffix, MUT (= GENE + "@" + suffix), FRSf`, plus the flags from `discern.grammar.classify(suffix)`:
   `IS_SYN, IS_NULL, is_stop, is_indel, is_prom, is_mis, pos` [VERIFIED: this is how `prep.load_mutations` and `km_cohort` build them].
6. Assert `classify()` succeeds for every suffix; list failures.
7. Isolates with no variants must still appear in the isolate table (they are "S" by default).

QC per isolate: median depth, lineage, number of variants, TB-Profiler version, tbdb version. Pre-registered exclusion: median depth < 20x (owner may change before running). Excluded isolates are counted and listed, not hidden.

## 7. Pre-registration (write and commit BEFORE Phase 8): `docs/PREREG_external_pakistan.md`

Contents (owner may edit thresholds now, never later):
* Datasets, label rules (Section 1), exclusion rule (depth), frozen FRSf rule, tool/DB versions.
* Primary outcomes, **per dataset, never pooled**: catalogue firm-call sensitivity and specificity (Clopper-Pearson 95% CI); fraction of U calls; model AUC with the frozen thresholds from `art["thresholds"]`.
* Secondary: U treated as S and as R; intermediate-as-R for Saeed; results by TB-Profiler lineage (counts shown; no stratum claim if n_R < 10).
* Pass/fail:
  * Specificity lower 95% bound < 0.90 in either dataset => "specificity does not transport" (reported as a failure).
  * Dataset with n_R < 20 => sensitivity reported descriptively only.
  * Phase 0 gate fails => no external claims.
* Not allowed: pooled PPV/NPV (Rashid is a 50/50 case-control), any India claim, the word "validated", any comparison that pools sources.
* Allowed: "Evaluated on N isolates from Pakistan (lineage-3 majority), assay X: sensitivity a (CI), specificity b (CI)".

## 8. Score with the frozen model

Copy `eval/ext_km_validate.py` to `eval/ext_pakistan_validate.py`. Keep everything from `piezo.ResistanceCatalogue` onward **unchanged**; replace only the loader (`load_km_cohort`) with yours. Core of the scoring (same as the KM script [VERIFIED]):

```python
rc = piezo.ResistanceCatalogue(str(models_dir / "catalogue_bdq.csv"))
calls = apply_catalogue(rc, ext[["UNIQUEID", "y"]], muts, frs_min=cfg["features"]["frs_min"])
tm = ternary_metrics(y, calls)
art = joblib.load(models_dir / "bdq_lr.joblib")
core, a = sample_tables(muts, ids, art["frs_min"], tuple(art["dbd"]))
X = design(core, a, art["vocab"], "F2").values
raw_p, P = score(art["final"], art["platt"], X)
```
Run once per dataset (and once more with intermediate-as-R for Saeed). Save to `results/external_pending/<dataset>/`:
`metrics.json`, `calls.csv` (isolate, y, lineage, catalogue call, model p, variants), `false_negatives.csv`, `false_positives.csv`, `u_zone.csv`.
Also report how many calls changed because of Phase 3 canonicalization (run with and without).

## 9. Lineage stand-in check (answers "is my mmpL5 haplotype a lineage proxy?")

From each isolate's mmpL5 calls build flags `I948V`, `T794I`, `D767N` (class V--, VI-, VIN as in the earlier analysis) and cross-tabulate with TB-Profiler `main_lineage`/`sub_lineage`. Report the table; **no conclusion beyond it**.
Caution [VERIFIED in a public TB-Profiler example]: an isolate called lineage 1.1.2 carried mmpL5 T794I + I948V (the VI- class), so a "VI- = lineage X" guess may be wrong.
The same check should later be repeated on CRyPTIC v3.4 `GENOMES.parquet` (has a `LINEAGE` column [VERIFIED in `km_cohort.py`]; Zenodo record 15680920).

## 10. Deliverables

`results/external_pending/` (all CSV/JSON above), `docs/PREREG_external_pakistan.md`, `LOG.md` (commands, versions, checksums, exclusions, every [CHECK] outcome),
a one-page summary with the allowed-claims sentence per dataset, and `make verify` output (14/14) after the last change.

## 11. If something goes wrong

| Symptom | Likely cause | Action |
|---|---|---|
| tb-profiler install fails (Apple Silicon) | no arm64 builds | `CONDA_SUBDIR=osx-64` or Docker |
| No `frequency/freq` below 1.0 in JSON | AF thresholds | read `tb-profiler profile -h`, lower the soft AF cutoff, re-run one isolate |
| Rv0678 frameshifts show `p.XxxNNfs` | `change` is protein-level | use the nucleotide field; assert every indel has a `c.` form |
| Known variant not in vocab after conversion | naming | check `IndelIndex` output `how` column; run `test_indel_canon.py` |
| Fewer isolates than papers state | merged/removed runs | reconcile run -> isolate table; report differences |
| Labels missing | supplement/metadata gap | Section 5 STOP rule |

## 12. Known unknowns (do not hide these in the final report)

* Contents of Saeed's supplementary XLSX files and Rashid's SRA metadata were not opened.
* TB-Profiler runtime, memory, exact JSON field names, and exact AF-threshold semantics were not tested.
* ENA API queries, EBI FTP root and the Rashid accession prefix were not run from the author's environment.
* Authors' pipelines (GATK / bcftools) differ from Clockwork used in training; call differences are expected and are what Phase 0 measures.
