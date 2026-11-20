# Progress — Discern

## G0 — Baseline reproduction (AGENTS.md E01–E13): PASS ✅
Run 2026-11-18/19, env `discern` (Python 3.12.16; sklearn 1.9.1, pandas 3.1.1, numpy 2.5.5,
catomatic 0.2.3, piezo 0.9.2). Data manifest 40a4ac904529c836… (138 files).
**14/14 result files pass_all=True.** Identical to the cloud run → cloud→local invariance confirmed.

| Exp | Model / test | Key result | Expected | Pass |
|---|---|---|---|---|
| E01 | _data_prep/e01_prep | 11,868 / 866; all 11 counts exact | exact | ✅ |
| E02 | who_v2/e02_who_as_is_frs0.75 | sens_U_as_S 0.4965 | 0.497 ±0.005 | ✅ |
| E02 | who_v2/e02_who_as_is_frs0.1 | sens_U_as_S 0.6790 | 0.679 ±0.005 | ✅ |
| E03 | catalogue_tierA/e03_s1_catalogue_cv | firm 0.9443 / sens 0.7924 / spec 0.9835; 661 U / 206 R | 0.944/0.792/0.984 | ✅ |
| E04 | lr_f2_baseline/e04_s1_patient_cv | AUC 0.9085, AUPRC 0.6836, sens 0.8349 / spec 0.9481; ternary 0.7347/0.9889; n_feat 120 | 0.909/0.684/0.835/0.948 | ✅ |
| E05 | lr_f2_baseline/e05_ablation | F0 0.791/0.503, F1 0.905/0.672, F2 0.909/0.684, F2@0.75 0.799/0.548 | same ±0.015 | ✅ |
| E06 | lr_f2_baseline/e06_cross_source | C→N 0.910/0.844/0.963; N→C 0.828/0.675/0.953 | 0.91/0.844/0.963; 0.83/0.675/0.953 ±0.03 | ✅ |
| E07 | lr_f2_baseline/e07_loso | pooled sens 0.7875 / spec 0.9537; per-site spec(n≥300) 0.944–0.979 | 0.788/0.954 | ✅ |
| E08 | structure_novel_variant/e08_novel_variant | AUC region 0.545 / +structure 0.593 (NEGATIVE); n=312, 143 sets | AUC 0.50–0.65 | ✅ |
| E09 | hybrid_exploratory/e09_hybrid | AUC-in-U 0.744; resolves 173 S(17 R)/126 R(51 S); 362 remain(114 R) | 0.744; 173/126/362 ±15 | ✅ |
| E10 | prevalence_correction/e10_prevalence_shift | N→C ECE 0.288→0.013, Brier 0.102→0.017; C→N ECE 0.568→0.101, Brier 0.501→0.124 | same | ✅ |
| E11 | lr_f2_baseline/e11_imbalance_learning_curve | balanced(1:1 resample) 0.905/0.659/0.606; LC 0.896/0.905/0.905/0.907/0.909 | 0.904/0.653/0.592; LC ±0.01 | ✅ |
| E12 | priority_list/e12_variant_priority | 332 sets; top5/10/20/50 = 28/34/40/50%; largest c-11a 102 (8.8% R) | 332; 28/34/40/50 ±3 | ✅ |
| E13 | predictive_values/e13_ppv_npv | model PPV@1/2/5/10% = 14/25/47/65%; cat 33/50/72/85% | same | ✅ |

### Method-fidelity corrections made during G0 (NOT tuning — logged in DECISIONS.md)
1. `vocab_from` used `core[mask]` (column-select on integer indices in pandas 3.x) → fixed to `.iloc`.
2. E10 ECE: hand-rolled binning gave 0.142 vs 0.101; switched to standard quantile binning
   (sklearn calibration_curve convention) → 0.101 exact. Brier matched throughout (probabilities correct).
3. E11 "1:1 balanced": first tried class_weight (wrong intervention); spec means 1:1 RESAMPLING →
   corrected → 0.905/0.659/0.606 within tolerance.
4. E05 result-writer used attribute access on dicts → fixed to key access.

### Notes / caveats
- E08 subset is 312 samples / 143 sets (spec said "about 328 / 170"); the acceptance criterion is
  the AUC range (0.50–0.65), which passes. Subset-count variance logged; negative result holds.
- Feature count is **120** all-data (91–97 per fold), not the spec's "about 135".
- sklearn 1.9.1 `penalty=` deprecation warning is harmless; cleanup deferred (do not change mid-run).

## G1 — Pre-registration: DONE ✅
docs/DECISIONS.md pre-registers E14/E15/E16/E17/EXT/G7 (timestamp 2026-11-18T18:38Z, before any
novelty result). No-tuning commitment + claims guardrails locked.

## Next: G2 (mechanism biology) → G3 → G4 → G5 → G6/G7, then product build.

## Acceptance verifier session (2026-11-19)

Added an offline verifier for saved reports and tests using synthetic records. This run
checked the existing G0 result files. Step A and the subsequent experiment and product
milestones remain pending.

```powershell
make verify PYTHON=C:/Users/shaah/miniconda3/envs/discern/python.exe 'PYTEST_ARGS=-q --basetemp=.pytest_cache/verify-tmp'
```

```text
"C:/Users/shaah/miniconda3/envs/discern/python.exe" eval/verify.py
Research use only
Status | Test | Checks
PASS | lr_model/e04_s1_patient_cv | 9 within tolerance
PASS | lr_model/e05_ablation | 8 within tolerance
PASS | lr_model/e06_cross_source | 6 within tolerance
PASS | lr_model/e07_loso | 2 within tolerance
PASS | lr_model/e11_imbalance_learning_curve | 8 within tolerance
PASS | structure_test/e08_novel_variant | 1 within tolerance
PASS | catalogue/e03_s1_catalogue_cv | 9 within tolerance
PASS | data_prep/e01_prep | 11 within tolerance
PASS | hybrid/e09_hybrid | 7 within tolerance
PASS | predictive_values/e13_ppv_npv | 16 within tolerance
PASS | prevalence_correction/e10_prevalence_shift | 8 within tolerance
PASS | priority_list/e12_variant_priority | 5 within tolerance
PASS | who_baseline/e02_who_as_is_frs0.1 | 4 within tolerance
PASS | who_baseline/e02_who_as_is_frs0.75 | 4 within tolerance
14/14 acceptance tests passed.
"C:/Users/shaah/miniconda3/envs/discern/python.exe" -m pytest tests/ -q --basetemp=.pytest_cache/verify-tmp
.......                                                                  [100%]
7 passed in 0.96s
```

## Step A — Model-class ablation: DONE ✅ (2026-11-19, Qoder solo)
`eval/e05b_model_comparison.py` on identical F2 features / S1 patient-grouped nested protocol.
Generalized runner reproduces golden E04 exactly (LR-L2 AUC 0.908539). A priori configs (D007),
not tuned. RF best AUC 0.9094 (<0.005 over LR, noise); SVM best sens@spec0.985 0.673.
**Decision: KEEP golden LR (F2).** Results: results/models/{lr_l1,rf_model,hgb_model,svm_model}/
e05b_model_cv (descriptive) + approaches/head_to_head/e05b_comparison. Added `eval/run_all.py`
(ordered orchestrator). Codex verifier restored intact; verify 14/14 + pytest 7/7 green.

## Next: G2 (mechanism biology, P-E16) → G3 → G4 → G5 → G6/G7, then product build.

## Novelty gates G2–G7 (2026-11-19, Qoder solo; all descriptive, no tuning; 14/14 acceptance intact)
- **G2 / P-E16 (mechanism):** Snobre susceptible-frameshift fraction REPLICATES: 0.2504 (162/647) vs
  reported 0.248. Mechanism features ΔAUC +0.0002 (negligible), sens@spec0.985 +0.0127. Epistatic veto
  (R→U) on S1: 41 flips = 39 true-S removed (benefit) + 2 true-R (harm); cross-source NICD→CRyPTIC 39
  flips. lof_effects.csv saved. Ship decision deferred to G5. `models/lr_mechanism/`, e16_lof_*.
- **G3 / P-E14 (conformal) + M3 shrinkage:** R-zone certifies (Bentkus bound 0.196 ≤ β0.20, achieved
  0.143); S-zone NOT certifiable at α0.015/δ0.05 (best achievable bound 0.0155) → certified policy
  abstains heavily (U 630→11617). Certificates BREAK under lab shift (C→N S-error 0.67 vs bound 0.0097;
  N→C R-error 0.98 vs bound 0.179) — motivates per-lab prevalence slider. Shrinkage: 1042 variants,
  mean|Δ| 0.143, 48 grade flips. `approaches/conformal/`, `models/lr_shrinkage/`.
- **G4 / P-E15 (VOI queue):** info-gain ≈ count-only (K10 Δ+0.000, K50 Δ+0.011), both ≫ random.
  Top-50 sets resolve ~53% of 661 referrals; top-10 ~38%. Valid negative on info-gain superiority.
  `approaches/voi_queue/`.
- **G5 / P-EXT (external):** KM cohort VERIFIED 792/35 (prev 0.0442), MEDIUM/undocumented assay, ZERO
  training overlap, lineage L2=559/L4=198. **KM external validation DONE** (downloaded CRyPTIC v3.4.0
  MUTATIONS.parquet, Zenodo 15680920 CC-BY 4.0; eval/ext_km_validate.py): production catalogue
  sens_firm 0.792 / spec_firm 0.988 (matches internal CV); model AUC 0.868, sens@spec0.95 0.686;
  lineage L2 sens 0.64 (weak spot), L1/L4 sens 1.0. Model raw binary threshold transfers poorly
  (spec 0.77) → validates catalogue-as-final-call design. Sobol 744 SA / Choi 49 Korea accessions
  inventoried; their sens/spec still BLOCKED on raw reads (#13). `approaches/external_km|external_sobol|external_choi/`.
- **G6 / P-E17 (MIC audit):** 454 borderline samples (191 R = 22% of positives); excl AUC 0.9085→0.9045.
  Sensitivity analysis only; production keeps all rows. `approaches/mic_audit/`.
- **G7 / P-G7 (head-to-heads):** WHO v2 ternary FRS0.25 sens_firm 0.768 > our model 0.735 (we do NOT beat).
  Adlard shuffled 0.7955 vs patient-grouped 0.7924 (grouping stricter). Ghosh infeasible (data not local).
  `approaches/head_to_head/g7_who_v2|g7_adlard`.

## Next: product build (production artefacts + decision/evidence + API + UI + tests + Docker + docs).

## E18 — Expanded-training richness test: NEGATIVE (2026-11-19; D010)
Downloaded full CRyPTIC v3.4.0 MUTATIONS (1.0 GB, Zenodo 15680920). Richer training (KM held out):
HIGH+MEDIUM (23,441/1,046 R) S1 AUC 0.838; HIGH-only CRyPTIC (16,923/201 R) AUC 0.667; both BELOW
golden 0.909. KM external AUC stable ~0.86. The golden AUC is NICD-composition-driven; more CRyPTIC
data dilutes it. More data != better metrics. Golden + production untouched. `models/lr_expanded/`.

## E21 — Tuned model-class comparison (2026-11-19)

Ran the D012 grids on the saved patient-grouped S1 folds and the fixed KM cohort.
Random forest led S1 AUC (0.910345) and histogram gradient boosting led KM AUC
(0.889224). Tuned L1 LR gained 0.000035 S1 AUC and 0.000830 KM AUC over golden.
The D013 rule retains golden LR (F2). All pre-existing result JSON and production
files checked matched their pre-run hashes.
[Rankings and full-training settings](../results/approaches/head_to_head/e21_comparison/comparison.csv)
and [decision](../results/approaches/head_to_head/e21_comparison/REPORT.md).

Ran `make verify PYTHON=C:/Users/shaah/miniconda3/envs/discern/python.exe`
with `PYTEST_ADDOPTS` pointing to a fresh workspace-local temporary directory for
the Windows sandbox. Verification output for E21 and the acceptance totals follows.

```text
DESCRIPTIVE | hgb_tuned/e21_km_external | no acceptance criteria
DESCRIPTIVE | hgb_tuned/e21_tuned_cv | no acceptance criteria
DESCRIPTIVE | lr_l1_tuned/e21_km_external | no acceptance criteria
DESCRIPTIVE | lr_l1_tuned/e21_tuned_cv | no acceptance criteria
DESCRIPTIVE | lr_l2_tuned/e21_km_external | no acceptance criteria
DESCRIPTIVE | lr_l2_tuned/e21_tuned_cv | no acceptance criteria
DESCRIPTIVE | rf_tuned/e21_km_external | no acceptance criteria
DESCRIPTIVE | rf_tuned/e21_tuned_cv | no acceptance criteria
DESCRIPTIVE | svm_tuned/e21_km_external | no acceptance criteria
DESCRIPTIVE | svm_tuned/e21_tuned_cv | no acceptance criteria
DESCRIPTIVE | head_to_head/e21_comparison | no acceptance criteria
14/14 acceptance tests passed.
"C:/Users/shaah/miniconda3/envs/discern/python.exe" -m pytest tests/ -q
........                                                                 [100%]
8 passed in 4.45s
```

## E19 — Final production decision (2026-11-19)

D011 passed. Ship the golden L2 LR (F2) and catalogue in `models/lr_model/`
with the audited Tier A veto enabled. KM had zero flips, with firm sensitivity
0.791667 and firm specificity 0.988136 in both modes. The S1 firm-sensitivity
reduction was 0.000631. The negative M1 result, descriptive M3 results, and M4
certificate and lab-shift results are recorded in the ladder summary.
[D011](DECISIONS.md#d011--e19-production-ship-decision-veto-rule-locked-2026-11-19t055121z) ·
[KM veto audit](../results/approaches/ship_decision/e19_km_veto/results.json) ·
[Ladder summary](../results/approaches/ship_decision/e19_ladder_summary/results.json)

All 53 pre-existing result JSON and model artifact hashes matched. The enabled
predictor reproduced the 41 recorded S1 downgrades without changing the model
outputs. With the flag off, all 51 checked prediction responses matched their
pre-change hashes. The smoke check used recorded inputs with variants. An empty
variant list already failed in the golden predictor before the E19 change.

Ran `make verify PYTHON=C:/Users/shaah/miniconda3/envs/discern/python.exe`
with `PYTEST_ADDOPTS` set to a fresh workspace-local temporary directory.

```text
DESCRIPTIVE | ship_decision/e19_km_veto | no acceptance criteria
DESCRIPTIVE | ship_decision/e19_ladder_summary | no acceptance criteria
14/14 acceptance tests passed.
"C:/Users/shaah/miniconda3/envs/discern/python.exe" -m pytest tests/ -q
........                                                                 [100%]
8 passed in 7.64s
```

## G8 — Offline research explorer (2026-11-19)

The CSV-first workflow accepts complete recorded profiles and lists of recorded
sample IDs. It presents the saved call, model score, and lab finding. Study
results, the review queue, and variant evidence use the archived data.
[Backend review](BACKEND_REVIEW.md) · [Demo](DEMO_SCRIPT.md)

The empty-variant regression now passes. The 51 prior flag-off response hashes
and 55 original result and model artifact hashes matched. The live HTTP and CLI
checks agreed on the saved outputs. Chrome selected the example file through the
file input and submitted it successfully. All seven pages rendered at desktop,
tablet, and mobile widths with the research banner and no page-level overflow.
Docker served the homepage and health endpoint with container networking disabled.

```text
14/14 acceptance tests passed.
"C:/Users/shaah/miniconda3/envs/discern/python.exe" -m pytest tests/ -q
23 passed, 2 warnings in 12.72s
```

## Laptop interface review (2026-11-19)

Replaced the sidebar with a compact header and centered the upload and real
example side by side. The result leads with the recorded call, lab finding, and
saved score, followed by per-variant study counts. The [backend review](BACKEND_REVIEW.md)
lists the data gaps that limit a call explanation.

Chrome rendered all seven pages at 1366 × 768, 1280 × 720, and 390 × 844 without
page-level horizontal overflow. The primary upload and result panels fit within
both laptop viewports. The file picker and drag-and-drop accepted the bundled CSV,
the file picker rejected a wrong extension, and the upload returned the matching
saved score. All 55 baseline result and model artifact hashes matched.

```text
14/14 acceptance tests passed.
23 passed, 2 warnings in 13.41s
```

## Backend hardening and predictive-value correction (2026-11-19)

Verified structured upload errors, CSV locations, HTTP status handling, and
omission of submitted values from errors. Replay now includes training-fold
prevalence and explicit statuses for missing call traces and individual score
uncertainty. The legacy engine defers catalogue processing failures and validates
complete mutation notation at its request schema boundary. [Backend review](BACKEND_REVIEW.md)

PPV and NPV now use actual call rates across the full cohort. Direct reweighting
of archived cases agreed with the explorer. The legacy engine uses the rates for
its selected archived policy. After substituting the old PPV and NPV solely for
comparison, all 51 golden response hashes matched. The enabled path still
reproduced all 41 downgrades, and all 55 baseline archived result and model
artifact hashes matched. The updated server passed live upload and replay checks.

```text
14/14 acceptance tests passed.
31 passed, 5 warnings in 14.03s
```

## Evidence-card browser verification

Installed Chromium and Playwright in the ignored workspace cache. Verified all
seven pages and their research-use banners in Chromium. The example CSV upload
matched the archived call, model score, lab result, and observed changes. The
prevalence control updated the summary estimates while keeping the saved score.
At widths of 1440, 768, and 390 pixels, the source breakdown opened by keyboard
and the result page had no horizontal overflow. Reduced motion disabled the
card entrance animation. The walkthrough detected no page errors or external
requests.

```text
14/14 acceptance tests passed.
31 passed, 5 warnings in 5.33s
```

## 5 Novel Clinical & Genomic Approaches: Completed & Verified (2026-11-20)

Implemented the 5 high-impact, test-driven additions without altering golden mode or retraining baseline models:
1. **Sequencing-Assay & tNGS Panel Blind-Spot Auditor (`src/discern/panels.py`, `eval/panel_audit.py`, `GET /panels`)**:
   - Profiles WGS, Deeplex Myc-TB, Deeplex Myc-TB XL, ONT AmPORE-TB, and ICMR-NIRT MinION.
   - Detects missing pump genes (epistatic veto blind), unsequenced core genes (atpE/pepQ), and homopolymer indels on Nanopore or at low VAF (<0.25).
   - Generates 5-gene coverage comparison table and inline 5-gene coverage strip.
2. **Unsupervised Batch Prevalence Estimation (MLLS / Saerens EM) (`src/discern/label_shift.py`)**:
   - Solves local prevalence estimation for regional labs lacking pDST using unlabeled isolate batches (N >= 30).
   - Computes 100-rep bootstrap 95% confidence intervals and deduplicates serial patient isolates.
   - Renders interactive inline SVG PPV curve (0.5%–25% prevalence) with shaded 95% CI band in `/queue` and dynamically updates slider.
3. **Longitudinal Serial-Isolate Tracking & Clonal-Ribbon Timeline (`src/discern/longitudinal.py`, `GET /patient/{patient_id}`)**:
   - Supports optional `patient_id` and `collection_date` (ISO or Day offset) columns in uploaded CSVs.
   - Computes clonal trajectories, detecting de novo emergence, rising minority subclones, and polyclonal competition.
   - Renders interactive inline SVG Muller plot with timeline scrubber and accessible fallback table.
4. **Native Lab-Export Importers & Mechanism-Stratified Reflex Guidance (`src/discern/importers.py`, `result.html`)**:
   - Ingests TB-Profiler JSON, Deeplex Myc-TB CSV, and AmPORE-TB TSV/CSV formats directly with HGVS-to-GARC variant normalization.
   - Maps resistance drivers to structural domains and provides actionable mechanism-stratified clinical reflex guidance (CFZ cross-resistance, pDST reflex, target protection).
5. **Lightweight Batch-Recurrence & Contamination Alert (`src/discern/surveillance.py`, `queue.html`)**:
   - Identifies identical rare variants across distinct patients to trigger WGS surveillance escalation.
   - Differentiates biological transmission from technical homopolymer sequencing artifacts (e.g., Rv0678 nt 192–198 poly-G runs).

```text
14/14 acceptance tests passed.
61 passed in 12.14s (0 warnings, 0 failures)
```

## Upload and interface review (2026-11-20)

Verified the five additions through the shared upload handler. The checks cover native export validation, panel report loading, dated timelines, independent browser sessions, and batch prevalence rendering. The release check read the saved evaluation artifacts.

```text
14/14 acceptance tests passed.
75 passed in 16.87s
```

Chrome checks at desktop, tablet, and phone widths exercised CSV and JSON intake, the queue filters, selectable timeline dates, and the performance comparison. The run recorded zero external requests or browser console errors. The local server responded on port 8000 under Python 3.12. The first container build timed out downloading dependencies, and a retry used the cached base layer.
