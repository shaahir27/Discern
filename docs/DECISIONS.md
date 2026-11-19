# Decisions & Pre-registration — Discern

Convention (AGENTS.md rule 12): every ambiguity is resolved to the simplest option and logged
here. Every novelty experiment is pre-registered BELOW **before** it is run. No model is tuned
to hit a number; all rungs and negative results are reported.

---

## D001 — Environment & reproducibility (locked 2026-11-18T18:38Z)
- Python 3.12.16 (conda env `discern`). Direct interpreter used:
  `C:/Users/shaah/miniconda3/envs/discern/python.exe` (`conda run` hit a plugin error; bypassed).
- Exact versions: scikit-learn 1.9.1, pandas 3.1.1, numpy 2.5.5, scipy 1.18.1, pyarrow 26.0.1,
  catomatic 0.2.3, piezo 0.9.2. These are pinned in requirements.txt for cloud→local invariance.
- **Feature count is 120** (111 vocabulary + 9 biology aggregates) all-data; 91–97 per fold.
  The spec's "about 135" was inaccurate prose; all docs/latency claims use the measured 120.
- pandas 3.1.1 / numpy 2.5.5 are newer than the spec floor (>=2.2 / >=1.26). E01–E04 reproduced
  EXACTLY on these versions, so no drift; versions are pinned to keep it that way.
- Known warning: sklearn 1.9.1 deprecates `LogisticRegression(penalty=...)` (removed in 1.10).
  Results unaffected now (l2 is default). Cleanup item: drop the explicit `penalty="l2"` arg
  before upgrading sklearn past 1.9. Do NOT change it mid-experiment (would alter the fit path).

## D002 — Cohort definition (fixed)
- Production/training cohort = **11,868 samples / 866 resistant** (CRyPTIC QUALITY=HIGH 10,704/83
  + NICD 1,164/783), per AGENTS.md §5. The repo DST files hold 22,347 R/S IDs (incl. LOW/MEDIUM
  quality and SEQTREAT/CLIRES sources); the QUALITY=HIGH + genome inner-join yields 11,868.
  We do NOT silently switch cohorts. External cohorts are additional, never merged into training.

## D003 — Enhanced-mode isolation (hard rule)
- All novelty changes (mechanism features, epistatic veto, shrinkage, conformal) live behind
  config flags in `config.yaml:enhanced`, **default OFF**. Golden mode (flags off) must stay
  byte-identical to E01–E04. A test asserts golden outputs unchanged when flags are off.
- "Enhanced" results are reported as separate model folders, never overwriting golden results.

## D004 — G0 status: COMPLETE, 14/14 PASS ✅
- E01–E13 all within tolerance (see docs/PROGRESS.md gate table). Identical to the cloud run →
  cloud→local invariance confirmed on pandas 3.1.1 / numpy 2.5.5 / sklearn 1.9.1.
- Data manifest: 138 files, manifest_sha256 40a4ac904529c836… (data/raw/MANIFEST.json).

## D005 — Method-fidelity corrections during G0 (NOT tuning)
- `vocab_from`: `core[mask]` selects columns on integer indices in pandas 3.x → use `.iloc[mask]`.
- E10 ECE: adopt standard quantile binning (sklearn calibration_curve convention). Hand-rolled
  edges gave 0.142 vs 0.101; standard gives 0.101 exact. Brier matched throughout → probabilities
  were correct; only the binning convention differed. This is method fidelity, not number-chasing.
- E11 "1:1 balanced" means 1:1 RESAMPLING (undersample majority), not class_weight. class_weight
  keeps all data (AUPRC 0.684); resampling drops data (AUPRC 0.659 ≈ spec 0.653). Corrected.
- E08 subset = 312 samples / 143 sets (spec "about 328/170"); acceptance is the AUC range
  (0.50–0.65) → 0.593 passes. Subset-count variance noted; negative result unchanged.
- None of these changed a model to hit a target; they corrected the experimental method to match
  the spec's intent. All rungs and the E08 negative result are reported as-is.

---

# PRE-REGISTRATION (timestamp 2026-11-18T18:38:14Z)

At this timestamp `results/` contains ONLY G0 baseline runs (e01–e04). No E14/E15/E16/E17/EXT
result file exists. Each experiment below is locked before execution. Primary metrics and the
decision rule are fixed now; we will report whatever comes out, including failures.

## P-E16 — Mechanism biology + epistatic veto (G2)
- **Hypothesis:** position/reading-frame-aware truncation features (Snobre 2024) and a
  high-confidence mmpL5/S5 epistatic veto improve correctness of Rv0678 loss-of-function calls.
- **Method:** build `models/lof_effects.csv` from Rv0678 CDS (NC_000962.3:778,990–779,487, +,
  GTG start, 165 aa; ~200 nt 3' flank for alternative-reading-frame continuation). Classify each
  observed stop/indel: first-stop region (DBD ≤99 / mid 100–160 / late >160), ARF-conserved
  (≥90% identity over dimerization region). Veto fires only on HIGH-confidence LoF (frameshift/
  stop, not single-nt homopolymer indel); NEVER on atpE-driven R.
- **Primary outputs (descriptive FIRST):** carriers per first-stop class with n_resistant and
  fraction — does Snobre's 24.8% susceptible-frameshift replicate in THIS lineage mix? Report
  even if it contradicts. THEN ablation: catalogue+LR with mechanism features OFF vs ON, veto
  OFF vs ON, on S1 and cross-source. Enumerate EVERY call flip with its true label.
- **Decision rule:** ship enhanced only if it improves pre-registered primary metrics under
  patient-grouped CV without degrading external validation; otherwise report and keep golden.
- **Possible negative outcome:** few truncation carriers → negligible metric move; or Snobre
  fraction does not replicate → report the tension. Either is a valid result.
- **Caveats:** mmpL5 LoF calls are noisy (repeat-rich, 84,815 rows); veto may remove true R —
  audit flips. Small n → correctness claim, not an accuracy revolution.

## P-E14 — Conformal risk-controlled certificates (G3)
- **Hypothesis:** LTT/Bentkus risk control gives distribution-free abstention certificates,
  re-weighted to user prevalence, replacing hand-picked NPV/PPV thresholds.
- **Method:** on inner OOF calibrated scores, grid (t_S, t_R); certify S-zone error ≤α (0.015)
  and R-zone error ≤β (0.20) at δ=0.05 via Bentkus p-values; class-constant prevalence weights
  (Tibshirani 2019). Lineage as covariate only.
- **Primary metrics:** empirical S-zone/R-zone error vs targets; certified bound vs achieved;
  U-zone size + resistant-in-U; across S1 / LOSO / cross-source. Report where certificates
  BREAK under lab shift (expected; a result, not a failure).
- **Honest constraint:** ~170 positives/fold → R-zone bound is wide (report Bentkus bound, not
  the 0.80 point target). U-zone likely widens vs baseline. Pre-registered as acceptable.

## P-E15 — Value-of-information lab queue (G4)
- **Hypothesis:** info-gain ranking resolves more catalogue-U samples per phenotyped variant set
  than count-only ranking.
- **Method:** score(set) = n_U_carriers × expected beta-binomial entropy reduction after one
  label. Retrospective label-reveal simulation: for K=1,2,5,10,20,50 reveal existing labels of
  top-K sets, rebuild catalogue decisions, measure %U resolved. Baselines: count-only, random.
- **Primary output:** %U-resolved vs K curves for all three rankings. Report if info-gain barely
  beats count-only (valid negative).
- **Caveat:** retrospective, assumes catalogued variants transfer across lineages; not prospective.

## P-EXT — External validation (G5, the headline)
- **Hypothesis:** the tool generalizes to non-overlapping external cohorts.
- **Cohorts (assert ZERO training-ID overlap before testing):**
  - KMVALIDATION2025: 792 isolates / 35 R (CRyPTIC v3.4.0). Assay/CC UNDOCUMENTED (MEDIUM) →
    stratify SA vs China vs Italy; **sensitivity analysis, NOT headline proof.**
  - Sobol 2026 KZN: 423 with BDQ pDST, MGIT 1.0 (assay-matched to NICD). Get EXACT pDST-R count
    from supplement (paper's "89" is GENOTYPIC, not pDST). Strongest external anchor.
  - Choi 2026 Korea: 49 isolates; breakpoint ≥0.5 ≠ WHO CC → reclassify at 1.0 / 0.25 first.
- **Primary metrics:** sens/spec per cohort + caveats; pooled only with assay stratification.
- **Hard limit:** NO Indian external resistant isolates exist in public data → NO Indian-
  sensitivity claim, ever. Indian audit = false alarms only.
- **Possible negative outcome:** external spec collapse (cf. Rojas-Fuertes 88%→47.5%). Report
  honestly; still the first external evaluation for a BDQ triage tool.

## P-E17 — Borderline-MIC label-noise audit (G6)
- **Pre-registered as a SENSITIVITY ANALYSIS, not a model change.** Samples with MIC within ±1
  two-fold dilution of CC (CRyPTIC plate 0.25; MGIT 1.0): n, n_resistant, share of positives;
  E04 metrics recomputed excluding them, side-by-side. Production model keeps ALL rows. A
  margin-weighted fit may be shown as a third column but is NEVER promoted from this audit alone.

## P-G7 — Same-data head-to-heads
- vs WHO v2 as-is on identical 11,868 (binary AND ternary at matched FRS conventions; note our
  ternary 0.735 does NOT yet beat WHO at FRS 0.25 = 0.769).
- vs Adlard under BOTH shuffled (their protocol) and patient-grouped (ours).
- vs Ghosh ONLY if their 632 isolates are reconstructable from Tables S9/S10; check training
  overlap first; if infeasible, say so. NEVER compare AUC across different cohorts.

---

## Claims guardrails (binding)
- ALLOW: "beat WHO v2 binary sensitivity on identical data"; "first phenotype-weighted tNGS
  panel coverage audit (to our knowledge)"; "validated on N non-overlapping external isolates
  [cohorts, assay caveats]"; "first predictive integration of truncation-position/ARF reasoning
  for BDQ (to our knowledge)"; "audited epistatic downgrade rule".
- FORBID: "first lineage-deconfounded" (Kulkarni/catomatic prior art); "first ML/minor-allele/
  probability model" (Ghosh/Adlard); "catalogues ignore context"; "works on Indian patients";
  "predicts new variants"; "doubled PPV" (realistic ~50→60 at 2%); any cross-dataset "beat";
  unqualified "first externally validated"; guaranteed NPV/PPV without bound+δ; rupee savings.
- Sensitivity ceiling ~84% (731/866 carry a known-gene variant; ~16% invisible). Never promise beyond.

## No-tuning commitment
We will not adjust models, thresholds, or feature sets to hit a target number after seeing
results. The ladder (M0 baseline → M1 mechanism → M2 veto → M3 shrinkage → M4 conformal) is
ablated rung-by-rung; the highest rung that improves the pre-registered primary metrics under
patient-grouped CV AND does not degrade external validation is shipped. If none beats baseline,
we report that and ship baseline + external validation + certificates. Every number in the UI,
README, and model card is read from results/*.json via results/index.json — never hand-typed.

---

## D006 — Per-model production artefacts (root `models/`)
Root `models/` holds the RUNNABLE production artefacts (metrics live in `results/`). Layout is
**per-model** so retraining never overwrites an earlier model:
`models/<model_name>/{<model>.joblib, metadata.json, catalogue_bdq.csv, catalogue_version.json}`.
`metadata.json` carries versions, manifest hash, training prevalence (866/11868), thresholds
`lo`/`hi`/`t_bin`, feature names, and hashes of cited results files. The product loads the chosen
folder via a metadata pointer; default production pick = golden LR (F2) unless a higher ladder rung
wins under the no-tuning rule. P-E16's `lof_effects.csv` lives at `models/lr_mechanism/lof_effects.csv`.

## D007 — Step A model-class ablation (exploratory, a priori configs, NOT tuned)
Same F2 features, S1 patient-grouped 5-fold, nested patient-grouped Platt + thresholds as E04.
Generalized runner `nested_oof_est` reproduces golden E04 exactly (LR-L2 AUC 0.908539 == E04).
Fixed a priori hyperparameters: LR-L1 C=0.3 liblinear; RandomForest 300 trees min_samples_leaf=5
max_features=sqrt; HistGradientBoosting 200 iters lr=0.1 15 leaves l2=1.0; linear SVM C=0.3
(decision_function → Platt). Results (results/models/<name>/e05b_model_cv, descriptive):

| model | AUC | AUPRC | binary sens/spec | sens@spec0.985 |
|---|---|---|---|---|
| rf_model | 0.9094 | 0.6818 | 0.8337/0.9548 | 0.6594 |
| lr_model (golden) | 0.9085 | 0.6836 | 0.8349/0.9481 | 0.6039 |
| lr_l1 | 0.9053 | 0.6753 | 0.8395/0.9455 | 0.5993 |
| hgb_model | 0.9052 | 0.6787 | 0.8372/0.9528 | 0.6259 |
| svm_model | 0.8953 | 0.6829 | 0.7979/0.9766 | 0.6732 |

**Decision:** no candidate beats golden LR by >0.005 AUC (RF +0.0008, within noise). KEEP golden
LR (F2) as production. SVM trades AUC for higher sens@spec0.985 (0.673) and specificity (0.977) —
noted as an operating-point alternative, not a winner. match != beat; no tuning.

## D008 — Concurrent-agent collision (resolved)
A Codex run built the verifier milestone (recomputing `eval/verify.py`, `tests/test_verify.py`
7 tests, `Makefile`) before Step A. Qoder accidentally overwrote `eval/verify.py` and restored it
from the git index (`git show :eval/verify.py`). Confirmed green: verify 14/14 acceptance pass,
pytest 7/7 (the earlier "7 errors" was a `--basetemp` permission artifact, not a test failure).
User confirmed Codex stopped; Qoder continues solo. Added `eval/run_all.py` (ordered pipeline
orchestrator; runs all or a named subset). Codex's verifier/Makefile/tests kept intact.

## D009 — G5 external validation: characterization DONE, genotype pipeline BLOCKED
Verified from local data (eval/ext_validate.py): KMVALIDATION2025 BDQ = 792 isolates / 35 R
(prevalence 0.0442), all MEDIUM quality, METHOD_3/METHOD_CC UNDOCUMENTED, ZERO training-ID overlap;
lineage L2=559/L4=198/mixed=20/L1=9/L3=6 (R: L2=22/L4=9/L1=4). Sobol = 744 SA run accessions
(423 pDST subset + exact pDST-R count NOT verifiable locally; paper's "89" is genotypic). Choi = 49
Korea accessions (breakpoint >=0.5 != WHO CC; reclassify first).
**BLOCKER:** external variant calls are NOT local. The CRyPTIC v3.4 parquets present are
phenotype/metadata only (DST_MEASUREMENTS, DST_SAMPLES, GENOMES, WGS_SAMPLES); the schema documents
a MUTATIONS table that was not downloaded. Sobol/Choi are accession lists only (need ENA raw reads
PRJNA1257273 / PRJNA1301050 + variant calling + supplement pDST counts). No documented fetch URL in
the repo. **No external sens/spec fabricated.** To finish G5: acquire CRyPTIC MUTATIONS (BDQ genes)
for the 792 KM isolates (and/or ENA reads for Sobol/Choi), then run catalogue+model. KM remains a
sensitivity analysis (assay undocumented), NOT headline proof. NO Indian external R exists ever.

## D010 — E18 expanded-training: NEGATIVE result (more data does NOT beat metrics)
Downloaded full CRyPTIC v3.4.0 MUTATIONS (Zenodo 15680920). Tested richer training (KM held out):
- HIGH+MEDIUM expanded (23,441 / 1,046 R, +182 new R): S1 AUC **0.838** (golden 0.909);
  NICD->CRyPTIC sens collapsed to 0.332 (golden 0.675); KM external AUC 0.864 (golden prod 0.868).
- HIGH-only expanded (16,923 / 201 R, CRyPTIC only — NICD is MEDIUM so excluded): S1 AUC **0.667**;
  KM external AUC 0.863.
**Conclusion:** the golden pooled AUC 0.909 is a COHORT-COMPOSITION effect — NICD (783 R, high
prevalence, clear resistance) dominates the signal. Adding CRyPTIC isolates (low prevalence, harder
borderline/invisible resistance) dilutes it; CRyPTIC-only resistance is intrinsically hard (AUC 0.667).
Naive data expansion does NOT improve metrics. KM external generalization (~0.86) is stable regardless.
Golden baseline + production UNTOUCHED. Genuine improvement levers: assay-matched resistant-rich data
(e.g. Sobol KZN MGIT 1.0, like NICD) for diversity, and mechanism biology for hard CRyPTIC cases —
not quantity. No tuning; negative result reported honestly. `models/lr_expanded/e18*`.

## Acceptance verifier

Use the index as the registry of saved reports. Recompute each check from the reported
observations, targets, and tolerances, then compare the outcome with the saved flags.
Show descriptive records separately and exclude them from the acceptance count. A run
with no indexed acceptance checks fails. Write a full-precision failure report when a
value exceeds its tolerance.

## D011 — E19 production ship decision (veto rule locked 2026-11-19T05:51:21Z)

> Adopt the epistatic veto as a Tier A production catalogue modifier ONLY IF on the
> KM external cohort: (a) spec_firm(veto) >= spec_firm(golden) [no external specificity
> degradation]; AND (b) every veto flip is R->U with ZERO true-R->S [never misses a
> resistant as susceptible]; AND (c) the S1 sens_firm reduction <= 0.005 [D007 noise
> margin; actual -0.00063 passes]. Otherwise KEEP the golden catalogue and report the
> veto as an audited descriptive safety rule. Model class stays golden L2 LR regardless.

### Recorded outcome

**Production uses golden L2 LR (F2) and the golden catalogue in `models/lr_model/`,
with the audited Tier A veto enabled in `config.yaml`.** Setting
`enhanced.epistatic_veto: false` selects golden mode. The enabled predictor uses the
existing E16 high-confidence annotations in `models/lr_mechanism/lof_effects.csv`
because live requests carry variant names and fractions. The modifier acts on
annotated variants and sends affected R calls for lab testing.

| Rung | Outcome |
|---|---|
| M0 golden | Retained. E21 found no candidate that cleared the joint D013 gate. Golden KM model AUC was 0.868183. |
| M1 mechanism | AUC gain 0.000249 fell below the S1 margin. Sensitivity at 98.5% specificity rose 0.012702. The S1 result resolved the ship decision. |
| M2 veto | Adopted under D011. KM firm specificity stayed 0.988136 and firm sensitivity stayed 0.791667. KM had zero flips, so the external safety condition passed without observing a downgrade. S1 firm sensitivity fell 0.000631, within the 0.005 margin. Of the 41 S1 R-to-U downgrades, 2 were truly resistant and 39 truly susceptible. |
| M3 shrinkage | Evidence-card analysis only. The saved audit covered 1,042 variants and recorded 48 grade flips. |
| M4 conformal | Report the certificates and observed lab-shift failures. At delta 0.05 the R-zone bound was 0.196307 against a 0.20 target. The best S-zone bound was 0.015516 against a 0.015 target. The full policy would raise unresolved calls from 630 to 11,617 and stays disabled. |

The KM cohort is a sensitivity analysis with 35 resistant isolates, medium-quality
phenotypes, and undocumented assay and critical concentration. The certificates
and observed transfer failures motivate per-lab prevalence reporting.
[KM veto audit](../results/approaches/ship_decision/e19_km_veto/results.json) ·
[Ladder summary](../results/approaches/ship_decision/e19_ladder_summary/results.json) ·
[E21 comparison](../results/approaches/head_to_head/e21_comparison/results.json)

Production SHA-256 values at E19

| File | SHA-256 |
|---|---|
| `models/lr_model/catalogue_bdq.csv` | `a514efec80a9a63aafb76032c669d5916ebd84db135683714154b6ce186a7e3d` |
| `src/discern/decision.py` | `83cc79bc0d7b3dffb3b3fd816c87c137c8317588b7bd2c1f441d4c84260a46c3` |
| `src/discern/mechanism.py` | `8157ed7aee13edee5cd698e0950b0c9122e5f8d1b883144025f44c1e30bf1730` |
| `config.yaml` | `4bf1717f2c90757725be2f2cbb5263efab8df63ec110e1105483affcc30cbfd2` |
| `models/lr_mechanism/lof_effects.csv` | `d8b4ea377199b28f2960d21be3577ed036b4793064cab4e0f1f90c9603629ac4` |

## D012 — E21 tuning grids and evaluation protocol (locked before fitting)

Compare five classifiers on the golden F2 features and saved S1 folds. In each outer
training fold, learn the vocabulary from variants with at least three training carriers.
Select by mean ROC AUC using `GridSearchCV` with patient-grouped `GroupKFold(3)`
and `n_jobs=-1`. The estimator seed is 1 where applicable.

| Model | Fixed settings | Grid |
|---|---|---|
| `lr_l2_tuned` | LogisticRegression, penalty=l2, max_iter=4000 | C: 0.03, 0.1, 0.3, 1.0, 3.0 |
| `lr_l1_tuned` | LogisticRegression, penalty=l1, solver=liblinear, max_iter=4000, random_state=1 | C: 0.03, 0.1, 0.3, 1.0, 3.0 |
| `rf_tuned` | RandomForestClassifier, random_state=1, n_jobs=-1 | n_estimators: 300, 600; min_samples_leaf: 1, 3, 5, 10; max_features: sqrt, 0.5 |
| `hgb_tuned` | HistGradientBoostingClassifier, random_state=1 | learning_rate: 0.05, 0.1; max_leaf_nodes: 15, 31; l2_regularization: 0, 1, 5; max_iter: 200, 400 |
| `svm_tuned` | SVC, probability=False, max_iter=4000, gamma=scale | C: 0.1, 0.3, 1, 3; kernel: linear, rbf |

Refit the selected setting on the same three inner training splits for held-out
raw scores. Fit Platt scaling (logistic regression with C=1e6) on those scores,
using logit probability or the SVM margin. Use the golden threshold procedure
on calibrated inner scores (Youden, NPV ≥0.985, PPV ≥0.80, minimum 20 rows).
Refit on the full outer training set and score the held-out fold. Pool S1 scores
for AUC, AUPRC and sensitivity at specificity 0.95, 0.97 and 0.985. Use each
fold's training thresholds for binary sensitivity/specificity and ternary coverage.

Repeat selection, calibration and refit on all golden training samples, then score
the fixed KM cohort after checking for zero isolate-ID overlap. Report KM AUC,
sensitivity/specificity at the training Youden threshold, and ROC sensitivity at
specificity 0.95. The latter describes the held-out ROC curve. KM has undocumented
assay and critical concentration and is a sensitivity analysis. All E21 metrics
are descriptive.

Execution note, added before inspecting held-out scores. The first launch
left liblinear's random state implicit. The comparison was rerun with the
registered seed explicitly applied to that solver. The grids and rule stayed fixed.

## D013 — E21 decision rule (locked before fitting)

Rank candidates independently by unrounded S1 AUC and KM AUC, alongside the golden
LR references from E04 and the production KM evaluation. A candidate qualifies if
its S1 AUC improvement exceeds the existing D007 margin of 0.005 and its unrounded
KM AUC is at least the golden KM AUC. Keep golden LR if no candidate qualifies.
Record every candidate meeting the rule for separate E19 review. Log all settings
and outcomes. Keep the production decision with E19.

## D014 — G8 research explorer

The browser and CLI retrieve saved results for the bundled held-out records. The
upload accepts either recorded IDs or complete archived variant profiles. Profile
matching uses the original rows and allele fractions and ignores row order.
The interface shows the catalogue call, saved score, and laboratory finding first.
Study provenance, carrier evidence, and prevalence estimates remain accessible on
the result page. The Docker image uses the archived data and excludes trained
inference artifacts.

The uncertainty report resamples patients from the existing S1 predictions and
computes Wilson intervals from the saved KM firm-call denominators. It preserves
all E01–E21 result records and uses the seed already registered in the configuration.

## D015 — Backend validation and predictive-value reporting

The explorer returns typed upload errors with locations when the CSV provides
enough information. Replay derives its reference prevalence from the other
folds of the saved S1 cohort and marks individual uncertainty and call traces as
unavailable. The legacy engine defers calls on catalogue processing errors.

The prevalence control uses the archived probabilities of each firm call across
the full cohort. The legacy engine uses E03 rates in golden mode and E16 rates
when the veto is enabled. This corrects PPV and NPV reporting while preserving
saved calls, scores, and evaluation records.
