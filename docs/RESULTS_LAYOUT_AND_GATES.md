# Handoff to Claude — BDQ-Triage novelty verification phase

Purpose: re-verify the baseline, then EMPIRICALLY run the new experiments so every novelty
claim is backed by a real number before we build the product. This is a verification +
experiment phase, NOT the full product build yet.

## Attach these files to Claude
1. `AGENTS.md` (original build spec — golden reproduction stays in force)
2. `FINAL_NOVELTY_APPROACH.md` (AUTHORITY for novelty; supersedes conflicting parts of #3)
3. `NOVELTY_UPGRADES.md` (detail for conformal certificates §1, VOI queue §2, LoF biology §3)
4. Claude's own `BDQ-Triage Final Technical Report.md` (the baseline it already produced)
5. `reference/` scripts (the golden implementation that produced E01–E13)

## Paste-ready prompt (copy everything in the block)

```
You are continuing BDQ-Triage. Read AGENTS.md, FINAL_NOVELTY_APPROACH.md, and
NOVELTY_UPGRADES.md fully before doing anything. Where FINAL_NOVELTY_APPROACH.md conflicts
with AGENTS.md or NOVELTY_UPGRADES.md, FINAL wins. AGENTS.md's golden reproduction and
claims policy still bind you.

This phase = verify + measure, NOT tune. Hard rules:
- NEVER tune a model to hit a number. Pre-register every experiment in docs/DECISIONS.md
  BEFORE running it. Report ALL results including failures and negative results.
- The baseline (golden mode) must stay BYTE-IDENTICAL. All new biology/estimator changes live
  in an "enhanced" mode behind a config flag (default OFF). Assert golden outputs unchanged.
- Stick to the spec's 11,868 samples / 866 resistant for reproducibility. The repo files
  contain 22,347 R/S IDs (broader, includes LOW quality); document the QUALITY=HIGH + genome
  filter that yields 11,868. Do not silently switch cohorts.
- "Beat" claims ONLY on identical or reconstructable data. Never compare our AUC to a
  different cohort's AUC.

Do these in order, stop and report at each gate:

G0 — Re-verify baseline. Re-run E01–E13 with the reference scripts. Confirm counts
     (11,868 / 866) and the anchors (WHO v2 U-as-S 0.497 / 0.679; catalogue CV firm 0.944 /
     sens 0.792 / spec 0.984; LR AUC 0.909). If anything drifted, write
     results/ACCEPTANCE_FAILURES.md and STOP.

G1 — Pre-register. Write docs/DECISIONS.md entries for E14 (conformal), E15 (VOI queue),
     E16 (LoF/mechanism ablation), E17 (MIC label-noise audit), and EXT (external validation),
     each with its primary metrics and the rule "ship the best pre-registered rung, report
     failures." Timestamp before any results file.

G2 — Mechanism biology (E16). Build models/lof_effects.csv from the Rv0678 CDS
     (NC_000962.3:778,990–779,487, + strand, GTG start, 165 aa; fetch ~200 nt 3' flank for
     alternative-reading-frame continuation). Classify every observed stop/indel: first-stop
     codon region (DBD <=99 / mid 100–160 / late >160), ARF-conserved (>=90% identity over
     dimerization region). FIRST write the descriptive table (carriers per class, n_resistant,
     fraction) — report whether Snobre's 24.8% susceptible-frameshift finding replicates in
     THIS lineage mix, even if it contradicts. THEN ablate: catalogue+LR with mechanism
     features OFF vs ON, and the epistatic veto OFF vs ON (HIGH-confidence mmpL5/S5 LoF only;
     NEVER veto atpE-driven R). Enumerate every call flip with its true label.

G3 — Estimator + certificates (E14). Implement beta-binomial shrinkage of per-variant
     evidence over mechanism classes (lineage as COVARIATE, not stratifier). Implement
     conformal LTT/Bentkus risk-controlled thresholds with prevalence weighting. Audit:
     empirical S-zone/R-zone error vs targets, certified bound vs achieved, U-zone size,
     across S1 / leave-one-site-out / cross-source. Report where certificates BREAK under
     lab shift — that is a result, not a failure.

G4 — VOI queue (E15). Rank unresolved variant sets by expected information gain
     (n_carriers x beta-binomial entropy reduction). Retrospective label-reveal simulation:
     for K=1,2,5,10,20,50 reveal existing labels of top-K sets, rebuild catalogue decisions,
     measure %U resolved. Compare vs count-only and random baselines. Plot curves.

G5 — External validation (EXT). Train on the 11,868; test on NON-OVERLAPPING external
     cohorts (assert zero training-ID overlap):
       * KMVALIDATION2025: 792 isolates / 35 R (CRyPTIC v3.4.0, Zenodo 15680920). ASSAY/CC
         UNDOCUMENTED (QUALITY=MEDIUM) -> stratify SA vs China vs Italy, flag as sensitivity
         analysis, NOT headline proof.
       * Sobol 2026 KZN: 423 with BDQ pDST, MGIT 1.0 (assay-matched to NICD). Get the EXACT
         pDST-resistant count from the supplement (the paper's "89" is GENOTYPIC, not pDST).
         Sequences PRJNA1257273.
       * Choi 2026 Korea: 49 isolates, PRJNA1301050; breakpoint >=0.5 is NOT WHO CC ->
         reclassify at 1.0 / 0.25 before comparing.
     Report sens/spec per cohort + caveats. There are NO Indian external resistant isolates;
     make NO Indian-sensitivity claim.

G6 — MIC audit (E17). Pre-registered sensitivity analysis: samples with MIC within +/-1
     two-fold dilution of CC (CRyPTIC plate 0.25; MGIT 1.0); metrics with/without them.
     Do NOT promote a margin-weighted model to production from this audit.

G7 — Same-data head-to-heads. vs WHO v2 as-is on the identical 11,868 (binary AND ternary at
     matched FRS conventions — note our ternary 0.735 does NOT yet beat WHO at FRS 0.25 =
     0.769). vs Adlard under BOTH their shuffled protocol and our patient-grouped one.
     Ghosh head-to-head ONLY if you can reconstruct their 632 from Tables S9/S10 accessions
     (repo PulmonomicsLab/bdqr-mtb-standalone); check training overlap first; if infeasible,
     say so and do NOT claim a cross-dataset win.

Return: results/*.json for every experiment, a summary table (metric | baseline | enhanced |
external | delta), the descriptive Snobre-replication table, every call-flip list, all
negative results, and any ACCEPTANCE_FAILURES. State plainly which novelty claims the data
SUPPORTS and which it does NOT. Do not build the UI/product yet.
```

## Data sources (all public)
| Need | Where |
|---|---|
| Baseline cohort | github.com/fowler-lab/tb-bdq-cat (tarball, per AGENTS.md §4) |
| KM external + lineage | Zenodo 15680920 (CRyPTIC v3.4.0): DST_MEASUREMENTS.parquet, GENOMES.parquet (LINEAGE/SUBLINEAGE). NOTE: we already downloaded these to `cryptic_v34/` locally — reuse if available |
| Sobol KZN | PMC12637405 supplement (pDST table); sequences PRJNA1257273 |
| Choi Korea | PMC13028822 Tables 1–4; PRJNA1301050 |
| Rv0678 CDS | NC_000962.3:778,990–779,487 (+, GTG, 165 aa); UniProt I6Y8F7 |
| Ghosh 632 | LSA 9(7):e202503539 Tables S9/S10; github.com/PulmonomicsLab/bdqr-mtb-standalone |
| WHO v2/v3 | github.com/GTB-tbsequencing/mutation-catalogue-2023; WHO 9789240124653 (v3) |
| tNGS panels (P2 only) | Deeplex/AmPORE manuals; Rosendal 2026 PMC12903939 supplement |

## Guardrails to enforce when results come back
- Golden mode byte-identical; enhanced behind a flag.
- No tuning-to-hit-numbers; pre-registration timestamps precede results.
- 11,868 cohort fixed; 22,347 discrepancy documented.
- KM = sensitivity analysis (assay undocumented); Sobol = strongest anchor; Choi reclassified.
- No Indian-sensitivity claim; no cross-dataset "beat"; no "doubled PPV" (realistic ~50→60).
- Snobre tension reported if it doesn't replicate; epistatic veto never on atpE; certificates'
  breakage under shift reported.
- Every UI/README number must trace to a results/ file (AGENTS.md rule 1).

## Results folder layout (BINDING — applies to every experiment)

Two sibling containers inside `results/`. Trained models (`type=model`) live in `models/`;
every other approach (baseline/analysis/data/method/validation/comparison) lives in
`approaches/`. Each model/approach folder carries a consolidated summary page
(`MODEL_REPORT.md` for models, `SUMMARY.md` for the rest) plus one subfolder per test.

```
results/
  INDEX.md                     # human master list, grouped by type (Models/Baselines/Analyses/...)
  index.json                   # machine map: approach -> container/test -> source file (UI reads this)
  ACCEPTANCE_FAILURES.md       # only if a tolerance is missed (rule 5)
  models/                      # ==== trained models only (type=model) ====
    lr_model/                  # MODEL_REPORT.md   (Tier B, LR F2)
      e04_s1_patient_cv/  e05_ablation/  e06_cross_source/  e07_loso/
      e11_imbalance_learning_curve/        # each: results.json, results.csv, REPORT.md
    structure_test/            # MODEL_REPORT.md
      e08_novel_variant/
    lr_mechanism/              # MODEL_REPORT.md   (G2)
      e16_lof_ablation/        # descriptive Snobre table + OFF/ON flips
    lr_shrinkage/              # MODEL_REPORT.md   (G3)
  approaches/                  # ==== everything else ====
    data_prep/                 # type=data        SUMMARY.md
      e01_prep/
    who_baseline/              # type=baseline    SUMMARY.md
      e02_who_as_is_frs0.1/  e02_who_as_is_frs0.75/
    catalogue/                 # type=baseline    SUMMARY.md   (Tier A)
      e03_s1_catalogue_cv/
    hybrid/                    # type=analysis    SUMMARY.md
      e09_hybrid/
    prevalence_correction/     # type=analysis    SUMMARY.md
      e10_prevalence_shift/
    priority_list/             # type=analysis    SUMMARY.md
      e12_variant_priority/
    predictive_values/         # type=analysis    SUMMARY.md
      e13_ppv_npv/
    conformal/                 # type=method      SUMMARY.md        (G3)
      e14_certificate_audit/
    voi_queue/                 # type=analysis    SUMMARY.md        (G4)
      e15_label_reveal_sim/
    external_km/  external_sobol/  external_choi/   # type=validation  SUMMARY.md  (G5)
      ext_km_792/  ext_sobol_423/  ext_choi_49/
    mic_audit/                 # type=analysis    SUMMARY.md        (G6)
      e17_borderline_mic/
    head_to_head/              # type=comparison  SUMMARY.md        (G7)
      g7_who_v2/  g7_adlard/  g7_ghosh/
```

Routing: `write_result(model, ...)` maps legacy names via `ALIAS`, then places the folder by
`type` — `model` → `models/<name>/`, otherwise → `approaches/<name>/` (see `_container()` in
`src/discern/results_io.py`). `migrate()` is idempotent and relocates any folder to its correct
container. Naming: model/approach folder = canonical snake_case name; test folder = experiment
ID + short description. The `type` drives INDEX grouping and the summary filename. Each test
folder contains exactly:
- `results.json` — machine-readable source of truth (observed, expected, tolerance, diff,
  pass flags, versions, data-manifest hash, seed, git commit, timestamp).
- `results.csv` — the table(s) (or named `*.csv` if several).
- `REPORT.md` — human-readable, shareable. Template:
  1. Title: `<approach> — <test>`; date, run ID, git commit, data hash, library versions.
  2. What was tested (1–2 plain lines).
  3. Method (brief + link to the spec section).
  4. Key results table: metric | value | expected | tolerance | pass.
  5. Interpretation (plain English: what it means; which novelty claim it SUPPORTS or NOT).
  6. Caveats / negative results (state failures plainly).
  7. How to reproduce (exact command).
  8. Source files (results.json, results.csv).

Compatibility with AGENTS.md rule 1 ("every number read from results/*.json"): the UI and
docs read numbers via `results/index.json`, which maps each metric to its nested
`results.json`. Keep `index.json` and `INDEX.md` regenerated on every run. No hand-typed
numbers anywhere; REPORT.md numbers are copied from results.json, never invented.

## After Claude returns
Paste its summary table + results here. I will verify the load-bearing numbers against the
raw files (as I did for KM 792/35), check the guardrails, and confirm which claims are
empirically supported. THEN we start the product build (AGENTS.md M0–M9 + the enhanced layer).