# BDQ-Triage — Novelty Upgrade Specification (addendum to AGENTS.md)

> **SUPERSEDED IN PART by `FINAL_NOVELTY_APPROACH.md` (2026-11-18, later).** A hostile
> prior-art re-audit killed the "lineage deconfounding" and "gene expansion" angles
> (Kulkarni 2025 + catomatic already do lineage adjustment; WHO v2 BDQ = Rv0678/atpE/pepQ
> only). The conformal certificates (§1), VOI queue (§2), truncation/ARF biology (§3) and
> epistatic veto (§4) here remain valid. Where the two files conflict, FINAL wins.
> Correction: Vargas 2021 is **AAC**, not Nat Comms.

Research date: 2026-11-18. This file is an **addendum**, not a replacement: everything in
`AGENTS.md` (golden reproduction, acceptance numbers E01–E13, claims policy, working rules)
stays in force. The upgrades below are **new layers and new experiments (E14–E17)** that sit
on top of the reproduced baseline. Rule 1 of AGENTS.md still applies: no number in this file
is a result; every new experiment defines *what to measure*, never an expected value.

---

## 0. Verdict on the problem statement

The statement is worth doing and the spec is unusually rigorous (patient-grouped splits,
nested calibration, prevalence transport, a pre-registered negative result). But as written,
**every headline component is already published**:

| Component in current spec | Already published by | Consequence |
|---|---|---|
| Catalogue on CRyPTIC+NICD, minor alleles to FRS 0.1, mmpL5-LoF masking | **Adlard et al., Microbial Genomics 2025** — the paper behind `fowler-lab/tb-bdq-cat` itself (same 11,867/866 cohort) | Tier A is a *reproduction*. Say so proudly; never claim it as new. |
| ML on BDQ with minor-allele read fractions + output probabilities + explainability | **Ghosh et al. 2026, Life Science Alliance ("BDQR-MTB")**: MLP/RF/SVM, AO/DP features, 632 isolates, AUC 0.87, web tool | "First ML with minor alleles/probabilities for BDQ" is dead. |
| Logistic regression for BDQ; PPV/NPV vs prevalence reasoning | **Nimmo et al. 2024**: Firth logistic regression, simulated PPV/NPV across prevalences, ~69% genotypic sensitivity ceiling | LR choice and prevalence awareness are not novel per se. |
| Regression instead of univariate catalogue grading | **Kulkarni et al. 2025, Nat Comms** (21,589 drug–variant pairs) | Do not claim. |
| MIC prediction / "beyond binary" | CRyPTIC 2024 PLOS Comput Biol; Wang 2025 | Do not claim. |
| CFZ phenotype helping interpret Rv0678 | **Sonnenkalb et al. 2023, Lancet Microbe** (230 variants graded with BDQ+CFZ MICs + structure) | CFZ-as-evidence is published; a CFZ multi-task model is also a data-join burden (repo has BDQ only; NICD has no CFZ). |
| Conformal prediction for AMR | mBio 2025, patient-level *E. coli* AST conformal (468k isolates) | Conformal-for-AMR exists. Conformal **risk-controlled abstention for TB triage at user prevalence** does not (searched, not found). |

**Where the defensible novelty lives:** the *triage system layer* — (1) sample-level
risk-controlled deferral with certificates, (2) a value-of-information queue that tells the
lab which variants to phenotype next (catalogue growth engine), (3) reading-frame-aware
loss-of-function biology that contradicts the WHO "any truncation = resistant" expert rule.
None of these three was found published for BDQ or TB. Plus the evaluation rigor itself
(patient-grouped + LOSO + cross-source + novel-variant holdout) exceeds Adlard's shuffled
5-fold and Ghosh's mixed BioProjects — that is the "Strongest Engineering" story.

---

## 1. UPGRADE 1 — Risk-controlled triage (conformal risk control + prevalence weighting)

**One-line pitch:** the tool abstains with a *mathematical certificate*, not a hand-picked
threshold — and the certificate is re-issued at the user's lab prevalence.

### What changes
AGENTS.md §9 picks `lo`/`hi` by empirical NPV≥0.985 / PPV≥0.80 on inner out-of-fold scores.
Keep that as the **baseline policy**. Add a second, certified policy:

- **Method:** Learn-then-Test / Risk-Controlling Prediction Sets (Bates et al. 2021;
  Angelopoulos et al. 2021/2022) with a null "abstain" option — exactly the S/U/R design.
- **Risks controlled (monotone in thresholds):**
  - S-zone risk = fraction of truly-resistant among S-calls ≤ α (default α = 0.015, i.e. NPV ≥ 98.5%).
  - R-zone risk = fraction of truly-susceptible among R-calls ≤ β (default β = 0.20, i.e. PPV ≥ 80%).
- **Algorithm (≈100 lines numpy/scipy):** on the inner out-of-fold calibrated scores, grid
  candidate (t_S, t_R) pairs; for each, compute empirical risk on calibration rows; apply
  fixed-sequence testing with **Bentkus p-values** (scipy.stats.binom) to certify the pair at
  level δ (default 0.05); pick the least-conservative certified pair (smallest U zone).
- **Patient grouping:** calibration and test *patients* are exchangeable under GroupKFold, so
  split-conformal marginal validity holds; within-patient correlation inflates variance
  (effective n ≈ patients). State this in the model card; do not claim i.i.d. guarantees.
- **Prevalence weighting:** under pure prior shift the conformal weights are class-constant
  (π_user/π_train for R rows, (1−π_user)/(1−π_train) for S rows) — Tibshirani et al. 2019
  weighted conformal. This certifies NPV/PPV **at the user's prevalence**, which the current
  §10 reweighting only estimates. If π_user is estimated rather than known, label the
  guarantee "approximate (Barber et al. 2022 coverage-gap bound applies)".
- **Honesty requirement:** with ~170 positives per outer fold, the R-zone certificate will be
  *wide* (Bentkus bound at empirical PPV 0.80 certifies roughly 0.69–0.75, not 0.80). Report
  the bound, never the point target, as the guarantee. The U zone will likely **widen** vs
  the baseline policy. That is the honest finding, not a failure — pre-register it.

### New experiment E14 — certificate audit
Measure, per policy (baseline vs certified), on S1 pooled OOF, on S2 (LOSO), on S3 (cross-source):
- empirical S-zone error rate and R-zone error rate vs their targets;
- certified bound vs achieved error (does the certificate hold out-of-site? out-of-source?);
- U-zone size and resistant-in-U count for both policies;
- coverage fluctuation across folds (report min/max, not just mean).
Write `results/e14_conformal.json` + CSV. **No acceptance tolerance** — this is new science;
report what happens, including certificate violations under site/source shift (expected and
interesting: exchangeability breaks across labs — that is a *result*, and it motivates the
per-lab prevalence slider).

### API/UI surface
- `/predict` response adds `policy` ("baseline" | "certified"), `certificate`
  (`alpha`, `beta`, `delta`, `bound_s`, `bound_r`, `prevalence_weighted`: true/false).
- Result page: a small "why we abstained" card showing the certified bound in plain language
  ("if we call this Susceptible, at most ~1.5 in 100 such calls would be wrong — certified on
  N held-out patients").
- `/performance`: E14 table next to E04.

### Claims
- ALLOWED: "first risk-controlled triage application to bedaquiline resistance (to our
  knowledge)", "distribution-free abstention certificates, re-weighted to your lab's
  prevalence", "certificates audited across sites and sources; they degrade under lab shift —
  we show by how much".
- FORBIDDEN: "guaranteed 98.5% NPV" without the bound/δ/conditions; "conformal prediction
  invented by us"; any guarantee claim on the R zone beyond the Bentkus bound.
- Cite in model card: Bates 2021 (arXiv:2101.02703), Angelopoulos LTT (arXiv:2110.01052),
  Tibshirani 2019 (arXiv:1904.06019), Barber 2022 (arXiv:2202.13415), Mozannar & Sontag 2020
  (learning-to-defer), and the mBio 2025 E. coli conformal paper as prior art in AMR.

---

## 2. UPGRADE 2 — Value-of-information lab queue (the catalogue grows itself)

**One-line pitch:** don't just send samples to the lab — tell the lab **which experiment
shrinks the unknown zone fastest**, and prove it retrospectively.

### What changes
E12 already ranks unresolved variant sets by carrier count. Upgrade the ranking to
**expected information gain (EPIG-style, Smith et al. 2023)**:

```
For each unresolved variant set s among catalogue-U samples:
  n_s   = number of U samples whose unresolved variants are exactly s
  prior = BetaBinomial posterior of P(resistant | s) from existing labelled carriers
          (pooled across the full dataset; a_s, b_s pseudo-counts, background 0.1 as in §7)
  gain_s = expected entropy reduction of that posterior after ONE new phenotype label
           (average over the two possible outcomes, weighted by their posterior predictive probs)
  score_s = n_s × gain_s
```
Pure scipy; no new data. Tie-break by n_s. Keep count-only and random rankings as baselines.

### New experiment E15 — retrospective label-reveal simulation (honest by construction)
1. Take the 661 catalogue-U samples (E03). Hide their labels.
2. For K = 1, 2, 5, 10, 20, 50: reveal the *existing* labels of all samples in the top-K
   variant sets (as if the lab had phenotyped them), rebuild the catalogue decision for the
   remaining U samples with the same §7 settings, and measure the fraction of U samples that
   now resolve to firm R/S.
3. Repeat for count-only and random rankings; plot %U-resolved vs K for all three.
4. Report the curve, not a single cherry-picked K. Write `results/e15_voi.json` + CSV + plot.
Limitation to state: the simulation assumes catalogued variants transfer across lineages;
it is retrospective, not a prospective trial.

### Product surface
- `/priority` page: rank column switches between "samples affected" (E12) and "information
  gain" (E15); each set gets an evidence card: carriers, current posterior, expected gain,
  "phenotyping this set would resolve ~N queued samples".
- `/queue`: unchanged sample-level ranking (Tier B probability), but each U sample links to
  its variant set's priority card — the sample queue and the experiment queue are connected.
- Demo line: "This is a to-do list for the wet lab, ordered by how much each experiment
  teaches the catalogue. The tool is designed to make itself less uncertain over time."

### Claims
- ALLOWED: "value-of-information prioritization of phenotyping experiments for resistance
  catalogue growth (not found published for TB)", "retrospective simulation: top-K sets
  resolve X% of the queue" (X read from results/e15_voi.json).
- FORBIDDEN: "prospective validation", "the lab saved money/time" (no cost claims, §17).
- Cite: Smith et al. 2023 EPIG (arXiv:2304.08151); lab-in-the-loop precedents exist in drug
  discovery (bioRxiv 2025; Serrano 2026) — cite as inspiration, claim only the catalogue-
  curation application.

---

## 3. UPGRADE 3 — Reading-frame-aware loss-of-function interpretation (the biology engine)

**One-line pitch:** a frameshift is not a verdict. **24.8% of mmpR5-frameshift isolates are
bedaquiline-SUSCEPTIBLE** (Snobre et al. 2024, AAC — verified). The WHO expert rule "any
loss-of-function = resistant" is wrong for a quarter of truncations, and our tool reads
*where* the protein breaks, not just *that* it breaks.

### Published basis (verified against the paper)
- Stop codons inside the DNA-binding domain (codons 34–99, esp. α2–α3 helix 36–62) → resistant.
- Late stops in the second dimerization domain (codons ~100–160) or at the C-terminus
  (Rv0678 is 165 aa) → protein retains near-WT structure/function → often susceptible.
- Alternative reading frames (ARFs) after a frameshift can restore a conserved sequence
  (40.7% of frameshifts in Snobre's analysis); late-stop in the original ORF occurred for 32%.
- No Rv0678–DNA experimental complex exists in the PDB (4NB5 is protein-only, 1.64 Å,
  Radhakrishnan 2014). The repo's `Rv0678_WT.pdb` DNA chains C/D are a **model** — the model
  card must say "distance-to-DNA features computed against a modeled complex" (see §7 below).

### Implementation (policy layer, golden reproduction untouched)
New module `src/bdq_triage/lof_policy.py`, applied in `decision.py` **after** Tier A aggregation,
only in "enhanced" mode (config flag `lof_policy: enhanced|off`, default off for reproduction):

1. **Build-time asset** (`eval/build_lof_table.py`, pure Python, no torch): fetch Rv0678 CDS
   from NC_000962.3 (coordinates in the repo/WHO catalogue; record source+hash in manifest).
   For every observed indel/stop variant in the dataset, translate the mutated ORF and the
   alternative frames; record: first stop codon position (original frame), whether an ARF
   yields a ≥90%-identity protein over the dimerization region, affected domain
   (DBD 34–99 / dimerization 100–160 / C-term >160). Store `models/lof_effects.csv`.
2. **Decision rules (enhanced mode):**
   - Truncation with first stop ≤ 99 (DBD destroyed) → keep catalogue R evidence (confident LoF).
   - First stop > 160, or ARF restores conserved sequence → **do not auto-R from this variant**;
     evidence card says "late truncation — protein likely functional (Snobre 2024)"; the
     sample's call falls back to its other variants, else U with warning `LATE_TRUNCATION_SUSPECT`.
   - Frameshift with first stop in 100–160 and no conserved ARF → U with warning
     `TRUNCATION_UNCERTAIN_POSITION` (this is the genuinely uncertain middle).
3. **Warnings added:** `LATE_TRUNCATION_SUSPECT`, `TRUNCATION_UNCERTAIN_POSITION`.
4. **Evidence cards** for truncations show: domain diagram position, first-stop codon,
   ARF-conservation flag, and the Snobre citation.

### New experiment E16 — LoF-policy ablation (descriptive first, then ablation)
- Descriptive table: all samples carrying Rv0678 stops/indels, stratified by first-stop region
  (≤99 / 100–160 / >160 / ARF-conserved), with n, n_resistant, fraction. **Report whatever it
  shows** — if our cohort contradicts Snobre (different lineage mix), that is a published-
  tension finding and must be reported, not hidden.
- Ablation: catalogue + Tier B metrics with policy off vs enhanced, on S1 and S3. Expect the
  change to concentrate in a small number of samples (Adlard found 54 mmpL5-masked samples;
  late-truncation carriers will be similarly few). Report per-sample call flips with labels.
- Write `results/e16_lof_policy.json` + CSV. No tolerance; no tuning.

### Claims
- ALLOWED: "position-aware interpretation of Rv0678 truncations, following Snobre et al. 2024",
  "we do not treat every frameshift as resistance — published evidence shows ~1 in 4
  frameshift isolates remain susceptible".
- FORBIDDEN: "we discovered frameshift tolerance" (Snobre did); "our rule fixes the WHO
  catalogue" (show the ablation, let the jury conclude).

---

## 4. UPGRADE 4 — Epistatic veto at the decision layer (small, grounded)

mmpL5/mmpS5 loss-of-function **nullifies** Rv0678-mediated resistance (pump broken → nothing
to derepress): Vargas et al. 2021 AAC (82 double-LoF isolates; "mmpR mutations are not
valid markers" when the pump is broken), Nimmo 2024 Genome Med, Adlard 2025 (54 samples).
The spec already computes `mmpL5_LoF`, `rv_x_L5LoF` and the `PUMP_DAMAGE_MODIFIER` warning,
but the final call can still be R from an Rv0678 variant while the pump is broken.

**Change (enhanced mode only):** if the catalogue R is driven *solely* by Rv0678 variants and
`mmpL5_LoF > 0` (or mmpS5 LoF), downgrade the call to **U** with `PUMP_DAMAGE_MODIFIER` and an
evidence card explaining the epistasis ("the pump this mutation would unleash is itself broken").
atpE-driven R calls are NOT downgraded (target-site resistance is pump-independent).
Count and report every flipped sample with its label in E16's CSV (same experiment file,
separate section). This is a decision-layer rule, not a model change; golden mode untouched.

---

## 5. UPGRADE 5 — Borderline-MIC label audit (pre-registered sensitivity analysis)

BDQ MGIT960 critical concentration = **1.0 µg/mL, WHO interim** (Kaniga 2020) — contested,
so labels near it are genuinely noisy. `METHOD_MIC` is in the DST tables; the repo also ships
`UKMYC_GROWTH_BDQ.csv.gz` (per-well growth) if deeper auditing is wanted.

**E17 (audit only, NOT a model change):**
- Table: samples with MIC within ±1 two-fold dilution of CC (0.5–2.0 µg/mL): n, n_resistant,
  share of all positives; metrics (E04 set) recomputed excluding them, reported side-by-side.
- Pre-register in `docs/DECISIONS.md` BEFORE running: "this is a label-noise sensitivity
  analysis; we will not select the better-looking number; the production model keeps all rows."
- Optional stretch (only if pre-registered): margin-weighted fit
  `w_i = min(1, |log2(MIC_i/CC)|)` reported as a third column, never promoted to production
  on the basis of this audit alone.
- Write `results/e17_mic_audit.json`. This converts a hidden weakness (noisy labels) into a
  credibility exhibit.

---

## 6. Optional stretch items (only if the core is green)

### 6a. ESM-2 zero-shot lookup for novel missense (P2, honest-null-ready)
- Build-time only: esm2_t6_8M on CPU, 165 masked forward passes → wt-marginal LLR for all
  ~3,135 single missense variants of Rv0678; ship as static `models/esm2_llr.tsv` (~100 KB).
  Runtime stays pure numpy (CSV join). Same for atpE/pepQ if time permits.
- Evaluate exactly like §11.1 (S4 novel-variant holdout, missense-only subset): does the LLR
  feature beat the structure features' AUC ~0.57? **Pre-commit to reporting either outcome.**
- Novelty caveat: a 2026 bioRxiv preprint already uses ESM-2 features for MTB resistance
  prediction (incl. BDQ). Claim only "zero-shot novel-variant triage lookup integrated into
  the catalogue-U evidence path", never "first protein-LM for MTB".
- Coverage caveat: most novel Rv0678 resistance alleles are indels/nonsense — the lookup
  cannot score them. Say so on the evidence card.

### 6b. WHO catalogue v3 freshness (timeliness, not a BDQ breakthrough)
- WHO 3rd edition released **10 Nov 2026** (>71,000 isolates, >75 countries). Verified: the
  expanded evidence headlines are **pretomanid and cycloserine, NOT bedaquiline**.
- Cheap win: download the v3 graded list, diff BDQ-gene gradings v2→v3, show a small
  "what changed" table on `/performance`, and state in the model card that the tool's
  catalogue is data-driven (Adlard-style) and can be rebuilt when new phenotypes arrive —
  which is exactly what Upgrade 2's queue feeds. Do not claim v3 improves BDQ.

---

## 7. Corrections to current docs (do these regardless of upgrades)

1. **Structure provenance:** `Rv0678_WT.pdb` chains C/D (DNA) are a **modeled complex** — no
   experimental Rv0678–DNA structure exists (4NB5 = protein-only, Radhakrishnan 2014 JBC).
   Model card + `docs/DECISIONS.md` must say distance-to-DNA is computed against a model.
   DBD range 34–99 is confirmed correct (winged HTH; key contacts R84, D90, R92).
2. **Invisible-resistance framing:** our ~16% (131/866 no Tier-1 variant, Adlard) vs Nimmo's
   ~29.5% no candidate-gene mutation are *different definitions*. Define both precisely in
   the glossary; never mix them.
3. **Minor-allele framing:** emergence of BDQ resistance during treatment is confirmed
   (Derendinger 2023 Lancet Microbe: 48% acquired on treatment; Ardizzoni 2025 AAC: acquired
   Rv0678 → unfavorable outcome OR≈132). BUT baseline minor alleles predicting later fixation
   is an **evidence gap** (minor variants fluctuate/transient; FRS does not correlate with MIC
   cross-sectionally per Adlard). Frame the FRS-0.1 feature as "detection sensitivity for
   emerging-resistance candidates worth monitoring", NEVER "predicts treatment failure".
4. **Lineage limitation:** cross-source performance drop (AUC 0.91→0.83) plausibly reflects
   lineage/founder differences (NICD L2/L4-enriched vs global CRyPTIC; Eswatini M146T outbreak
   shows founder markers masquerading as resistance; Kulkarni 2026 Nat Comms shows lineage
   SNPs improve BDQ MIC prediction). The repo has no lineage calls — state this as a known
   limitation and future work; do NOT invent lineage features from BDQ-gene variants alone.
5. **WHO v3 exists** (10 Nov 2026) — model card baseline section must acknowledge it.
6. **Prior-art citations** (model card "related work" section, new): Adlard 2025 (baseline),
   Ghosh 2026 (ML prior art), Nimmo 2024 (LR + prevalence prior art), Kulkarni 2025
   (regression grading), Sonnenkalb 2023 (CFZ + structure), Snobre 2024 (frameshifts),
   mBio 2025 (conformal AMR prior art). A related-work section that names its collisions is
   itself a credibility signal to a technical jury.

---

## 8. Rejected ideas (do not build; reasons recorded so nobody re-proposes them)

| Idea | Verdict | Reason |
|---|---|---|
| CFZ multi-task / auxiliary-label model | REJECT for hackathon | Sonnenkalb 2023 already used CFZ to grade Rv0678; repo has BDQ-only tables (CFZ join = full day, CRyPTIC EBI FTP, NICD missingness not-at-random); CFZ CC also interim; Rv0678-independent CFZ mechanisms inject noise. |
| Extra epistasis features (Rv0678×atpE, multi-hit Rv0678) | REJECT | atpE variants ~1–2% of isolates; double carriers single digits among 866 positives → noise fitting. Keep the one grounded mmpL5 veto (§4). A descriptive double-mutant count table + Fisher test (underpowered, labeled) is acceptable in E16. |
| Deep learning (CNN/transformer on variants) | REJECT | 866 positives; Ghosh 2026 MLP already occupies the slot; spec non-goal. |
| MIC regression / ordinal MIC model | REJECT | CRyPTIC 2024 + Wang 2025 published; scope creep. |
| "Heteroresistance early-warning predicts failure" claim | REJECT as claim | Evidence gap (§7.3). Keep minor-allele detection; frame as monitoring candidates. |
| Lineage de-confounding features | REJECT for now | No lineage data in repo; would require external SNP calls — future work, honest limitation. |
| Margin-weighted training as production change | REJECT | Reads as tuning-to-hit-numbers; allowed only as pre-registered E17 sensitivity column. |

---

## 9. Pitch narrative (5-minute demo spine, extends docs/DEMO_SCRIPT.md)

**Positioning sentence:** "BDQ-Triage is not another resistance predictor. It is a
risk-controlled triage system that knows what it doesn't know, tells the lab exactly which
experiment to run next, and reads the biology of each mutation — built for the BPaLM rollout
where bedaquiline DST is the bottleneck."

**Why India / why now (AI for Bharat + HealthTech tracks):**
- India rolled out the all-oral BPaLM regimen under NTEP from Sept 2024 (PIB, 9 Aug 2024) —
  bedaquiline is now frontline for DR-TB, and resistance to it is the regimen's failure mode.
- Of 27 supranational TB labs surveyed, only 9 performed BDQ DST (EID 2021); MGIT needs
  culture + ~2 weeks. Sequencing is faster than phenotyping — triage is the missing layer.
- Pooled BDQ resistance among DR-TB globally: 5.7% (95% CI 3.6–8.3; Hu 2025 meta-analysis,
  31 studies). India estimate ~7.1% (Sakhi 2026, single-source — label it). Use these as
  *slider presets labelled "literature estimate, replace with your lab's prevalence"*,
  never as Indian validation claims (§17 still forbids "works on Indian patients").

**Demo spine (extends the existing script):**
1. Clear R: `Rv0678@141_ins_c` → catalogue R, evidence card (existing).
2. **Biology engine:** a late-stop Rv0678 truncation → enhanced mode shows "protein likely
   functional (Snobre 2024)", call U instead of blind R; side-by-side with WHO-rule behavior.
   Then a pump-broken double mutant → epistatic veto card.
3. **Certainty engine:** unresolved `c-11a` case → U with the conformal certificate card
   ("at most ~X in 100 S-calls wrong, certified on N held-out patients"); slide prevalence
   7%→2% and watch the certificate re-weight and PPV fall (existing slider, now certified).
4. **Priority engine:** `/priority` sorted by information gain; show the E15 simulation curve
   — "phenotype these 10 variant sets and ~X% of the queue resolves". The tool hands the lab
   a ranked experiment list, not just a sample list.
5. **Honesty moat:** `/performance` — patient-grouped CV, LOSO, cross-source, the E08
   novel-variant NEGATIVE result (AUC ~0.57) shown proudly, E14 certificate audit including
   where certificates break under lab shift. "Every number on this page is regenerated by
   `make verify` from public data."

**Prize mapping:** Strongest Engineering → evaluation rigor + certificates + reproducibility;
Best User Impact → lab queue + BPaLM bottleneck; originality → the three engines, each with a
verified white-space claim and a related-work section that names its collisions.

### Claims-policy diff (extends AGENTS.md §17)
ADD to allowed: "first risk-controlled triage application to BDQ (to our knowledge)";
"certified abstention at your lab's prevalence (approximate when prevalence is estimated)";
"retrospective simulation of catalogue growth"; "position-aware truncation interpretation
following Snobre 2024".
ADD to forbidden: "first ML/minor-allele model for BDQ"; "predicts treatment failure from
minor alleles"; "guaranteed NPV/PPV" without bound+δ+conditions; "WHO v3 improves BDQ";
"frameshift tolerance discovered by us"; any Indian-sensitivity claim (unchanged).

---

## 10. Build order and effort (3-day on-site budget)

| Order | Item | Effort | Gate |
|---|---|---|---|
| 0 | AGENTS.md M0–M9 baseline green (unchanged) | per spec | `make verify` E01–E13 |
| 1 | §7 doc corrections + related-work section | 0.5 day | grep tests still green |
| 2 | Upgrade 4 epistatic veto + Upgrade 3 LoF policy + E16 | 1 day | descriptive table written; flips enumerated; golden mode byte-identical |
| 3 | Upgrade 1 conformal + E14 | 1 day | certificate audit written (violations reported, not hidden) |
| 4 | Upgrade 2 VOI queue + E15 | 0.5–1 day | simulation curves vs both baselines |
| 5 | E17 MIC audit (pre-registered) | 0.25 day | decision logged before run |
| 6 | Stretch: ESM-2 lookup / WHO v3 diff | only if 1–5 green | honest null pre-committed |

New tests to add (extend AGENTS.md §13): LoF-table build determinism; enhanced-vs-golden
mode isolation (golden outputs byte-identical with flag off); conformal certificate math
unit tests (Bentkus p-value vs scipy binom; weighted conformal reduces to unweighted at
π_user = π_train); VOI simulation never reads hidden labels before reveal (leakage test);
E17 pre-registration file exists before results file (timestamp check).

---

## 11. References (verified 2026-11-18)

- Adlard et al. 2025, Microbial Genomics 10.1099/mgen.0.001429 (baseline paper; bioRxiv 2025.01.30.635633)
- Ghosh et al. 2026, Life Science Alliance 9(7):e202503539 (BDQR-MTB; PMC13135273)
- Nimmo et al. 2024, Genome Med / PMC11072239 (BDQ genomics limits, Firth LR, prevalence PPV/NPV)
- Nimmo et al. 2024, PMC10877763 (historic Rv0678 reservoir; mmpL5 LoF negation)
- Kulkarni et al. 2025, Nat Comms (PMC11876447; regression re-grading of WHO catalogue)
- Kulkarni et al. 2026, Nat Comms (PMC13315309; lineage-aware CNN for BDQ MIC)
- Sonnenkalb et al. 2023, Lancet Microbe 4:e358 (BDQ+CFZ catalogue, 230 variants, structure)
- Snobre et al. 2024, AAC (PMC11619236; frameshifts, 24.8% susceptible, late stops/ARFs)
- Vargas et al. 2021, AAC (PMC8522733; mmpR+mmpS5/L5 double LoF, marker invalidity)
- Radhakrishnan et al. 2014, JBC (PMC4047419; MmpR5 structure 4NB5, DBD 34–99)
- Derendinger et al. 2023, Lancet Microbe (PMC10842724; 48% on-treatment acquisition)
- Ardizzoni et al. 2025, AAC (PMC12057352; acquired Rv0678 → outcome OR≈132; transient minorities)
- Brankin & Fowler 2023/2025 (PMC10072237, PMC12509863; minor alleles in catalogues, FRS 0.90→0.05)
- Hu et al. 2025, BMC Infect Dis (PMC12067902; pooled BDQ-R 5.7% among DR-TB)
- Sakhi et al. 2026 (PMC13325593; India 7.14% — single-source, label it)
- Kaniga et al. 2020 (PMC7098739; BDQ MGIT CC 1.0 µg/mL interim)
- EID 2021, 27(3) (DST capacity: 9/27 supranational labs do BDQ DST)
- PIB 9 Aug 2024 (BPaLM rollout under NTEP); NTEP BPaLM overview (ntep.in/node/6930)
- WHO catalogue v3, 10 Nov 2026 (who.int news; 9789240124653; expanded: pretomanid/cycloserine)
- Bates et al. 2021, RCPS (arXiv:2101.02703); Angelopoulos et al., LTT (arXiv:2110.01052),
  Conformal Risk Control (arXiv:2208.02814)
- Tibshirani et al. 2019, weighted conformal under covariate shift (arXiv:1904.06019)
- Barber et al. 2022, conformal beyond exchangeability (arXiv:2202.13415)
- Mozannar & Sontag 2020, learning-to-defer (arXiv:2006.01862)
- Smith et al. 2023, EPIG (arXiv:2304.08151)
- mBio 2025, confidence-based patient-level AMR conformal (10.1128/mbio.03431-25)
- Meier et al. 2021, NeurIPS (ESM zero-shot LLR); Brandes et al. 2023 (PMC10484790; ESM1b lookup)
- Vambe et al. 2026, Nat Comms (s41467-026-73551-w; Eswatini M146T founder outbreak)
- CRyPTIC compendium 2022, PLOS Biol (pbio.3001721; 12,289 isolates × 13 drugs incl. CFZ)