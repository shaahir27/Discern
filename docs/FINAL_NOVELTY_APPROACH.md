# BDQ-Triage — FINAL Novelty Approach & Metric-Beating Procedure

Research + verification date: 2026-11-18. Supersedes the novelty framing in `NOVELTY_UPGRADES.md`
(where they conflict, THIS file wins). Every number below is tagged:
**[verified]** = recomputed by us from downloaded public files this session;
**[reported]** = from a cited paper; **[unverified]** = must check before quoting.

Hard rule inherited from AGENTS.md: **no tuning to hit numbers.** Every experiment here is
pre-registered (write it into `docs/DECISIONS.md` BEFORE running), reports all rungs including
failures, and claims "beat" ONLY on identical or reconstructable data.

---

## 1. The one idea

> **We don't count mutations. We read HOW each one breaks the protein, check WHAT ELSE is
> broken in the same germ, and then prove the answer on germs the tool has never seen —
> scored on the metric that actually decides treatment: the false-alarm rate at the lab's
> own resistance prevalence.**

Three pillars, one engine:

1. **Mechanism-aware interpretation (biology).** Published work examines why some
   truncations have differing laboratory results and how gene interactions affect
   interpretation. Discern reports its archived evaluation of a downgrade to laboratory
   referral. [Snobre et al.](https://doi.org/10.1128/aac.00854-24) · [Schön et al.](https://doi.org/10.1128/aac.00266-25)

2. **External proof (evaluation).** No bedaquiline genotypic tool has been validated on a
   genuinely external, non-overlapping multi-cohort set. The "external" benchmark in the
   Fowler lab's own 2025 preprint is a same-consortium holdout; the one variant-level model
   that tried an external test saw specificity **collapse from 88% to 47.5%** (Rojas-Fuertes
   2026). We found and **verified** ready external cohorts (§4). This is the credibility
   unlock the whole field lacks.

3. **Certified, prevalence-aware triage (ML).** Beta-binomial shrinkage over mechanism
   classes stabilizes rare-variant evidence; conformal risk control turns the hand-picked
   "send to lab" thresholds into mathematical certificates, re-weighted to the user's
   prevalence; an information-gain queue tells the lab which experiment shrinks the unknown
   zone fastest. (Carried from NOVELTY_UPGRADES.md; still verified unpublished for TB.)

**Why it matters (India, now):** BPaLM rolled out under NTEP Sept 2024; bedaquiline DST
exists at ~2 national reference labs and takes up to 8 weeks (MSF Mumbai pilot); ICMR-NIRT's
Indian Mutation Catalogue v2 covers 14 drugs but **not bedaquiline** [reported, Claude report
§2, primary-source check pending]. At realistic prevalence (~2–7%), half of naive
"resistant" calls are false alarms — and a false R call can wrongly exclude a patient from
the best regimen. Mechanism-aware grading + external proof + prevalence-certified abstention
attacks exactly that failure mode.

---

## 2. What is NOT novel (retractions — do not claim these)

Fresh hostile prior-art audit killed several claims, including parts of my own earlier file:

| Claim | Status | Killer prior art |
|---|---|---|
| "First lineage/phylogeny-deconfounded grading" | **DEAD** | Kulkarni 2025 Nat Comms: kinship-PC-adjusted multivariable regression across 15 drugs incl. BDQ. **catomatic itself ships lineage/lab fixed effects + cluster random effects [verified from README].** Kulkarni 2026 CNN conditions BDQ MIC on a 62-SNP lineage barcode. |
| "Catalogues ignore biological context" | **DEAD (strawman)** | Miotto 2017 homoplasy criterion; Adlard 2025 manual phylogenetic-marker filter; Sonnenkalb 2023 reclassified 7 lineage-marker variants. |
| "Expand to new BDQ genes (mtrAB/lpqB/hadC)" | **DEAD + factually wrong** | WHO v2 BDQ graded list = **89 variants in Rv0678/atpE/pepQ ONLY** [verified GARC file]. Our `MUTATIONS_BDQ.csv.gz` has **only 5 genes** (mmpL5 84,815 rows; pepQ 5,235; Rv0678 3,080; atpE 1,059; mmpS5 887) [verified]. mtrAB/lpqB/hadC are not WHO v2 BDQ genes. |
| "First ML / minor-allele features / probabilities for BDQ" | **DEAD** | Ghosh 2026 BDQR-MTB (MLP/RF/SVM, AO/DP minor-allele features, probabilities, SHAP); Adlard 2025 + Brankin 2023 (minor alleles in catalogues). |
| "First prevalence-aware BDQ reasoning" | **DEAD** | Nimmo 2024 simulated PPV/NPV across prevalence; Saerens 2002 correction math. |
| "First Bayesian AMR grading" | **DEAD as invention** | ClinGen Bayesian framework; catomatic mixed-effects; Grazian 2023 Bayesian MIC mixtures. Beta-binomial mechanism-class shrinkage for TB is unpublished but is a *standard estimator choice*, not an invention — ablate it, don't headline it. |
| Per-lineage variant grades ("X=R in L4, U in L2") | **RACING WHO** | Kulkarni 2025 states lineage-distribution grading is under investigation for WHO v3. Sparse cells (most variant×lineage resistant counts 0–3) make it unidentifiable. Exploratory only. |

**Correction to NOVELTY_UPGRADES.md:** Vargas 2021 is **AAC** (PMC8522733), not Nat Comms.
Rv0678 coordinates: **NC_000962.3:778,990–779,487, + strand, 498 nt, GTG start, 165 aa**
[verified vs UniProt I6Y8F7] — not 7697–8194.

---

## 3. What IS novel (verified white space, defensible wording)

1. **Position/reading-frame-aware truncation features in a predictive model + decision layer.**
   Snobre 2024 is descriptive/structural; Sonnenkalb graded truncations by MIC not position;
   Rojas-Fuertes 2026 is missense-only (explicitly cannot handle LoF); Kulkarni pools LoF.
   Discern evaluates truncation position and alternative reading frame annotations in
   its archived research.

2. **Archived evaluation of an epistatic downgrade to laboratory referral.**
   WHO's third mutation catalogue and TB-Profiler are relevant prior art. Discern's E16 report
   documents the outcomes of its archived rule. Similar read fractions alone do not establish
   that variants are on the same genome. [WHO catalogue](https://www.who.int/publications/i/item/9789240124653) · [TB-Profiler paper](https://doi.org/10.1186/s13073-026-01767-y)
   Require HIGH-confidence LoF (frameshift/stop, excluding homopolymers); audit every flip;
   NEVER veto atpE-driven R (target-site resistance is pump-independent).

3. **True multi-cohort external validation of a BDQ triage tool** (§4). *Wording:* "validated
   on N non-overlapping external isolates from [cohorts], the first such evaluation for a
   bedaquiline genotypic triage tool (to our knowledge)." NOT "first externally validated"
   unqualified — KM assay is undocumented and the SA subset shares the training country.

4. **The benchmark nobody occupies:** calibration (ECE/Brier), PPV/NPV at 1–7% prevalence,
   cross-lineage sens/spec, patient-grouped CV as primary, abstention/coverage curves,
   external sens/spec. No competitor reports these (§5).

5. **Archived analyses of abstention and queue evidence.** E14 did not certify susceptible
   calls at its target, and bounds failed across sources. E15 compared the same number of
   variant sets, which exposed different numbers of sample labels across strategies.
   [E14](../results/approaches/conformal/e14_certificate_audit/REPORT.md) · [E15](../results/approaches/voi_queue/e15_label_reveal_sim/REPORT.md)

---

## 4. External validation data — VERIFIED this session

Recomputed from CRyPTIC v3.4.0 parquets (Zenodo 15680920) saved in `cryptic_v34/`:

| Cohort | n (R) | Genome | Lineage | Overlap w/ training | Assay / CC | Role |
|---|---|---|---|---|---|---|
| **KMVALIDATION2025** | **792 (35 R)** [verified] | all 792 | L2 559 / L4 198 / mixed 20 / L1 9 / L3 6 [verified] | **0** [verified] | **UNDOCUMENTED** (QUALITY=MEDIUM, METHOD_3 empty, CC all NaN) [verified] | Primary non-overlapping external set. Stratify SA vs China vs Italy; flag assay risk; sensitivity analysis only, not headline proof. Provenance: PRJEB58017 (SA MDR-TB), PRJNA436454 (China national survey), PRJNA903660 (Italy). |
| **Sobol 2026, KZN ZAF** (PMC12637405) | 632 sequenced; **423 with BDQ pDST**; 89 with resistance *mutations* [reported] | PRJNA1257273 public [verified accession] | available via calls | external | **MGIT 1.0 µg/mL** [reported] — **assay-matched to NICD training** | Strongest external anchor (same assay as 90% of our resistant training). **pDST-resistant count ≠ 89** (89 is genotypic); extract exact pDST-R from supplement [unverified]. 67% of resistant in transmission clusters. |
| **Choi 2026, Korea** (PMC13028822) | 49; 12 "elevated MIC" [reported] | PRJNA1301050 public [verified 49 runs] | — | external | 7H9 broth microdilution; elevated = ≥0.5 mg/L (**≠ WHO CC**) | Reclassify at 1.0 (MGIT-equiv) or 0.25 before comparing; small, supplementary tables. |
| New CRyPTIC-source isolates | 680 raw / 10 R [verified]; ~183 HIGH-quality / 4 R [reported] | partial | — | absent from training freeze | plate ECOFF 0.25 | Weak independence (same consortium, temporal holdout). Includes 41 Indian isolates, **all susceptible**. |

**Verified negatives (do not overclaim):**
- MUMBAI2018 / WHO2019 / NEJM2018 / LID2015 / PHE2019 have **ZERO bedaquiline rows** in
  v3.4.0 [verified]. **No independent Indian BDQ-resistant isolates exist in public data.**
  All 13 Mumbai BDQ-R are already in training. Indian claim = false-alarm audit only
  (site 5, ~5 R) + future Singla/NTEP collaboration. Singla 2025 (117 patients, 42 R) has
  **no WGS** [verified from PDF] — unusable without author collaboration.
- Training files contain **22,347 unique R/S IDs** (21,178 CRyPTIC + 1,169 NICD) [verified] —
  broader than the spec's **11,868 / 866** (which is QUALITY=HIGH + genome-join filtered).
  **Reconcile and stick to 11,868 for reproducibility**; document the filter.
- **46 isolates flip R↔S** between the training freeze and v3.4.0 (TMAS re-reading);
  54 cross-source R/S conflicts [reported by data agent] → label-noise exhibit (E17).
- v3.4.0 binary calls are **per-source at each source's own CC, NOT harmonized**
  (CRyPTIC plate 0.25 vs NICD MGIT 1.0 = 4-fold gap). Never pool external cohorts blindly;
  stratify by assay/CC. DATA_SCHEMA.pdf has no harmonization guidance [verified].

---

## 5. Competitor leaderboard & legitimate "beat" claims

| Competitor | Data n(R) | Protocol | BDQ metrics [reported] | Head-to-head? |
|---|---|---|---|---|
| WHO v2 / TB-Profiler | our 11,868 | as-is rules | binary sens **68.0** / spec 98.2 (Adlard anchor); U-as-S 49.7/67.9 | **YES, identical data** |
| Adlard 2025 catomatic-2b | 11,867 (866) | **shuffled 5-fold** | firm 94.2±0.4 / sens **79.4±1.8** / spec **98.5±0.3** | **YES, identical data** |
| Ghosh 2026 BDQR-MTB | 632 (**294 R = 46.5%**) | **random 80:20**, threshold tuning (not calibration), no CIs | MLP acc 83.6 / sens 86.3 / spec 81.3 / AUC 0.87; blind acc 84 / AUC 87.4 | **CONDITIONAL**: accessions in Tables S9/S10; code `PulmonomicsLab/bdqr-mtb-standalone` (FASTQ/VCF→pred+SHAP) [verified exists]. Re-call variants, reuse split, check training overlap first. |
| Kulkarni 2025 re-grading | 41,141 WHO subset | 50 kinship PCs, 5-fold | BDQ sens **58.9** (55.9–61.9) after allele-freq lowering | YES (farhat-lab/who-analysis public) but different cohort |
| Kulkarni 2026 CNN | 9,127 BDQ MIC | lineage-stratified 80/10/10 | 89% MICs within 1 dilution; **BDQ binary only TIES catalogue (p>0.05)** | partial (MtbQuantCNN public) |
| Nimmo 2024 Firth LR | 975 (21.9% R) | random 70/30 | AUROC **0.73** (0.70–0.76) | reconstructable, different cohort |
| MTB++ / AURA | CRyPTIC 6,224 / 12,185 | random splits | BDQ **F1 8.51 / 23.08** | YES (public code+data); low-value beat |
| Rojas-Fuertes 2026 | variant-level 62 missense | stratified 5-fold + WHOv2 external | AUC 0.826; **external spec 47.5%** | NO (variant-level ≠ isolate-level) |
| GenTB / Mykrobe | — | — | **no BDQ prediction** (GenTB explicitly excludes BDQ) | n/a |

**LEGITIMATE beats:**
- **vs WHO v2 on identical 11,868:** our binary sens 83.5 > 68.0. VALID. *Caution:* our
  ternary sens 73.5 does **not** beat WHO v2 at FRS 0.25 (76.9) — match FRS/abstention
  conventions before any ternary claim; mechanism features (M1) target closing this gap.
- **vs Adlard:** same data. Report BOTH their shuffled protocol and our harder patient-grouped
  one. A win = mechanism model beats their catalogue under both.
- **vs Ghosh:** only after reconstructing their 632 and replicating their split; expose
  random-split clonal leakage + 46.5% R enrichment + threshold-tuning-as-calibration.
- **vs MTB++/AURA:** trivial F1 beat, low value; optional.

**INVALID (never claim):** our AUC 0.909 (7.3% R) vs Ghosh 0.87 (46.5% R); our sens vs
Kulkarni 58.9 or Nimmo 0.73 (different cohorts); anything vs Rojas-Fuertes (different unit).

**Unreported axes = our leaderboard:** ECE/Brier, PPV/NPV @1–7%, cross-lineage sens/spec,
patient-grouped primary, coverage/abstention curves, external sens/spec. Nobody occupies it.

---

## 6. The metric-beating procedure (pre-registered ladder)

**Baselines (anchors, identical data):** WHO v2 as-is; Adlard catomatic-2b (reproduced);
current LR F2 (AUC 0.909, ternary firm 94.7 / sens 73.5 / spec 98.9).

**Method ladder — each rung ablated, NONE selected post-hoc:**
- **M0** baseline LR F2 (reproduction gate).
- **M1** + mechanism-aware truncation features. Build-time `models/lof_effects.csv` from
  Rv0678 CDS (778,990–779,487, GTG, 165 aa; fetch ~200 nt 3′ flank for ARF continuation).
  Classes: DBD-stop (≤99) / mid-stop (100–160) / late-stop (>160) / ARF-conserved (≥90%
  identity over dimerization region) / frameshift-position. Replace the single `rv_indel`/
  `rv_stop` max-FRS features with class-resolved features.
- **M2** + epistatic decision-layer veto (HIGH-confidence mmpL5/S5 LoF only; Rv0678-only R
  calls → U + `PUMP_DAMAGE_MODIFIER`; atpE-R never vetoed). Enumerate every flip with label.
- **M3** + beta-binomial shrinkage of per-variant evidence over mechanism classes (lineage as
  **covariate**, NOT stratifier). Targets the LOW_EVIDENCE rare-variant noise and the
  cross-source sensitivity drop.
- **M4** + conformal risk-controlled thresholds (LTT/Bentkus; prevalence-weighted
  certificates). Changes abstention quality, not AUC.
- **M5 (exploratory, pre-registered as such)** + lineage fixed effect; per-lineage grades
  only as a sparsity audit (expect unidentifiable; report as negative).

**Primary metrics (pre-register targets BEFORE running):**
1. Sensitivity at **fixed specificity 98.4%** (catalogue anchor) — beat 79.2%? Ceiling ~84%
   (731/866 carry a known-gene variant [verified run]; ~16% invisible caps it).
2. Specificity at **fixed sensitivity 80%** — beat 98.4% → drives PPV@2% from ~50% toward 70%+.
3. **PPV at 2% prevalence** (Bayes from CV sens/spec) — the clinical money metric.
4. **ECE / Brier** (we have cross-source ECE 0.288→0.013 [run]; nobody else reports).
5. **Cross-source NICD→CRyPTIC sensitivity** — beat 67.5% via M3 shrinkage.
6. **External sens/spec** per cohort: KM 792/35 (stratified, assay-flagged), Sobol 423
   (assay-matched, exact pDST-R from supplement), Choi 49 (reclassified). Headline first.
7. **Coverage/abstention curves** + novel-variant abstention (keep E08 AUC ~0.57 negative result).

**Decision rule:** ship the highest rung that improves the pre-registered primary metrics
under patient-grouped CV AND does not degrade external validation. If no rung beats baseline,
**report that honestly** and ship the baseline + external validation + certificates (still a
first). Never pick the best-looking number after the fact.

**Honesty guardrails:**
- Snobre's 24.8% susceptible-frameshift fraction may not replicate in our lineage mix —
  report the tension (E16 descriptive table stratified by first-stop region).
- KM assay/CC undocumented → sensitivity analysis, not proof; Sobol (MGIT 1.0) is the
  stronger anchor. Get Sobol's exact pDST-R count before quoting.
- Epistatic veto on noisy mmpL5 LoF can veto TRUE resistance → high-confidence LoF only,
  audit flips, never atpE.
- No Indian external R → no Indian-sensitivity claim, ever.
- "Beat" only on same/reconstructable data; external = "first to validate", not "beat".

---

## 6b. Optional P2 module — tNGS panel-coverage auditor ("audit the test, not just the germ")

**Status: bounded, optional, drop-in.** Reuses the existing variant table; no new model, no
new cohort. If time runs out, drop it — the core (§§1–6) stands alone. It widens the problem
statement from "interpret one genome" to "make the whole diagnostic pipeline smarter," and
gives the jury one memorable number.

**The broken process.** India is scaling **targeted NGS (tNGS)** — cheap, fast sequencing of
a *few chosen genes* instead of the whole genome (ICMR-NIRT pilot 2025; national tNGS workshop
Mar 2026). Panels are designed by **static catalogue mapping** (Rosendal 2026, PMC12903939) or
**manual literature curation** (Qu 2024, PMC11978391); one patent exists (WO2026074580A1).
WHO's 2023 tNGS guidance found evidence **insufficient for BDQ, pretomanid, LZD, CFZ, DLM**
(Jung 2026 review). Nobody scores a panel against **observed phenotypic resistance** in real
patients, weighted by local lineage. That is the white space.

**What it computes (pure arithmetic over existing tables):**
- Input: a panel's covered genomic regions (public amplicon/gene coordinates) + our resistant
  germs (866 R; report the HIGH-quality subset too) + catalogue R-grade variant positions
  (`MUTATIONS_BDQ` POSITION / AMINO_ACID_NUMBER → NC_000962.3 coordinates).
- **Panel sensitivity** = fraction of resistant germs carrying ≥1 R-grade variant *inside*
  the panel's covered regions. **Missed %** = 1 − panel sensitivity.
- **Variant-level coverage** = fraction of distinct R-grade variants whose position is covered.
- **Local-weighted missed %** = reweight by Indian-site / lineage-3 variant frequencies
  (lineage labels available via v3.4.0 join, §4) → "misses N% of *Indian-lineage* BDQ resistance."
- **Marginal gene gain** = extra resistant germs covered by adding gene g (per added amplicon).
  Rank genes → "add these 1–2 genes to catch M more."

**Demo number plan (pre-register; report whatever comes out):** score Deeplex Myc-TB vs
AmPORE-TB vs the Rosendal 2026 panel on our 866 R; print missed % and the top genes to add.
**Honest possible outcome:** panels may already cover Rv0678 well and miss little — then the
finding is "coverage is adequate; the real gap is *interpretation* and novel variants," which
is itself useful and feeds the core triage story. **Do not promise a large miss %.**

**Detection ≠ coverage caveat (state it):** tNGS depth limits minor-allele detection, and
Rv0678 is indel/homopolymer-rich (Hall 2024, eLife 98300) — a covered position can still be
*missed* by the caller. Flag indel-rich regions; coverage is an upper bound on detection.

**Public panel sources:** Deeplex Myc-TB (Genoscreen) 18 gene regions / 24-plex — coordinates
in the user manual (Illumina/Stop TB info sheet); AmPORE-TB (Oxford Nanopore) gene list;
Rosendal 2026 (PMC12903939) 30 amplicons / 19 regions / 20 drugs incl. Rv0678 (supplement).

**Novelty wording (defensible):** "first phenotype-weighted, local-lineage coverage audit of
commercial tNGS panels against 11,868 real genomes (to our knowledge)." Cite the static-design
adjacency above. **FORBID:** "we invented tNGS"; "our panel is clinically validated"; "WHO
endorses our panel"; any rupee/cost-saving claim; "guaranteed detection" (coverage ≠ detection).

**Build:** ~0.5 day. `eval/panel_audit.py` (region-overlap join + marginal-gain ranking) →
`results/panel_audit.json` + a `/panels` report page. New test: a panel covering everything
returns 100% sensitivity; a panel covering nothing returns 0% (negative controls).

---

## 7. Build plan (3-day on-site + pre-work)

**Pre-work (now, before Dec 18):** download KM parquets (done, in `cryptic_v34/`); fetch
Sobol supplement + PRJNA1257273 and Choi Tables 1–4; build `lof_effects.csv`; reconcile
22,347 vs 11,868; pre-register ladder + metrics in DECISIONS.md; (optional) pull Ghosh S9/S10
accessions and test variant re-calling on demo SRR28888046.

| Day | Work | Gate |
|---|---|---|
| 1 | M1 mechanism features + M2 veto + E16 ablation (descriptive first); external-validation harness on KM + Sobol | golden reproduction byte-identical with enhancements OFF; flips enumerated |
| 2 | M3 shrinkage + M4 conformal; same-data head-to-heads (WHO v2, Adlard shuffled+grouped); new axes (ECE, PPV@2%, coverage) | pre-registered metrics written; all rungs reported |
| 3 | Ghosh head-to-head IF data ready; UI (mechanism evidence cards, certificate card, external-validation page, info-gain queue); docs/demo; `make verify` | demo runs offline; every number traces to results/ |
| P2 | **tNGS panel-coverage auditor (§6b)** — only if Days 1–3 green; `eval/panel_audit.py` + `/panels` page | negative controls pass; missed-% reported as-is |

New tests: LOF-table determinism; enhanced-vs-golden isolation; veto never fires on atpE-R;
external cohorts have zero training-ID overlap (assert); conformal math unit tests;
pre-registration file timestamp precedes results.

---

## 8. Pitch (5-minute spine)

1. **Problem (30 s):** BPaLM is India's new DR-TB backbone; BDQ DST = 2 labs, 8 weeks;
   ICMR's own catalogue omits BDQ. At 2% prevalence half of naive R calls are false alarms —
   and a false R can wrongly deny a patient the best regimen.
2. **Biology engine:** a late-stop Rv0678 truncation → "protein likely functional (Snobre
   2024)", U not blind R; a pump-broken double mutant → epistatic veto card. Side-by-side
   with the pool-all-LoF behavior of current catalogues.
3. **External proof (the headline):** "Every BDQ tool until now validated on its own
   training consortium. We tested on 792 isolates from South Africa/China/Italy it never
   saw, plus 423 MGIT-tested KZN isolates — here are the numbers, including where it fails."
   Show the Rojas-Fuertes 88%→47.5% external collapse as why this matters.
4. **Certainty + prevalence:** unresolved c-11a → U with conformal certificate; slide
   prevalence 7%→2%, watch PPV fall and the certificate re-weight.
5. **Priority engine:** info-gain queue + simulation curve ("phenotype these 10 sets →
   X% of the queue resolves").
6. **Honesty moat:** performance page — patient-grouped CV, cross-source swing, E08
   novel-variant NEGATIVE result, external cohorts, every number regenerated by `make verify`.
7. **(Optional P2) Audit the test itself:** "Commercial tNGS panel X misses N% of Indian-lineage
   BDQ resistance; add gene Y to catch M more — arithmetic from 11,868 genomes." Shows you
   understand the whole pipeline, not just one step (§6b).

**Positioning sentence:** "The first mechanism-aware, externally-validated bedaquiline
triage engine — it reads how a mutation breaks the protein, proves itself on strangers, and
certifies when to say 'I don't know'."

**Claims policy additions:** ALLOW "first predictive integration of truncation-position/ARF
reasoning for BDQ (to our knowledge)"; "audited epistatic downgrade rule"; "validated on
N non-overlapping external isolates [cohorts, assay caveats]"; "beat WHO v2 binary sensitivity
on identical data". FORBID "first lineage-deconfounded" (Kulkarni/catomatic); "first ML/
minor-allele/probability model" (Ghosh/Adlard); "first deployed epistatic veto" (TB-Profiler v6.6.5/WHO v3); "catalogues ignore context"; "works on Indian
patients"; "predicts new variants"; any cross-dataset "beat"; unqualified "first externally
validated"; guaranteed NPV/PPV without bound+δ.
- Panel auditor (§6b): ALLOW "first phenotype-weighted, local-lineage coverage audit of
  commercial tNGS panels (to our knowledge)". FORBID "we invented tNGS"; "our panel is
  clinically validated"; "WHO endorses our panel"; rupee/cost-saving claims; "guaranteed
  detection" (coverage ≠ detection).

---

## 9. Verified references (this session)

- CRyPTIC v3.4.0: Zenodo 15680920 (DST_MEASUREMENTS/GENOMES parquets downloaded + recomputed).
- KMVALIDATION2025: 792/35, zero overlap, MEDIUM/undocumented CC [verified]; PRJEB58017,
  PRJNA436454, PRJNA903660.
- Sobol 2026 AJRCCM PMC12637405: 632 sequenced / 423 pDST / 89 genotypic-R; MGIT 1.0;
  PRJNA1257273; 67% transmission clusters [reported+accession verified].
- Choi 2026 Pathogens PMC13028822: 49 isolates, 12 elevated-MIC (≥0.5, ≠WHO CC);
  PRJNA1301050 [verified].
- Ghosh 2026 LSA 9(7):e202503539, PMC13135273; code PulmonomicsLab/bdqr-mtb-standalone
  [verified]; 632 (294 R), random 80:20, AUC 0.87.
- Adlard 2025 Microb Genom 10.1099/mgen.0.001429, PMC12177157; shuffled CV 94.2/79.4/98.5;
  follow-up bioRxiv 2025.10.02.679941, repo cryptic-catalogues-2025.
- Kulkarni 2025 Nat Comms PMC11876447 (50 kinship PCs, BDQ sens 58.9); Kulkarni 2026
  PMC13315309 (lineage CNN, BDQ binary ties catalogue); farhat-lab/who-analysis; MtbQuantCNN.
- catomatic README: lineage/lab fixed effects + agglomerative-cluster random effects [verified].
- Snobre 2024 AAC PMC11619236 (24.8% frameshift susceptible; late stops >160; ARFs; DBD 34–99).
- Vargas 2021 **AAC** PMC8522733 (82 double-LoF; mmpR markers invalid when pump broken).
- Sonnenkalb 2023 Lancet Microbe PMC10156607 (265 variants; 7 lineage-marker reclassifications).
- Rojas-Fuertes 2026 bioRxiv 10.64898/2026.02.12.705543 (missense-only RF, AUC 0.826,
  external spec 47.5%).
- Nimmo 2024 PMC11072239 (Firth LR AUROC 0.73); MTB++ Bioinformatics 40:i39 (BDQ F1 8.51);
  AURA bioRxiv 2026.04.25.720842 (BDQ F1 23.08).
- WHO v2 GARC: BDQ = 89 variants, Rv0678/atpE/pepQ only [verified]; WHO v3 10 Nov 2026
  (9789240124653; adds pretomanid/cycloserine, NOT BDQ headline).
- Rv0678 = NC_000962.3:778,990–779,487, +, 498 nt, GTG, 165 aa; UniProt I6Y8F7 [verified].
- Miotto 2017 ERJ PMC5898944 (homoplasy criterion); Brankin 2023 PMC10072237 (minor alleles).
- Conformal: Bates 2021 arXiv:2101.02703; Angelopoulos LTT arXiv:2110.01052; Tibshirani 2019
  arXiv:1904.06019; mBio 2025 10.1128/mbio.03431-25 (E. coli prior art). EPIG arXiv:2304.08151.
- India: BPaLM NTEP Sept 2024 (PIB 9 Aug 2024); ICMR-NIRT catalogue v2 omits BDQ; MSF Mumbai
  pilot ($230 seq vs $119 DST, 8 weeks); Singla 2025 CID (117/42 R, no WGS) [verified PDF].
- tNGS panel auditor (§6b): WHO tNGS guidance 2023 (iris.who.int 18b79164-04da-4a25-90ec-ccf944c57164);
  Rosendal 2026 PMC12903939 (static catalogue-mapped panel, 30 amplicons/19 regions/20 drugs);
  Qu 2024 PMC11978391 (manual curation); Jung 2026 review (MDPI 2227-9059/14/1/93; BDQ evidence
  insufficient); Deeplex Myc-TB info sheet (stoptb.org); patent WO2026074580A1; Hall 2024
  eLife 98300 (Rv0678 homopolymer/indel calling errors); ICMR-NIRT tNGS pilot 2025
  (Microbiol Spectrum 10.1128/spectrum.02917-25).
