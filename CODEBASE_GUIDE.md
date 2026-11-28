# Discern — Complete Codebase Guide & Zero-to-Expert Reference

> **Research use only.** Discern is an offline research and clinical-triage workbench for studying *Mycobacterium tuberculosis* resistance to **bedaquiline (BDQ)** from DNA mutations. It is **not** a certified diagnostic device and does not prescribe treatment. Every call must be interpreted alongside culture drug-susceptibility testing (DST) and national TB guidelines.

---

## Quick Navigation: How to Read This Guide

- **New to TB genetics or machine learning?** Read **Sections 1, 2, and 6** first — they use plain-English analogies and step-by-step arithmetic to explain every biological and mathematical idea from scratch.
- **Want to trace what happens when you upload a file or call `/predict`?** Read **Sections 3, 5, and 8** (end-to-end architecture + two complete worked walkthroughs).
- **Looking for a specific Python module, HTML template, evaluation script, test, or result file?** Jump straight to **Section 10**, which documents **every single file in the entire repository** one by one.
- **Ready to run the server, CLI, Docker container, or test suite?** See **Section 11**.

---

## 1. The Big Picture: What Problem Does Discern Solve?

### 1.1 Why Tuberculosis Labs Need Speed *and* Humility

Bedaquiline (BDQ) is the backbone antibiotic of modern all-oral regimens for multidrug-resistant tuberculosis (MDR-TB). To check whether a patient's TB strain resists bedaquiline, laboratories traditionally grow the live bacteria in liquid culture (such as **MGIT960**) or microtiter plates (**UKMYC**) and expose them to the drug. Culture DST is trustworthy, **but *M. tuberculosis* divides roughly once every 18–24 hours**, so culture takes **6 to 8 weeks**. During those weeks, a patient can stay on a failing drug regimen or receive toxic backup drugs unnecessarily.

DNA sequencing (reading the bacteria's genome or targeted resistance genes directly) takes **hours to days instead of months**. However, turning raw DNA mutations into a clinical resistance call has a dangerous catch:

- **Lopsided cost of error**:
  - Falsely calling a harmless mutation **Resistant (`R`)** pulls the patient off bedaquiline onto longer, harsher, less effective regimens.
  - Falsely calling a real resistance mutation **Susceptible (`S`)** leaves the patient on a failing drug and risks further resistance and transmission.
- **Many mutations are rare or never seen before**: In a single TB regulatory gene (`Rv0678`), hundreds of distinct mutations exist, and new ones keep appearing across clinics. An ordinary machine-learning classifier that forces every DNA sample into a binary `R` or `S` box will blindly guess on mutations it has never seen.

Discern solves this with a strict **two-tier triage architecture** and three possible outputs:

| Call | Label Shown in UI | Meaning & Action |
|---|---|---|
| **`R`** | **Resistant** | A statistically supported catalogue resistance variant is present (and not vetoed by an inactivated efflux pump). |
| **`S`** | **Susceptible** | Only wild-type sequence or benign variants are observed (still includes a caveat when no core variants are present). |
| **`U`** | **Send for lab test** *(Referral)* | Variants are novel, rare, contradictory, or low-confidence. **Discern refuses to guess.** It sends the sample to the wet lab and uses calibrated machine learning to rank which `U` samples deserve culture testing first. |

### 1.2 Two Non-Negotiable Architecture Rules

1. **Tier A (the curated rule catalogue) always makes the final call (`R` / `S` / `U`).** The machine-learning model (Tier B) **never overrides** the catalogue call. It only provides:
   - a calibrated probability of resistance adjusted to the local clinic's resistance prevalence,
   - a priority percentile ranking for unresolved (`U`) samples in `/queue`, and
   - a `MODEL_DISAGREES` warning flag when Tier A and Tier B strongly disagree.
2. **100% offline, deterministic, and non-logging at runtime.** There are zero external HTTP calls, telemetry pings, CDNs, external Google Fonts, or request-body access logs (`--no-access-log`). It is designed to run inside an air-gapped TB reference lab or district hospital server.

### 1.3 Live Runtime vs. Offline Evaluation vs. Exploratory R&D (Avoiding Confusion)

A common mistake when opening this repository is confusing modules that power the **live web app** with scripts that belong strictly to **offline scientific audits** or **unwired research prototypes**. Keep these three lanes separate:

| Lane | What Runs Here | Where It Lives |
|---|---|---|
| **1. Live Production Workbench** | Loads frozen `models/lr_model/*` and `models/lr_mechanism/lof_effects.csv`. Runs native lab-export import (`importers.py`), indel canonicalization (`indel_canon.py`), panel blind-spot audit (`panels.py`), Tier A catalogue + live high-confidence epistatic veto (`decision.py`), Tier B L2 logistic regression (`features.py`, `model.py`), Platt & Bayes prevalence math (`prevalence.py`, `label_shift.py`), evidence cards (`evidence.py`), domain diagrams (`domains.py`), quality checks (`quality.py`), batch surveillance (`surveillance.py`), and longitudinal timelines (`longitudinal.py`). | [src/discern/decision.py](src/discern/decision.py), [src/discern/api/app.py](src/discern/api/app.py), [models/lr_model/](models/lr_model/) |
| **2. Offline Evaluation & Verification Suite** | Rebuilds cohorts from raw CRyPTIC/NICD/KM tables, validates splits `S1–S4`, executes benchmarks `E01–E24`, writes immutable JSON/CSV/MD archives into `results/`, regenerates `docs/MODEL_CARD.md` and `docs/RESULTS.md`, builds production artifacts into `models/`, and verifies acceptance gates (`14/14` in `eval/verify.py`). | [eval/](eval/), [results/](results/), [eval/verify.py](eval/verify.py) |
| **3. Exploratory / Preserved Research Modules** | Scientifically evaluated and deliberately **not shipped** in live calling after negative or non-superior benchmark outcomes — or kept as standalone R&D prototypes. Includes conformal certificates ([src/discern/conformal.py](src/discern/conformal.py)), empirical-Bayes shrinkage ([src/discern/shrinkage.py](src/discern/shrinkage.py)), Beta entropy value-of-information ([src/discern/voi.py](src/discern/voi.py)), 3D PDB structural/biophysical exploration ([src/discern/structure.py](src/discern/structure.py), [eval/e22_biophysics_full_cohort.py](eval/e22_biophysics_full_cohort.py)), and standalone lineage/MIC prototyping ([eval/e23_novel_pipeline.py](eval/e23_novel_pipeline.py)). | [src/discern/conformal.py](src/discern/conformal.py), [src/discern/shrinkage.py](src/discern/shrinkage.py), [src/discern/voi.py](src/discern/voi.py), [src/discern/structure.py](src/discern/structure.py), [eval/e23_novel_pipeline.py](eval/e23_novel_pipeline.py) |

> **Important Configuration Nuance (`config.yaml` vs Live `Predictor`):** In [config.yaml](config.yaml), the offline benchmark reproduction section sets `enhanced: {lof_policy: off, mechanism_features: false, epistatic_veto: true, shrinkage: false, conformal: false}` so that `eval/verify.py` reproduces baseline acceptance numbers (`E01–E13`) without modification. Meanwhile, in the **live application**, `Predictor` in [src/discern/decision.py](src/discern/decision.py) loads `models/lr_mechanism/lof_effects.csv` and applies an **audited static high-confidence epistatic veto** (`veto_high_conf`) that flips `Rv0678`-only `R` calls to `S` when accompanied by a pre-audited truncating `mmpL5`/`mmpS5` mutation at comparable frequency (updating live catalogue sensitivity/specificity to `0.7918` / `0.9872` from `e16_lof_ablation`).

### 1.4 Beginner's Glossary

| Term | Plain-English Meaning |
|---|---|
| **Isolate** | One biological TB sample collected from a patient, cultured or sequenced. One patient (`SUBJID`) can contribute multiple isolates over time. |
| **GARC Variant Format** | Standardized mutation string `GENE@MUTATION`, e.g., `Rv0678@141_ins_c` or `atpE@A63P`. |
| **FRS (Fraction of Reads Supporting) / VAF** | Between `0.0` and `1.0`. Measures what fraction of sequencing reads carry the mutation. `FRS = 0.25` means 25% mutant and 75% wild-type bacteria in the same sputum specimen (**subclonal heteroresistance**). |
| **Wild-type (`WT`)** | The normal, unmutated reference TB strain (**H37Rv**, GenBank accession `NC_000962.3`). |
| **WHO v2 Catalogue / `piezo`** | The World Health Organization's standardized mutation lookup table (`piezo` is the Python parser library; `catomatic` builds catalogues from raw binomial counts). |
| **Firm Call vs. `U`** | `R` and `S` are **firm calls**. `U` ("Undetermined / Send for lab test") is an intentional abstention when evidence is missing, scarce, or ambiguous. |
| **Sensitivity / Specificity** | **Sensitivity**: fraction of truly resistant isolates caught as `R`. **Specificity**: fraction of truly susceptible isolates safely cleared as `S`. |
| **PPV / NPV** | **Positive Predictive Value**: if Discern says `R`, the probability the germ is actually resistant given your clinic's true prevalence. **Negative Predictive Value**: if Discern says `S`, the probability it is actually susceptible. |
| **ROC AUC vs. AUPRC** | **AUC** (`0.5` = coin flip, `1.0` = perfect) measures ranking across all thresholds. **AUPRC** (area under precision-recall curve) focuses specifically on the rare positive (`R`) class, making it much stricter on imbalanced clinical cohorts (`7.3%` resistant). |
| **Platt Scaling** | Fitting a 1-D logistic curve on top of raw log-odds so model output numbers behave as true calibrated probabilities. |
| **Epistasis** | When a mutation in Gene B suppresses or cancels the biological effect of a mutation in Gene A (here, a broken `mmpL5` pump neutralizes a broken `Rv0678` repressor). |

---

## 2. Biological Foundations (No Biology Degree Needed)

### 2.1 The Five Bedaquiline-Relevant Genes

Discern reads variants across five loci on the *M. tuberculosis* H37Rv (`NC_000962.3`) chromosome:

| Gene | Biological Role | How Mutations Affect Bedaquiline | Role in Discern |
|---|---|---|---|
| **`Rv0678`** (`mmpR5`) | Transcriptional **repressor** protein (`165` amino acids; DNA-binding domain spanning codons `34–99`, recognition helix α4 at `62–68`). Normally sits like a padlock on DNA to keep the `mmpS5-mmpL5` pump turned *off*. | Loss-of-function mutations (frameshifts, stops, DBD missense, promoter `c-11a`) break the padlock, turning the pump *on* and pumping bedaquiline out of the cell → **Resistance (`R`)**. | **Primary driver** of BDQ resistance in clinical surveillance (`~88.5%` of known-locus resistant strains carry `Rv0678` mutations). |
| **`atpE`** | Encodes the c-subunit of **ATP synthase** — the actual molecular motor bedaquiline binds and shuts down. | Target-site missense changes prevent bedaquiline binding → high-level **Resistance (`R`)**. | **Core feature gene**. Never vetoed by `mmpL5` (bypass resistance does not need the efflux pump). |
| **`pepQ`** | Cytoplasmic peptidase gene. | Specific mutations cause low-to-moderate cross-resistance alongside `Rv0678` pathways. | **Core feature gene**. Never vetoed by `mmpL5`. |
| **`mmpL5`** | Large transmembrane **efflux pump** channel that physically expels bedaquiline from the bacterium. | Frame-disrupting **loss-of-function (LoF)** indels/stops destroy the pump itself, **re-sensitizing** the germ even if `Rv0678` is broken! | **Excluded from Tier A positive rules**; used for Tier B `mmpL5_LoF` interaction features and the **Epistatic Veto**. |
| **`mmpS5`** | Accessory membrane fusion partner essential for the `mmpS5-mmpL5` pump complex. | Truncating/frameshift mutations likewise abolish efflux activity. | Co-monitored for `mmpL5/mmpS5` epistatic LoF and panel coverage auditing. |

### 2.2 Epistasis ("The Broken Pump Rule")

Imagine bedaquiline is water entering a basement, **`mmpL5-mmpS5` is the bilge pump**, and **`Rv0678` is the padlock keeping the pump's switch OFF**:

1. **Normal germ**: Padlock (`Rv0678`) is intact → pump stays off → bedaquiline builds up → bacterium dies (**Susceptible `S`**).
2. **Broken padlock (`Rv0678` LoF alone)**: Switch springs open → pump turns high → drug pumped out → bacterium survives (**Resistant `R`**).
3. **Broken padlock AND broken pump (`Rv0678` LoF + `mmpL5`/`mmpS5` LoF)**: Even with the switch wide open, a snapped pump pipe cannot pump water out → bedaquiline stays inside → bacterium dies (**Susceptible `S`**)!

**Real biological catch implemented in [src/discern/mechanism.py](src/discern/mechanism.py):**
- In mixed infections, one bacterial subclone might carry `Rv0678` at `90%` frequency while an unrelated minor subclone carries `mmpL5` LoF at `15%`. For epistasis to apply in the same cell line, the pump-breaking mutation must co-exist on a matching clonal fraction:
  `FRS(mmpL5 or mmpS5 LoF) >= FRS(Rv0678) - 0.20` (WHO 3rd Ed. 2026 / Snobre et al. 2024 linkage rule).
- In the live predictor ([src/discern/decision.py](src/discern/decision.py)), only pre-audited **high-confidence** truncating variants recorded in `models/lr_mechanism/lof_effects.csv` can trigger the live veto, and **only** when `Rv0678` is the sole resistance call (any `R` variant in `atpE` or `pepQ` is preserved).

### 2.3 Why Low Allele Frequency (`FRS >= 0.10`) Changes Everything

Routine bioinformatics pipelines historically filtered mutations unless at least `75%` of sequencing reads carried them (`FRS >= 0.75`). In tuberculosis lungs, patients often harbor heterogeneous subpopulations where a new `Rv0678` mutant is still expanding at `15%–50%` of reads.
- Raising the threshold to `0.75` blinds the WHO v2 catalogue (sensitivity among firm calls collapses from **`0.801`** at `FRS >= 0.1` down to **`0.561`** at `FRS >= 0.75`, and `U-as-S` sensitivity drops from `0.679` to `0.497`).
- Lowering the inclusion floor to **`FRS >= 0.10`** recovers subclonal resistance, while triggering a clinical **`HETERORESISTANCE_SUBCLONAL`** alert for `0.10 <= FRS < 0.75` (explaining why fast 14-day plate tests sometimes look susceptible before resistant colonies outgrow in liquid MGIT960 culture).

### 2.4 The ~15.8% "Invisible Resistance" Floor

In Discern's audited training cohort of `866` bedaquiline-resistant isolates (`results/approaches/panel_audit/e24_panel_coverage/results.json`), **`729`** carry at least one mutation across the five known loci, while **`137` isolates (`15.8%`) are phenotypically resistant despite having ZERO variants in `Rv0678`, `atpE`, `pepQ`, `mmpL5`, or `mmpS5`**. Because no five-gene assay can see mechanisms outside those five genes, calling `S` with zero variants always attaches `NO_VARIANT_INVISIBLE_RESISTANCE_POSSIBLE`.

### 2.5 Sequencing Assay & tNGS Panel Blind Spots ([src/discern/panels.py](src/discern/panels.py))

Not all hospital sequencers read all five genes. Targeted next-generation sequencing (tNGS) panels trade whole-genome breadth for sputum speed and cost:

| Panel Key | Assay / Platform | Genes Covered | VAF Floor | Known-Locus Ceiling | Overall BDQ-R Ceiling | Main Blind Spots |
|---|---|---|---|---|---|---|
| `wgs` | Whole Genome Sequencing | `Rv0678, atpE, mmpL5, mmpS5, pepQ` | `0.10` | `100.0%` | **`84.18%`** | None (beyond the 15.8% invisible loci) |
| `deeplex_xl` | GenoScreen Deeplex Myc-TB XL (beta) | All 5 loci | `0.03` | `100.0%` | **`84.18%`** | None among known 5 loci |
| `ampore_tb` | ONT AmPORE-TB / `wf-tb-amr` | `Rv0678, atpE` | `0.15` | `94.51%` | **`79.56%`** | Blind to `pepQ` (`5.49%` loss) and `mmpL5`/`mmpS5` (`EPISTATIC_VETO_BLIND`, homopolymer risk) |
| `deeplex_myc_tb` | GenoScreen Deeplex Myc-TB | `Rv0678` only | `0.03` | `88.48%` | **`74.48%`** | Blind to `atpE`, `pepQ`, `mmpL5`, `mmpS5` (`11.52%` known loss) |
| `icmr_nirt_minion` | ICMR-NIRT India MinION Pilot | `Rv0678` only | `0.15` | `88.48%` | **`74.48%`** | Blind to `atpE`, `pepQ`, `mmpL5`, `mmpS5` + Nanopore homopolymer indel caution |

Nanopore platforms (`ampore_tb`, `icmr_nirt_minion`) also trigger **`HOMOPOLYMER_ARTIFACT_RISK`** when encountering indels inside poly-G/poly-C stretches such as the well-known `Rv0678` `138_ins_g` / `141_ins_c` / `192_ins_g` / `193_del_g` homopolymer hotspots ([src/discern/mechanism.py](src/discern/mechanism.py), [src/discern/panels.py](src/discern/panels.py)).

---

## 3. End-to-End System Architecture

```text
====================================================================================================
PART A: OFFLINE TRAINING, EVALUATION & VERIFICATION PIPELINE  (eval/ -> data/, models/, results/)
====================================================================================================

  [data/raw/tb-bdq-cat-main/] (CRyPTIC + NICD + WHO v2 + wildcards + Rv0678_WT.pdb + MANIFEST.json)
               │
               ▼
       eval/e01_prep.py + src/discern/prep.py
               │   • Keep CRyPTIC QUALITY=="HIGH" + all NICD (MGIT960)
               │   • Deduplicate UNIQUEID (prefer PHENOTYPE=="R", non-null MIC), inner join GENOMES
               │   • Build 11,868 isolates (866 R, 11,611 SUBJID) + 25,554 mutation rows
               ▼
  [data/processed/samples.parquet & mutations.parquet]
               │
               ├──► splits.py (S1 GroupKFold by SUBJID, S2 LOSO by SITEID, S3 cross-source, S4 variant-set)
               │
               ├──► eval/e02..e24 & ext_km_validate.py ──► results/{models,approaches}/**/results.json
               │                                               └──► eval/verify.py (14/14 acceptance gates)
               │                                               └──► eval/generate_docs.py -> docs/{MODEL_CARD,RESULTS}.md
               │
               └──► eval/build_production.py (fit on ALL 11,868 samples; nested inner calibration)
                           │
                           ├──► models/lr_model/catalogue_bdq.csv + catalogue_version.json
                           ├──► models/lr_model/bdq_lr.joblib + metadata.json (120 features; lo/hi/t_bin)
                           ├──► models/lr_model/evidence_store.csv (698 variant summaries w/ CRyPTIC vs NICD)
                           └──► models/lr_mechanism/lof_effects.csv (audited LoF lookup for live veto)


====================================================================================================
PART B: OFFLINE RUNTIME WORKBENCH  (src/discern/api/app.py & cli.py -> port 8000 / batch CSV)
====================================================================================================

  [Lab Inputs: Standardized 4/6-col CSV · JSON · TB-Profiler JSON · Deeplex CSV · AmPORE-TB TSV · Replay ID]
               │
               ├──► Lane 1 (Public Replay Archive): results/oof_predictions.parquet & demo/held_out_samples.csv
               │    (Shows held-out CV score so demo samples never display in-sample training optimism!)
               │
               └──► Lane 2 (Live Predictor + Session-Scoped Upload Store — zero disk logging):
                           │
                    importers.py + indel_canon.py + grammar.py + quality.py
                    (Auto-detect vendor export, normalize HGVS -> GARC, canonicalize indels, score upload QC)
                           │
                    ┌──────┴──────────────────────────────────────────────────────┐
                    ▼                                                             ▼
         TIER A: CATALOGUE + EPISTATIC VETO                          TIER B: CALIBRATED L2 LOGISTIC MODEL
    • Match non-syn, non-null muts (FRS >= 0.1)                 • Build F2 (111 vocab + 9 aggregates = 120 cols)
    • piezo rules (R > U > S; mmpL5/S5 excluded)                • LogisticRegression(C=0.3) -> raw score
    • Epistatic Veto: if Rv0678-only R AND high-conf            • Nested Platt scaling (C=1e6) -> p_train (at 7.30%)
      mmpL5/S5 LoF in lof_effects.csv w/ FRS >= Rv-0.20         • Odds-ratio shift -> p_adj (at user or MLLS pi)
      -> flip call to S                                         • Percentiles + model/bin thresholds (t_bin, lo, hi)
                    └──────┬──────────────────────────────────────────────────────┘
                           ▼
             POST-PROCESSING & DECISION ENVELOPE (decision.py, evidence.py, domains.py, panels.py)
    • Call is ALWAYS Tier A's call (Model NEVER overrides R/S/U)
    • Bayes PPV/NPV at prevalence pi using audited CV catalogue rates (sens=0.7918, spec=0.9872)
    • Assay blind-spot warnings, clinical reflex cards, SVG domain tracks, batch recurrence & timeline
                           │
                           ▼
        Jinja2 HTML Templates (/result, /queue, /patient/{id}, /panels, /priority, /performance, /model-card)
```

---

## 4. Dataset, Splitting Hygiene & Variant Grammar

### 4.1 The Training Cohort (`E01`) & External Benchmarks

| Cohort | Samples | Resistant (`R`) | Preval. (`% R`) | Notes |
|---|---|---|---|---|
| **CRyPTIC** (`QUALITY == "HIGH"`) | `10,704` | `83` | `0.78%` | Global surveillance (UKMYC 14-day microtiter broth); multi-site including Mumbai (`Site 04`: `950`, `5 R`) and Peru (`Site 05`: `3,035`, `5 R`). |
| **NICD** (South Africa, `MGIT960`) | `1,164` | `783` | `67.27%` | Heavily resistance-enriched programmatic collection (`SITEID == 10`). |
| **Pooled Training Cohort (`samples.parquet`)** | **`11,868`** | **`866`** | **`7.30%`** | `11,611` unique patients (`SUBJID`). `25,554` total mutation rows (`1,115` minor-allele rows, median `FRS ~ 0.38`; `11` null `Z` rows across `4` samples). |
| **External KM (`KMVALIDATION2025`)** | `792` | `35` | `4.42%` | Strictly held-out CRyPTIC v3.4.0 genotypes, never used during training or threshold calibration. |

**Why the huge difference between CRyPTIC (`0.78%` R) and NICD (`67.3%` R) matters:**
Because `783` of `866` resistant isolates (`90.4%`) come from NICD's enriched archive, an unadjusted classifier trained on the pooled cohort assumes roughly `1 in 14` TB patients has BDQ resistance (`pi_tr = 7.30%`). In a typical surveillance clinic at `1%–5%` prevalence, raw probabilities over-predict risk unless re-calibrated via **odds-ratio prevalence shift** (Section 6). Similarly, every variant evidence card in `evidence_store.csv` reports CRyPTIC and NICD counts separately so a user is not fooled by high resistance percentages driven by NICD enrichment.

### 4.2 Four Anti-Leakage Cross-Validation Splits ([src/discern/splits.py](src/discern/splits.py))

One patient can have 2–5 serial isolates sampled months apart with identical resistance mutations. Splitting randomly by isolate would put Tuesday's sputum isolate in Training and Friday's identical isolate in Test. Discern defines four strict splits:

| Split ID | Definition | What Stress Test It Measures |
|---|---|---|
| **`S1`** | `GroupKFold(n_splits=5)` grouped by `SUBJID` | **Main benchmark**: same patient never appears in both train and test. (`g7_adlard` proved naive shuffled 5-fold inflated sensitivity from `0.7924` to `0.7955`). |
| **`S2`** | Leave-One-Site-Out (`SITEID`, `23` distinct sites; `11` main sites in `e07_loso`) | **Geographic generalization**: hold out an entire hospital/country. (`NICD` shares `SITEID == 10` with CRyPTIC South Africa.) |
| **`S3`** | Cross-Source (`CRyPTIC -> NICD` and `NICD -> CRyPTIC`) | **Severe shift test**: train strictly on global UKMYC and evaluate on enriched SA MGIT960, and vice versa. |
| **`S4`** | Variant-Set Grouping (`"\|".join(sorted(MUT))`) | **Novel-variant test**: test samples share zero exact mutation signatures with training samples (`eval/e08_novel_variant.py`). |

> **Golden Rule of Leakage Prevention:** Inside every training fold of `S1`, `S2`, and `S3`, the variant vocabulary (`vocab_from`), the logistic regression fit, the inner calibration logistic model (Platt scaling), and the probability cutoffs (`t_bin`, `lo`, `hi`) are fit **100% inside the training indices** (`eval/verify.py` and `tests/test_archive_checks.py` verify ID disjointness).

### 4.3 GARC DNA Variant Grammar & Indel Canonicalization ([src/discern/grammar.py](src/discern/grammar.py), [src/discern/indel_canon.py](src/discern/indel_canon.py))

Every variant enters the core engine in **GARC** format (`GENE@MUTATION`), parsed into exact boolean flags by [`grammar.classify`](src/discern/grammar.py):

| Example String | Meaning | Flags Set | Role in Discern |
|---|---|---|---|
| `Rv0678@R50Q` | Missense substitution at amino acid `50` (`R` → `Q`) | `is_mis=True`, `pos=50` | Evaluated in Tier A & Tier B |
| `Rv0678@141_ins_c` | Insertion of cytosine at nucleotide `141` (frameshift) | `is_indel=True` | Evaluated in Tier A & Tier B (`rv_indel`) |
| `Rv0678@193_del_g` | Deletion of guanine at nucleotide `193` (frameshift) | `is_indel=True` | Evaluated in Tier A & Tier B (`rv_indel`) |
| `Rv0678@135_minorindel` | Minor subclonal indel (`MINOR_MUTATION`) | `is_indel=True` | Evaluated with its `FRSf` |
| `Rv0678@Q131!` | Nonsense premature stop codon at residue `131` | `is_stop=True` | Evaluated in Tier A & Tier B (`rv_stop`) |
| `Rv0678@c-11a` | Promoter point mutation `11` bp upstream of start | `is_prom=True` | Evaluated in Tier A & Tier B (`rv_prom`) |
| `Rv0678@S68S` | Synonymous (first char == last char, e.g., `S` → `S`) | **`IS_SYN=True`** | Ignored (silent change; excluded from Tier A/B & timelines) |
| `Rv0678@S68Z` | Uncallable / low-depth locus (`Z` suffix) | **`IS_NULL=True`** | Triggers `NULL_CALL` warning; excluded from calling |

**Why sequence-based indel canonicalization ([src/discern/indel_canon.py](src/discern/indel_canon.py)) is critical:**
Inside a poly-G track like `GGGGGG` spanning coding positions `193–198`, an extra base anywhere inside the run produces the **exact same physical DNA sequence**, even though different pipelines name it at opposite ends (CRyPTIC/GARC left end `192_ins_g` vs HGVS/TB-Profiler right end `c.198dupG`). `IndelIndex` applies the variant to the actual `498`-nt H37Rv `Rv0678` coding sequence (`RV0678_CDS`), indexes the resulting mutant DNA string, and resolves equivalent HGVS/GARC spellings (such as `c.139dupG`, `c.141_142insC`, `c.209delA`, `c.193_195del`) to the **exact matching training variant name** when their edited sequences are identical (falling back to left-alignment only for genuinely novel coding indels).

**Native Lab Importers ([src/discern/importers.py](src/discern/importers.py)):** Users do not have to hand-convert outputs from clinical sequencers. `detect_and_parse_upload()` auto-detects and converts:
1. Standard **4-column Discern CSV** (`sample_id,gene,mutation,frs`) and **6-column Longitudinal CSV** (`...,patient_id,collection_date`),
2. **TB-Profiler JSON** (`dr_variants` + `other_variants` with `hgvs_to_garc` normalization and VCF frequency units),
3. **GenoScreen Deeplex Myc-TB CSV** (percentage frequency auto-scaled by `/100`, even when a fractional percent like `0.45%` is below `1.0`),
4. **ONT AmPORE-TB / `wf-tb-amr` TSV & CSV exports** (auto-naming single-file uploads after the filename if no `sample_id` column is included).

---

## 5. Two-Tier Prediction Engine & Feature Engineering

### 5.1 Tier A: Curated Binomial Catalogue + Epistatic Veto ([src/discern/catalogue.py](src/discern/catalogue.py))

The catalogue is built on the training cohort using **`catomatic==0.2.3`** (`BinaryBuilder`) and interpreted by **`piezo==0.9.2`**:
- Filters applied before catalogue building: remove `IS_SYN`, remove `IS_NULL`, remove `mmpL5` & `mmpS5`, keep rows with **`FRS >= 0.1`**.
- Statistical rule per variant (`catomatic.build`): **two-tailed Binomial test** against an assumed background solo-resistance rate of `10%` (`background=0.1, tails="two", strict_unlock=True`) combined with wildcard rules (`bdq_wildcards_no_mmpL5.json`).
- Outputs a GARC-format lookup table (`models/lr_model/catalogue_bdq.csv`, `331` lines) with three wildcard defaults (`Rv0678@*=`, `pepQ@*=`, `atpE@*=` → `S` unless listed as `R` or `U`).
- **Sample-level aggregation (`R > U > S`)**:
  - Look up every non-synonymous, non-null core variant (`FRS >= 0.1`). Any unknown variant returning a dictionary error or `"U"` contributes `U`.
  - If **any** variant is `R` → **`R`**; else if **any** variant is `U` → **`U`**; else → **`S`**.
  - **Live Epistatic Veto (`decision.py`):** If the initial Tier A call is `R`, check whether any `atpE` or `pepQ` variant caused an `R`. If **`Rv0678` is the only source of `R`**, AND the sample carries a truncating `mmpL5` or `mmpS5` LoF variant present in the audited high-confidence lookup (`models/lr_mechanism/lof_effects.csv`, `veto_high_conf`) whose `FRS >= max_Rv0678_FRS - 0.20`, flip the final call to **`S`**. (And if Tier A evaluation raises a runtime exception, `CATALOGUE_EVALUATION_FAILED` safely forces a deferral to **`U`** rather than crashing or mis-calling.)

### 5.2 Tier B: The Exact 120-Feature Logistic Regression Model (`F2`)

In [src/discern/features.py](src/discern/features.py) and verified in [models/lr_model/metadata.json](models/lr_model/metadata.json), every isolate is turned into a vector of **`120` continuous numbers** bounded in `[0, 1]` using the maximum allele frequency (`FRSf >= 0.10`) observed for each feature:

1. **Feature Set `F0` — Core Variant Vocabulary (`111` features in production):**
   One column per specific `GENE@MUTATION` in `{Rv0678, atpE, pepQ}` that appeared in **at least `3` training isolates** (`vocab_min_count: 3`). Value = max `FRSf` (`0.0` if absent). *Because `vocab_from` is run strictly inside each training fold, rare singletons never leak.*
2. **Additional `F1` Features — Gene & Domain Max-FRS Aggregates (`+7` features):**
   Allow the model to capture general domain disruption even when a specific substitution is rare (`2` counts instead of `3`):
   - `rv_missense`: max `FRSf` across all `Rv0678` missense variants,
   - `rv_mis_dbd`: max `FRSf` across `Rv0678` missense variants inside the **DNA-binding domain** (residues `34 <= pos <= 99`),
   - `rv_indel`: max `FRSf` across `Rv0678` frameshift/in-frame indels,
   - `rv_stop`: max `FRSf` across `Rv0678` premature stop codons (`!`),
   - `rv_prom`: max `FRSf` across `Rv0678` promoter variants (`c-11a`, `t-8c`, etc.),
   - `atpE_any`: max `FRSf` across any non-synonymous `atpE` variant,
   - `pepQ_any`: max `FRSf` across any non-synonymous `pepQ` variant.
3. **Additional `F2` Features — Epistatic Efflux Aggregates (`+2` features → Total `120`):**
   - `mmpL5_LoF`: max `FRSf` over `mmpL5` indels or stop codons,
   - `rv_x_L5LoF`: multiplicative interaction **`rv_any * mmpL5_LoF`** (where intermediate `rv_any = max(rv_missense, rv_indel, rv_stop, rv_prom)` is computed inside `sample_tables`), giving the linear classifier a negative coefficient when both genes fail simultaneously!

**Why L2 Logistic Regression (`C=0.3`, NO `class_weight`) Won Head-to-Head Comparisons (`E05`, `E05b`, `E11`, `E21`):**
Under pre-registered criterion **D013** (`eval/e21_tuned_comparison.py`), a more complex classifier may only replace golden L2 Logistic Regression if it gains **`> +0.005` AUC on S1** AND keeps **external KM AUC >= golden (`0.8682`)**.
- **Feature ablation (`E05`)** proved the aggregate terms are essential: `F0` alone gets S1 AUC `0.791` (AUPRC `0.503`), `F1` jumps to `0.905` (`0.672`), and `F2` reaches `0.909` (`0.684`). Using `FRS >= 0.75` collapses `F2` back to `0.799` / `0.548`.
- **Balancing hurts (`E11`)**: `1:1` class balancing drops AUPRC (`0.684 -> 0.659`) and destroys baseline probability calibration.
- **Model comparison (`E05b` & `E21`)**: Tuned Random Forest (`0.9103` S1 AUC, `+0.0018`), Tuned HistGradientBoosting (`0.9058` S1 AUC, `0.8892` KM AUC), Tuned L1 LR (`0.9086`), and RBF SVM (`0.8985`) all failed to clear `+0.005` S1 improvement while preserving KM performance, so interpretable **L2 Logistic Regression** remains the production model.

### 5.3 Protein Domain Architecture & 3D Biophysics ([src/discern/domains.py](src/discern/domains.py), [src/discern/structure.py](src/discern/structure.py))

- **Domain Tracks ([src/discern/domains.py](src/discern/domains.py))**: Encodes exact verified H37Rv coordinates for all five loci (e.g., `Rv0678` promoter `-49..-1`, N-terminus/dimerization, Winged-helix DBD `34–99` with recognition helix `62–68`, C-terminal dimerization) and builds interactive inline SVG lollipop tracks rendered on `/result`.
- **3D Structural & Biophysical Exploration ([src/discern/structure.py](src/discern/structure.py), `E08` & `E22`)**: Parses raw coordinates from `Rv0678_WT.pdb` (chains A/B protein residues `15–163`, chains C/D DNA) to compute minimum heavy-atom distance to DNA (`min_dna_dist`) and 10 Å Cα burial (`max_burial`), plus Grantham physicochemical distance (`max_grantham`), charge delta (`max_charge`), recognition helix flag (`rec_helix`, `62–68`), and proline/glycine helix breakers (`helix_breaker`).
  - *Honest scientific finding*: On the strictly novel all-missense subset `E08` (`n=312`, S4 variant-set split), region+count got AUC `0.545`, adding 3D distance/burial reached `0.593` (baseline negative result), and exploratory biophysics reached `0.688`. However, when rigorously tested across all `11,868` isolates on S1 in [`eval/e22_biophysics_full_cohort.py`](eval/e22_biophysics_full_cohort.py), full-cohort S1 AUC was **`0.9058` vs baseline `0.9085`** — no improvement — so biophysical features remain excluded from production, confirming why Discern safely refers novel variants to culture!

---

## 6. Probability Calibration, Prevalence Adjustment & Decision Policy (Step-by-Step Math)

Raw outputs of regularized classifiers (`C=0.3`) and prevalence-skewed cohorts cannot be read directly as clinical risk. Discern runs a three-step mathematical pipeline ([src/discern/model.py](src/discern/model.py), [src/discern/prevalence.py](src/discern/prevalence.py)).

### Step 1: Nested `GroupKFold` Platt Scaling (`score_raw` → `p_train`)

To prevent calibration overfitting, [`fit_nested`](src/discern/model.py) runs a **5-fold patient-grouped cross-validation *inside* the training set**:
1. Collect honest out-of-fold raw scores `raw_oof` (clipped to `[1e-6, 1 - 1e-6]`),
2. Convert to log-odds: `z = logit(raw) = ln(raw / (1 - raw))`,
3. Fit an unregularized 1-D Logistic Regression (`C=1e6`) mapping `z -> y`, producing intercept `b_cal` and slope `w_cal`,
4. Fit the main `C=0.3` model on all training isolates (`11,868`). At prediction time, a raw probability `score_raw` becomes `p_train = sigmoid(b_cal + w_cal * logit(score_raw))`.

*Concrete check from the live model:* If `score_raw = 0.840964`, `logit(0.840964) = ln(0.840964 / 0.159036) = +1.6654`. Passing through the fitted Platt scaler stretches confident log-odds to **`p_train = 0.900611` (`90.1%`)** at the training cohort's baseline prevalence **`pi_tr = 866 / 11868 = 0.072969` (`7.30%`)**.

### Step 2: Odds-Ratio Prevalence Correction (`p_train` → `p_adj`)

In clinical deployment, local resistance prevalence `pi` rarely equals `7.30%` (it may be `1%` in low-burden regions or `5%` default). Using Bayes' rule in odds form ([`prevalence.adjust`](src/discern/prevalence.py)):

```text
odds(p)  = p / (1 - p)
odds_adj = odds(p_train) / odds(pi_tr) * odds(pi)
p_adj    = odds_adj / (1 + odds_adj)
```

**Worked Numerical Example (`Rv0678@141_ins_c`, user clinic `pi = 0.05`):**
1. Training prevalence odds: `odds(0.072969) = 0.072969 / 0.927031 = 0.078713`
2. Target clinic odds: `odds(0.05) = 0.05 / 0.95 = 0.052632`
3. Model odds at training prevalence: `odds(0.900611) = 0.900611 / 0.099389 = 9.06148`
4. Multiply by the prior-odds shift factor (`0.052632 / 0.078713 = 0.66865`):
   `odds_adj = 9.06148 * 0.66865 = 6.05897`
5. Convert odds back to probability:
   **`p_adj = 6.05897 / (1 + 6.05897) = 0.858336` (`85.8%`)** — matching the live `Predictor` output to 6 decimal places!

On cross-source prevalence shock (`E10`), this zero-parameter adjustment collapses Expected Calibration Error (**ECE**) from **`0.288 -> 0.013`** (`NICD -> CRyPTIC`) and from **`0.568 -> 0.101`** (`CRyPTIC -> NICD`). Furthermore, [src/discern/ui/static/app.js](src/discern/ui/static/app.js) evaluates this exact formula live in the browser as you drag the prevalence slider.

### Step 3: Bayes Positive & Negative Predictive Values (`PPV` & `NPV`)

Predictive values shown on `/result` use **cross-validated ternary catalogue call rates** from `e16_lof_ablation` (`sens_firm = 0.7918`, `spec_firm = 0.9872`, reweighted across all archived isolates via `archived_predictive_values()` in [src/discern/prevalence.py](src/discern/prevalence.py) using `sens_U_as_S` / `spec_U_as_S` for PPV and `spec_U_as_R` / `sens_U_as_R` for NPV):
- **At `pi = 5%`:** **`ppv_if_R = 72.1%`** (`0.7207`) and **`npv_if_S = 99.1%`** (`0.9912`).
- **At `pi = 2%`** (`E13`, standard catalogue): **`50.2%`** PPV — meaning half of firm `R` calls in a `2%`-prevalence setting can be false alarms purely from Bayes' law, underscoring why laboratory confirmation remains critical.

### Step 4: Unsupervised Batch Prevalence Estimation via MLLS ([src/discern/label_shift.py](src/discern/label_shift.py))

When a lab uploads a batch CSV onto `/queue` (`>= 20` unique deduplicated patients), it often does not know its own local resistance rate `pi`. `estimate_prevalence_mlls()` runs **Saerens Expectation-Maximization (Maximum Likelihood Label Shift)** using the calibrated training probabilities `p_train` (deduplicated to each patient's earliest chronological isolate) to estimate `pi_hat` plus a `200`-replicate bootstrap 95% confidence interval (`n_bootstrap=200`) — with **zero true lab labels required** — and lets the user apply that estimated prevalence to the batch PPV/NPV curve in one click.

### Step 5: Calibration Thresholds & Complete Warning Code Matrix

Production model cutoffs frozen in [models/lr_model/metadata.json](models/lr_model/metadata.json):
- **`t_bin = 0.0230`**: Youden index (`sens + spec - 1`) optimal binary operating threshold.
- **`lo = 0.0549`**: Largest threshold where inner OOF isolates below `lo` (`n >= 20`) achieve `NPV >= 0.985` (`98.5%`).
- **`hi = 0.7273`**: Smallest threshold where inner OOF isolates above `hi` (`n >= 20`) achieve `PPV >= 0.80` (`80.0%`).

| Warning Code | Trigger Condition | Why It Matters |
|---|---|---|
| `MINOR_ALLELE_ONLY` | All contributing variants have `FRS < 0.90` | Resistance/susceptibility depends solely on non-fixed subclonal reads. |
| `HETERORESISTANCE_SUBCLONAL` | Any non-syn/non-null core variant has `0.10 <= FRS < 0.75` | Mixed population; can cause discordance between fast plates (UKMYC) and slow MGIT960. |
| `NOVEL_VARIANT` | Variant absent from `evidence_store.csv` (`0` prior carriers) | Never seen in training cohorts (`Rv0678@Q131!` → forces `U`). |
| `LOW_EVIDENCE` | Variant has `< 10` total carriers in `evidence_store.csv` | Small sample count (`Rv0678@R134!` has high point estimate, wide Wilson CI). |
| `NULL_CALL` | Any input ends in `Z` (e.g., `S68Z`) | Uncallable base position; re-sequence if critical region is missing. |
| `PUMP_DAMAGE_MODIFIER` | `mmpL5_LoF > 0` co-occurs with any non-syn/non-null `Rv0678` variant | Pump-loss-of-function may alter or nullify `Rv0678` efflux resistance. |
| `EPISTATIC_PUMP_INACTIVATION` | Live Epistatic Veto fires (audited high-conf `mmpL5`/`mmpS5` LoF matching `Rv0678` FRS) | Explains why a Tier A `Rv0678` `R` variant was overridden to `S`. |
| `NO_VARIANT_INVISIBLE_RESISTANCE_POSSIBLE` | Final call is `S` with zero non-syn/non-null variants | Reminds clinician that `~15.8%` of BDQ resistance lies outside these 5 loci. |
| `MODEL_DISAGREES` | Tier A says `S` while `p_train >= hi` (`0.7273`), or Tier A says `R` while `p_train <= lo` (`0.0549`) | Flags catalogue vs logistic-model disagreement for human inspection. |
| `CATALOGUE_EVALUATION_FAILED` | `piezo` raises an unexpected error during rule evaluation | Safe fail-closed behavior forcing `U` deferral instead of silent crash. |
| `EPISTATIC_VETO_BLIND` | Selected assay (`deeplex_myc_tb`, `ampore_tb`, `icmr_nirt_minion`) does not cover `mmpL5`/`mmpS5` | Alerts lab that epistatic pump suppression cannot be detected by that panel. |
| `UNSEQUENCED_CORE_LOCI` | Selected assay misses core resistance genes (`atpE` and/or `pepQ`) | Alerts lab to uncovered resistance loci (`E24` blind spots). |
| `HOMOPOLYMER_ARTIFACT_RISK` | Nanopore platform + indel inside a `G`/`C` homopolymer stretch (e.g., `141_ins_c`) | Guards against false Nanopore insertion/deletion basecalling errors. |
| `PANEL_BELOW_LOD` | Input `FRS` is below the selected assay's validated VAF limit of detection (`vaf_floor`) | Variant sits near/below technical sequencing noise floor for that panel. |

---

## 7. Clinical Workbench Features & Privacy Boundaries

### 7.1 Upload Quality, Longitudinal Tracking, and Batch Contamination Surveillance

- **Fresh Upload Quality Check ([src/discern/quality.py](src/discern/quality.py))**: Evaluates CSV/native uploads independently of resistance inference, scoring completeness across required columns, flagging duplicate `(sample_id, gene, mutation)` rows, malformed FRS (`<=0` or `>1`), synonymous rows, and bad optional dates — displayed under **"Upload quality"** so data hygiene never masquerades as clinical confidence.
- **Longitudinal Serial-Isolate Tracking ([src/discern/longitudinal.py](src/discern/longitudinal.py))**: When a CSV includes optional `patient_id` and `collection_date` columns (e.g., `Day 0`, `Day 28`, `Day 84` in [demo/longitudinal_example.csv](demo/longitudinal_example.csv)), `analyze_patient_timelines()` links visits chronologically, resolves overlapping variants by taking maximum FRS, tracks `DE_NOVO_EMERGENCE`, `SUBCLONAL_EXPANSION`, and `POLYCLONAL_COMPETITION`, and renders an interactive Muller-style clonal-ribbon SVG chart with playback scrubber on `/patient/{patient_id}`.
- **Batch Recurrence & Contamination Audit ([src/discern/surveillance.py](src/discern/surveillance.py))**: Inspects uploaded batches after deduplicating by `patient_id` so one patient's serial visits don't look like an outbreak. If `>= 2` **distinct patients** in a batch carry the same homopolymer indel, it raises a **QC contamination warning**; if `>= 2` distinct patients carry an identical rare non-homopolymer variant, it recommends **WGS outbreak / transmission investigation**.

### 7.2 Privacy, Isolation, and Error Safety ([src/discern/api/app.py](src/discern/api/app.py), [src/discern/research.py](src/discern/research.py))

1. **Public Replay Archive vs. Private Session Store**:
   - The `Archive` class loads **frozen, immutable** evaluation artifacts (`demo/held_out_samples.csv`, `results/oof_predictions.parquet`) and serves public replays labeled `"held-out in cross-validation"` so recorded examples never display overfit training predictions.
   - Live CSV uploads enter an ephemeral per-browser `UploadedBatch` scoped by a random `HttpOnly`, `SameSite=Strict` session cookie (`discern_upload`), served with `Cache-Control: no-store`, and strictly excluded from `/replay/{sample_id}`.
2. **Zero Payload Echo & Zero Access Logs**: `Makefile` and `serve.ps1` start Uvicorn with `--no-access-log`. CSV parser errors (`UploadError`) report exact row/column locations without echoing back raw user-submitted text (`test_backend_hardening.py`). JSON endpoints reject free-text genes/mutations with HTTP `422` via Pydantic strict schema validation (`schemas.py`).

---

## 8. Two Complete Live Walkthroughs (Verified Predictor Outputs)

### Walkthrough A: Clear Resistance (`Rv0678@141_ins_c`, `frs=0.92`, `pi=0.05`, `assay='ampore_tb'`)

```text
Input: [{"gene": "Rv0678", "mutation": "141_ins_c", "frs": 0.92}], prevalence=0.05, assay="ampore_tb"
```

1. **Validation & Grammar**: Pydantic accepts `gene="Rv0678"` and `grammar.valid_mutation("141_ins_c") == True`. No `IS_SYN` or `IS_NULL`. `FRS = 0.92 >= 0.10`.
2. **Tier A (`catalogue.py` + `decision.py`)**: Matches `Rv0678@141_ins_c` → **`R`** (`evidence`: `234` carriers, `187` resistant [`79.9%`], split `CRyPTIC 11/55` vs `NICD 176/179`). No `mmpL5`/`mmpS5` LoF is present, so Epistatic Veto cannot fire. Final **`call = "R"` (`"Resistant"`)**.
3. **Tier B (`features.py` + `model.py`)**: Sets `Rv0678@141_ins_c = 0.92`, `rv_indel = 0.92`, `rv_any = 0.92`, other 117 features `0.0`. Produces:
   - `score_raw` = **`0.840964`**
   - `probability_train_prevalence` (at `7.30%`) = **`0.900611`** (`90.1%`)
   - `probability_at_prevalence` (shifted to `5.0%`) = **`0.858336`** (`85.8%`)
   - `queue_rank_percentile` = **`98.4871`** (98.5th percentile across archived OOF scores)
4. **Predictive Values & Panel Guardrails**:
   - At `5%` prevalence using audited catalogue metrics (`sens=0.7918`, `spec=0.9872`): **`ppv_if_R = 0.720710` (`72.1%`)**, **`npv_if_S = 0.991246` (`99.1%`)**.
   - `panels.py` inspects `assay="ampore_tb"` (`Nanopore`, covers only `Rv0678, atpE`) + notices `141_ins_c` sits in a poly-C homopolymer run (`>= 3` identical bases), appending **`["EPISTATIC_VETO_BLIND", "UNSEQUENCED_CORE_LOCI", "HOMOPOLYMER_ARTIFACT_RISK"]`**.

### Walkthrough B: Novel Referral (`Rv0678@Q131!`, `frs=0.92`) vs. Low-Evidence Call (`Rv0678@R134!`, `frs=0.92`)

```text
Input B1: [{"gene": "Rv0678", "mutation": "Q131!", "frs": 0.92}], prevalence=0.05, assay="ampore_tb"
```

1. **Tier A**: `Rv0678@Q131!` is absent from `catalogue_bdq.csv` rules → returns **`U` ("Send for lab test")**. Tier B does **not** override it.
2. **Tier B**: Although `Rv0678@Q131!` is missing from the `111`-variant vocabulary (`NOVEL_VARIANT`), its mutation suffix ends in `!` (`is_stop=True`), so the aggregate column **`rv_stop = 0.92`** (alongside internal `rv_any = 0.92`) automatically fires!
   - `score_raw` = **`0.345453`**, `probability_train_prevalence` = **`0.393114`** (`39.3%`), `probability_at_prevalence` = **`0.302223`** (`30.2%`), ranking at **`60.6657`** percentile (`queue_rank_percentile`).
   - Warnings emitted: **`["NOVEL_VARIANT", "EPISTATIC_VETO_BLIND", "UNSEQUENCED_CORE_LOCI"]`**.
3. **Contrast with `Rv0678@R134!`** (`frs=0.92, pi=0.05, ampore_tb`): Unlike `Q131!`, `R134!` appeared in training just often enough to enter the catalogue as **`R`**, yet with fewer than `10` total carriers (`LOW_EVIDENCE`). Its Tier B score is `score_raw=0.499383`, `p_train=0.573117`, `p_adj=0.473049`, percentile `77.9123`, warnings `["LOW_EVIDENCE", "EPISTATIC_VETO_BLIND", "UNSEQUENCED_CORE_LOCI"]`. Notice how Discern calls `R134!` conservatively with `LOW_EVIDENCE` while routing `Q131!` straight to culture!

---

## 9. Complete Evaluation Ledger & Honest Negative Results

### 9.1 Understanding `14/14` vs `14/49`

- Running **`eval/verify.py`** (or `make verify`) checks the **14 pre-registered acceptance records** that have numerical tolerances defined in `expected` / `tolerance` inside `results/index.json` → **`14/14 acceptance tests passed`**.
- Reading **`results/INDEX.md`** displays **`14/49 tests pass`** across *all* 49 indexed folders. Why aren't all 49 green? Because the remaining 35 entries consist of (a) descriptive/audit records with `pass_all: null` (such as `E05b`, `E14`, `E15`, `E16`, `E17`, `E19`, `E21`, `G7`, `G8`), and (b) negative scientific findings or blocked external downloads (`E18` expanded CRyPTIC v3.4 degradation, `Sobol`/`Choi` ENA access blockers) preserved intentionally for scientific transparency!

### 9.2 Core Pre-Registered Benchmarks (`E01–E13`)

| Exp | Benchmark Description | Key Observed Result(s) | Pass / Role |
|---|---|---|---|
| **E01** | Data prep exact counts (`eval/e01_prep.py`) | `11,868` samples (`866` R), CRyPTIC `10,704`/`83 R`, NICD `1,164`/`783 R`, `11,611` SUBJID, `25,554` mut rows | **PASS** (Exact) |
| **E02** | WHO v2 catalogue as-is (`eval/e02_who.py`) | `FRS>=0.75`: firm `0.968`, sens_firm `0.561`, spec `0.987`, `U-as-S sens 0.497`. `FRS>=0.1`: firm `0.961`, sens_firm `0.801`, spec `0.981`, `U-as-S sens 0.679` | **PASS** (Baseline) |
| **E03** | Tier A S1 patient-grouped CV (`eval/e03_catalogue_cv.py`) | Firm `0.9443`, sens_firm `0.7924` (CI `0.760–0.824`), spec_firm `0.9835`, `661 U` (`206 R`); `U-as-S sens 0.6039`, `U-as-R sens 0.8418 / spec 0.9428` | **PASS** (Gate G1) |
| **E04** | Tier B L2 LR `F2` on S1 nested CV (`eval/e04_lr_f2.py`) | **AUC `0.9085`** (CI `0.894–0.923`), **AUPRC `0.6836`** (`0.647–0.719`), bin sens `0.8349`, spec `0.9481`, ternary firm `0.9469` (`0.735`/`0.989`), sens@spec95/97/98/98.5 = `0.830 / 0.782 / 0.652 / 0.604` | **PASS** (Gate G2) |
| **E05** | Feature ablation on S1 (`eval/e05_ablation.py`) | AUC/AUPRC: `F0` `0.791/0.503` → `F1` `0.905/0.672` → `F2` `0.909/0.684`; `F2@FRS0.75` `0.799/0.548` | **PASS** |
| **E06** | Cross-source S3 (`eval/e06_cross_source.py`) | `CRyPTIC -> NICD` AUC `0.9105`, sens `0.8442`, spec `0.9633`; `NICD -> CRyPTIC` AUC `0.8277`, sens `0.6747`, spec `0.9532` | **PASS** |
| **E07** | Leave-One-Site-Out S2 (`eval/e07_loso.py`) | Pooled binary sens `0.7875`, spec `0.9537`; site spec (`n>=300`) `0.9441–0.9785` (`per_site.csv`: `11` sites; Mumbai `950/5R`, Peru `3,035/5R`). No pooled AUC reported. | **PASS** |
| **E08** | Novel-variant S4 (`eval/e08_novel_variant.py`) | `n=312` (`143` variant sets, `169 R`): region+count AUC `0.5451`; `+structure` (`min_dna_dist`, `max_burial`) **`0.5932`** (negative baseline); exploratory `+biophysics` `0.6876` | **PASS** (Gate: abstains on new muts) |
| **E09** | Hybrid exploration inside `661 U` (`eval/e09_hybrid.py`) | ML AUC inside `U`: `0.7442`; resolves `173` `S` (`17 R`), `126` `R` (`51 S`), leaves `362` `U` (`114 R`). Exploratory only — never auto-called. | **PASS** |
| **E10** | Prevalence shift calibration (`eval/e10_prevalence_shift.py`) | `NICD -> CRyPTIC`: ECE **`0.288 -> 0.013`**, Brier `0.102 -> 0.017`; `CRyPTIC -> NICD`: ECE **`0.568 -> 0.101`**, Brier `0.501 -> 0.124` | **PASS** |
| **E11** | Class imbalance & learning curve (`eval/e11_imbalance.py`) | `1:1` balanced AUC `0.9051`, AUPRC `0.6587` (worse than uncalibrated LR `0.6836`); fraction `0.10 -> 1.00`: AUC `0.8963 -> 0.9085` | **PASS** |
| **E12** | Priority list inside U (`eval/e12_priority.py`) | `332` distinct U variant sets; top-5 cover `27.7%`, top-10 `34.3%`, top-20 `39.8%`, top-50 `49.6%`; largest: `Rv0678@c-11a` (`102`, `8.8% R`) | **PASS** |
| **E13** | Bayes PPV/NPV matrix (`eval/e13_ppv_npv.py`) | At `1/2/5/10%` prevalence: Cat PPV `33.3/50.2/72.2/84.6%`, Model PPV `14.4/25.3/46.6/64.8%`, NPV `97.7%–99.8%` | **PASS** |

### 9.3 Extended Evaluations, Ship Decisions & Negative Audits (`E05b`, `E14–E24`, `G7–G8`, `KM`)

| Study ID | What Was Tested | Honest Outcome & Decision |
|---|---|---|
| **E05b / E21** | Model-class comparison (`lr_l1`, `rf`, `hgb`, `svm`, raw + tuned) on identical S1 & KM | **Kept Golden L2 LR**: `rf_tuned` (`0.9103` S1, `0.8675` KM) and `hgb_tuned` (`0.9058` S1, `0.8892` KM) failed D013 (`>0.005` S1 AUC gain + preserve KM). |
| **E14** | Learn-then-Test conformal abstention (`conformal.py`) & Beta-Binomial shrinkage (`shrinkage.py`) | **Rejected conformal thresholds** (missed `<= 0.015` S-error target: `best=0.0155`, `CRyPTIC->NICD` S-error `0.671`). Shrinkage pulled `1,042` variants (`48` flipped a `0.5` grade; descriptive only). |
| **E15** | Beta entropy Value-of-Information queue (`voi.py`) | `info_gain_minus_count_K10 = 0.0` (ranking by expected information gain ties simple carrier count). Kept offline. |
| **E16 / E19** | `Rv0678` + `mmpL5/S5` epistatic LoF & D011 Ship Rule (`e16_mechanism.py`, `e19_ship_decision.py`) | Snobre co-occurrence replication: **`25.04%` (`162/647`)** vs `24.8%` published. On S1, veto flipped `41` calls (`39 S`, `2 R`), raising catalogue `spec` to **`0.9872`** (`sens` `0.7918`). **Shipped** audited high-conf lookup veto in `decision.py`. |
| **E17** | Borderline MIC label-noise audit (`eval/e17_mic_audit.py`) | Found **`454` isolates (`191 R`)** right on critical concentration boundaries, quantifying inherent lab assay noise. |
| **E18 / E18b** | Training on expanded CRyPTIC v3.4 (`23,441` isolates, or `16,923` HIGH-only) | **Negative result**: S1 AUC degraded to `0.8381` (all) and `0.6668` (HIGH-only) due to cross-protocol heterogeneity; `models/lr_model` stays on clean `11,868` set. |
| **E22** | Full-cohort evaluation of `structure.py` biophysical features (`eval/e22_biophysics_full_cohort.py`) | S1 AUC was **`0.90578`** vs baseline `0.9085` — biophysics helps a small novel missense subset (`E08`) but does not beat vocabulary+F2 globally. Excluded from live model. |
| **E23** | Experimental standalone pipeline (`eval/e23_novel_pipeline.py`) | Combines `mmpL5` rare lineage filtering, L2/L4 lineage indicators, and continuous `log2(MIC)` regression; reaches KM external AUC **`0.89355`** in standalone research runs (unindexed, outside production). |
| **E24** | Sequencing panel blind-spot audit (`eval/panel_audit.py`) | Computes coverage ceilings (`84.18%` WGS/Deeplex XL, `79.56%` AmPORE-TB, `74.48%` Deeplex/ICMR-NIRT) across the `866` resistant isolates (`137` invisible = `15.8%`). |
| **Ext-KM** | Held-out CRyPTIC v3.4.0 external cohort (`792` isolates, `35 R`, `eval/ext_km_validate.py`) | Catalogue firm rate `0.7753` (`178 U`, `11 R`), `cat_sens = 0.7917`, `cat_spec = 0.9881`, `model_auc = 0.8682`. |
| **Ext-Sobol / Choi** | Independent SA (`744`) and Korea (`49`) cohorts (`eval/ext_validate.py`) | Documented accession inventories; **honestly flagged as blocked** by offline/ENA FASTQ retrieval access. |
| **G7 / G8** | Adlard/WHO v2 head-to-heads + Uncertainty Context (`e18_head_to_head.py`, `g8_context.py`) | Confirms Discern does not claim to beat WHO at matched FRS (`0.768` vs `0.735` ternary) and publishes exact Wilson intervals across evaluation cells. |

---

## 10. Complete File-by-File Inventory of the Repository

Every file across the entire repository is listed below so nothing in `Discern` is a mystery.

### 10.1 Root Directory (`16` files)

| File | Purpose |
|---|---|
| [AGENTS.md](AGENTS.md) | Master build specification and hard scientific constraints (patient grouping, zero leakage, exact tolerances, offline guarantees). |
| [README.md](README.md) | High-level project overview, verified numbers from `results/`, architecture Mermaid diagram, and run instructions. |
| [CODEBASE_GUIDE.md](CODEBASE_GUIDE.md) | This comprehensive guide covering concepts, math, walkthroughs, and all repository files. |
| [DESIGN.md](DESIGN.md) | UI/UX design specification (light neutral tokens, offline system fonts, semantic HTML tables, accessible contrast, timeline/motion design). |
| [CHANGELOG.md](CHANGELOG.md) | Chronological log of engineering milestones, bug fixes, clinical guardrails, importers, and evaluation updates. |
| [config.yaml](config.yaml) | Global configuration: random seeds (`1, 2, 0`), paths, S1 `n_splits: 5`, `frs_min: 0.1`, `vocab_min_count: 3`, `dbd_range: [34, 99]`, `C: 0.3`, `npv_target: 0.985`, `ppv_target: 0.80`, `default_prevalence: 0.05`, and offline `enhanced` flags. |
| [pyproject.toml](pyproject.toml) | PEP 621 package build metadata for `discern` (`0.1.0`), registering the `discern` console script, UI package-data globs, and pytest settings. |
| [requirements.txt](requirements.txt) | Runtime and dev dependencies (`pandas>=2.2`, `numpy>=1.26`, `scipy>=1.11`, `scikit-learn==1.9.1`, `pyarrow>=15`, `joblib>=1.3`, `catomatic==0.2.3`, `piezo==0.9.2`, `fastapi>=0.110`, `uvicorn>=0.29`, `pydantic>=2.6`, `jinja2>=3.1`, `pyyaml>=6`, `pytest>=8`, `httpx>=0.27`, `ruff>=0.4`). |
| [Makefile](Makefile) | Convenience commands: `make verify` (acceptance checks + pytest), `make test`, `make docs` (regenerate `MODEL_CARD.md`/`RESULTS.md`), and `make app` (`uvicorn` on port `8000` with `--no-access-log`). |
| [Dockerfile](Dockerfile) | Offline `python:3.12-slim` image copying `src/`, `config.yaml`, `demo/`, `models/lr_model/`, `models/lr_mechanism/lof_effects.csv`, `results/`, and `docs/`, running Uvicorn as non-root user `discern`. |
| [serve.ps1](serve.ps1) | Windows PowerShell one-click launcher starting `discern.api.app:app` at `http://127.0.0.1:8000` with `--no-access-log`. |
| [.gitignore](.gitignore) | Prevents committing `.venv/`, `data/raw/`, `data/processed/`, `data/external_cryptic_v34/`, `*.tar.gz`, `*.mp4`, `*.zip`, `build_discern_demo.py`, and caches. |
| [.dockerignore](.dockerignore) | Keeps raw data, dev caches, zips, and unneeded model dirs out of Docker build contexts while explicitly whitelisting `models/lr_model/*` and `models/lr_mechanism/lof_effects.csv`. |
| [.gitattributes](.gitattributes) | Enforces LF line endings on key templates, Python server modules, JS, and `MODEL_CARD.md`. |
| `build_discern_demo.py` | Local demo video composition script (`1,807` lines, git-ignored; uses FFmpeg/PIL/NumPy to edit/render 4K60 screen capture demos with lower-third captions). |
| `discern.zip` | Local backup archive (`3.3 MB`, git-ignored) containing a portable packaged snapshot. |

### 10.2 Core Library: `src/discern/` (`29` modules + `api/` subpackage)

#### A. Data Prep, Parsing & Normalization
- [src/discern/__init__.py](src/discern/__init__.py): Package initializer (`Discern: mechanism-aware, externally-validated bedaquiline resistance triage`).
- [src/discern/prep.py](src/discern/prep.py): Builds clean `samples.parquet` and `mutations.parquet` from raw CRyPTIC (`QUALITY == "HIGH"`) and NICD (`MGIT960`) tables, joins patient IDs (`SUBJID`), builds `GENE@MUT`, and validates E01 exact counts (`load_samples`, `load_mutations`, `prepare`, `e01_counts`).
- [src/discern/splits.py](src/discern/splits.py): Deterministic cross-validation split generators preventing patient/clonotype leakage (`s1_folds`, `loso_folds`, `cross_source_folds`, `variant_set_keys`).
- [src/discern/grammar.py](src/discern/grammar.py): Regex-validated GARC mutation grammar checker and flag classifier (`valid_mutation`, `suffix_of`, `classify` → `IS_SYN`, `IS_NULL`, `is_stop`, `is_indel`, `is_prom`, `is_mis`, `pos`).
- [src/discern/indel_canon.py](src/discern/indel_canon.py): Sequence-based indel canonicalizer matching HGVS and shifted GARC homopolymer coordinates directly against the H37Rv `Rv0678` reference ORF (`apply_garc`, `del_matches_reference`, `hgvs_indel_to_garc`, `left_align`, `IndelIndex`).
- [src/discern/importers.py](src/discern/importers.py): Auto-detects uploaded lab formats (4/6-col CSV, TB-Profiler JSON, Deeplex Myc-TB CSV, AmPORE-TB TSV/CSV) and converts HGVS mutations and frequency units into standardized Discern inputs (`hgvs_to_garc`, `detect_and_parse_upload`).
- [src/discern/quality.py](src/discern/quality.py): Fresh CSV hygiene auditor scoring completeness, missing/out-of-range FRS, duplicate variants, synonymous entries, and invalid dates separately from clinical inference (`check_csv_quality`).

#### B. Biology, Features & Domain Visualizations
- [src/discern/features.py](src/discern/features.py): Extracts sparse core variant tables (`vocab_from`), builds `F0`/`F1`/`F2` design matrices (`120` production features: `111` vocabulary variants + `9` biological max-FRS & `mmpL5_LoF` aggregates) (`sample_tables`, `design`).
- [src/discern/mechanism.py](src/discern/mechanism.py): Implements `Rv0678` DNA translation (`165` aa), protein domain classification (`DBD <= 99`, `mid 100–160`, `late > 160`), poly-G/C `_is_homopolymer` detection, frame-restoring vs truncating logic, and Snobre/WHO `FRS >= Rv - 0.20` epistatic linkage (`load_cds`, `translate`, `first_stop_aa`, `region`, `apply_indel`, `classify_truncating`, `variant_effects`, `veto_flags`, `sample_mech_flags`, `apply_veto`).
- [src/discern/domains.py](src/discern/domains.py): Defines verified H37Rv protein structural regions across all five loci (`DomainSegment`, `GeneArchitecture`) and generates SVG mutation lollipop tracks for `/result` (`build_domain_track`, `attach_domain_tracks`).
- [src/discern/structure.py](src/discern/structure.py): Parses `Rv0678_WT.pdb` coordinates (protein chains A/B vs DNA chains C/D), computes minimum residue-to-DNA distances, 10 Å Cα burial, Grantham chemical distance, charge deltas, α4 recognition-helix flags, and helix breakers (`grantham_distance`, `charge_delta`, `parse_pdb`, `sample_features`, `sample_biophysical_features`).
- [src/discern/panels.py](src/discern/panels.py): Defines tNGS/WGS assay profiles (`PanelProfile`: `wgs`, `deeplex_myc_tb`, `deeplex_xl`, `ampore_tb`, `icmr_nirt_minion`), emits blind-spot / below-LOD / homopolymer warnings (`audit_sample_panel`), and computes known-locus sensitivity ceilings (`compute_cohort_panel_coverage`).

#### C. Tier A Catalogue, Tier B ML, Prevalence & Decision Engine
- [src/discern/catalogue.py](src/discern/catalogue.py): Builds two-tailed Binomial rule catalogues via `catomatic` + `wildcards` (`cat_muts`, `build_catalogue`), evaluates isolates (`apply_catalogue`: `R > U > S`), and scores official WHO v2 GARC rules (`who_as_is`).
- [src/discern/model.py](src/discern/model.py): Golden Tier B estimator (`LogisticRegression(C=0.3)`), nested inner `GroupKFold` Platt logit calibration (`fit_nested`, `logit`, `score`), probability cutoff finder (`thresholds` → `t_bin`, `lo`, `hi`), and operating-point helper (`sens_at_spec`).
- [src/discern/runner.py](src/discern/runner.py): Shared nested out-of-fold execution loop reusable across arbitrary splits (`S1`, `S2`, `S3`) and arbitrary sklearn estimators (`nested_oof`, `fit_nested_est`, `score_est`, `nested_oof_est`).
- [src/discern/prevalence.py](src/discern/prevalence.py): Log-odds prevalence shifting (`odds`, `adjust`, `prevalence_weights`), Brier score & 10-bin ECE (`brier`, `ece`), and ternary Bayes PPV/NPV formulas (`ppv_npv`, `archived_predictive_values`).
- [src/discern/label_shift.py](src/discern/label_shift.py): Unsupervised batch prevalence estimation via Maximum Likelihood Label Shift (Saerens EM) with chronological patient deduplication and bootstrap 95% CIs (`estimate_prevalence_mlls`, `deduplicate_patient_probabilities`).
- [src/discern/evidence.py](src/discern/evidence.py): Serves per-variant evidence cards from `models/lr_model/evidence_store.csv`, pairing Wilson confidence intervals with split `CRyPTIC` vs `NICD` carrier counts (`EvidenceStore`).
- [src/discern/decision.py](src/discern/decision.py): Production inference orchestrator (`Predictor`) combining Tier A `piezo` calls, audited high-conf epistatic veto (`lof_effects.csv`), Tier B `120`-feature probability & percentile rank, Bayes PPV/NPV, `EvidenceStore`, panel blind-spot audit, and all clinical warning codes.

#### D. Longitudinal, Surveillance, Reporting & CLI
- [src/discern/longitudinal.py](src/discern/longitudinal.py): Aggregates dated multi-visit samples per patient (`parse_timepoint_date`), detects clonal replacement/expansion patterns (`analyze_patient_timelines`), and supplies data for the `/patient/{patient_id}` Muller ribbon UI.
- [src/discern/surveillance.py](src/discern/surveillance.py): Audits uploaded multi-sample batches (deduplicated by `patient_id`) to separate multi-patient homopolymer lab contamination alerts from rare-variant WGS transmission escalations (`batch_recurrence_audit`).
- [src/discern/metrics.py](src/discern/metrics.py): Computes firm coverage, sensitivity/specificity among firm calls, `U-as-S`, and `U-as-R` metrics across experiments (`ternary_metrics`).
- [src/discern/research.py](src/discern/research.py): Read-only loader for frozen JSON/CSV reports (`Reports`), public held-out replay cases (`Archive`, backed by `oof_predictions.parquet` & `demo/held_out_samples.csv`), CSV error validation (`UploadError`), and local file helpers (`workspace`, `read_json`, `csv_text`, `read_batch`).
- [src/discern/results_io.py](src/discern/results_io.py): Standardized experiment writer (`write_result`) generating `results.json`, `REPORT.md`, approach `SUMMARY.md`/`MODEL_REPORT.md`, and `results/INDEX.md` + `results/index.json` (`render_report`, `render_approach_summary`, `render_index`, `regenerate_all`, `migrate`).
- [src/discern/documentation.py](src/discern/documentation.py): Builds structured `docs/MODEL_CARD.md` sections dynamically from registered evaluation reports (`model_card`).
- [src/discern/cli.py](src/discern/cli.py): Entrypoint for the `discern` CLI supporting four subcommands: `replay --sample-id`, `batch <input> <output>`, `eval --all`, and `verify`. *(Note: CLI takes input via `replay` or `batch` CSV; single-sample JSON prediction is served via `POST /predict`.)*

#### E. Preserved / Exploratory Statistical R&D Modules (Offline Only)
- [src/discern/conformal.py](src/discern/conformal.py): Learn-then-Test / Bentkus p-value risk-controlled threshold auditor (`bentkus_pvalue`, `certified_upper_bound`, `certified_thresholds`, `apply_policy`). Audited in `E14` and left inactive after missing the strict S-error cap.
- [src/discern/shrinkage.py](src/discern/shrinkage.py): Empirical-Bayes Beta-Binomial shrinkage pulling noisy single-carrier variant proportions toward their biological mechanism-class prior (`wilson`, `beta_prior_moments`, `variant_shrinkage`). Audited descriptively in `E14`.
- [src/discern/voi.py](src/discern/voi.py): Beta entropy information-gain queue ranker estimating how much testing one `U` isolate resolves ambiguous variants (`beta_entropy`, `expected_info_gain`, `variant_unlock`, `rank_sets`). Found to tie carrier counts in `E15`.

#### F. FastAPI Web Application Subpackage: `src/discern/api/` (`3` files)
- [src/discern/api/__init__.py](src/discern/api/__init__.py): Subpackage marker for the FastAPI interface.
- [src/discern/api/schemas.py](src/discern/api/schemas.py): Strict Pydantic v2 validation schemas (`Variant` enforcing allowed `GENES` + GARC `valid_mutation` + `0 <= frs <= 1`; `PredictRequest` enforcing safe `sample_id` regex `^[A-Za-z0-9_-]{1,40}$`, `0.001 <= prevalence <= 0.5`, and valid `assay` keys).
- [src/discern/api/app.py](src/discern/api/app.py): FastAPI web server (`634` lines) providing all `21` HTML, JSON, CSV, health, replay, private upload queue, batch surveillance, longitudinal timeline, panel audit, evidence, and model-card endpoints with private cookie-isolated `UploadedBatch` storage.

### 10.3 Web User Interface: `src/discern/ui/` (`14` HTML Templates + `3` Static Assets)

#### Templates (`src/discern/ui/templates/`)
- [base.html](src/discern/ui/templates/base.html): Global page shell, research-use banner header, local CSS/JS links with cache-busting query versions, and responsive sidebar/content container.
- [navigation.html](src/discern/ui/templates/navigation.html): Left navigation rail links (`/`, `/queue`, `/performance`, `/priority`, `/panels`, `/model-card`, `/about`) with inline SVG icons.
- [home.html](src/discern/ui/templates/home.html): Primary `/` sample-intake screen featuring CSV/native file upload, prevalence slider, optional sequencing-assay selector, held-out replay dropdown, and one-click illustrative timeline launcher.
- [assay_options.html](src/discern/ui/templates/assay_options.html): Shared `<option>` dropdown partial listing all `5` sequencing assays (`wgs`, `deeplex_myc_tb`, `deeplex_xl`, `ampore_tb`, `icmr_nirt_minion`).
- [result.html](src/discern/ui/templates/result.html): Detailed case evaluation screen (`596` lines) displaying final R/S/U verdict, calibrated score gauge & prevalence recalculator, Bayes PPV/NPV grid, mechanism reflex guidance, warning alerts, variant evidence cards, and SVG gene-domain lollipop tracks.
- [queue.html](src/discern/ui/templates/queue.html): Batch prioritization table (`/queue`, `254` lines) showing MLLS unsupervised batch prevalence estimate + inline SVG PPV curve, live quality audit card, recurrence/contamination surveillance panel, status filter, and CSV export.
- [patient_timeline.html](src/discern/ui/templates/patient_timeline.html): Interactive serial-isolate timeline (`/patient/{id}`) featuring an inline SVG Muller-style clonal-ribbon chart, scrubber slider, playback toggle, and selectable per-visit variant breakdown.
- [panels.html](src/discern/ui/templates/panels.html): Sequencing-Assay & tNGS Panel Coverage benchmark page (`/panels`) comparing blind spots and detection limits from `e24_panel_coverage`.
- [priority.html](src/discern/ui/templates/priority.html): Variant Referral Priority List page (`/priority`) showing the top unresolved variant combinations (`Rv0678@c-11a`, `pepQ@G197R`, etc.) and cumulative U-sample coverage from `E12`.
- [performance.html](src/discern/ui/templates/performance.html): Study validation dashboard (`/performance`, `277` lines) aggregating catalogue, model, cross-source, LOSO, novel-variant, prevalence-shift, MIC-noise, ship-decision, and external KM metrics.
- [model_card.html](src/discern/ui/templates/model_card.html): Renders [docs/MODEL_CARD.md](docs/MODEL_CARD.md) sections inside the local web chrome (`/model-card`).
- [about.html](src/discern/ui/templates/about.html): Research-use governance statement, offline operation summary, and methodological glossary page (`/about`).
- [table.html](src/discern/ui/templates/table.html): Reusable semantic HTML table partial with scrollable wrapper and human-readable column headers.
- [error.html](src/discern/ui/templates/error.html): Clean HTTP `4xx`/`5xx` error template displaying located upload error categories without leaking rejected values.

#### Static Assets (`src/discern/ui/static/`) — 100% Self-Contained Offline
- [app.css](src/discern/ui/static/app.css): Complete stylesheet (`1,597` lines) using system fonts and semantic light tokens (`Slate`/`Sky`/`Emerald`/`Red`/`Amber` WCAG AA contrast >= 7:1).
- [app.js](src/discern/ui/static/app.js): Pure vanilla JavaScript (`1,064` lines) driving real-time prevalence-slider log-odds & Bayes PPV/NPV updates, responsive nav toggling, queue filtering, file-picker feedback, and timeline scrubber playback (respecting `prefers-reduced-motion`).
- [discern.svg](src/discern/ui/static/discern.svg): Local vector brand icon (`#245744` rounded square emblem).

### 10.4 Evaluation, Benchmark & Build Scripts: `eval/` (`35` scripts)

- [eval/make_manifest.py](eval/make_manifest.py): Scans `data/raw/tb-bdq-cat-main/{data,catalogues,wildcards}` and records SHA-256 hashes + byte sizes into `MANIFEST.json`.
- [eval/e01_prep.py](eval/e01_prep.py): Executes cohort preparation and checks exact E01 row/patient/mutation counts.
- [eval/e02_who.py](eval/e02_who.py): Benchmarks the published WHO v2 GARC catalogue as-is at `FRS >= 0.75` and `FRS >= 0.1`.
- [eval/e03_catalogue_cv.py](eval/e03_catalogue_cv.py): Builds and cross-validates Tier A `catomatic` catalogues per S1 fold (and computes optimistic in-sample baseline).
- [eval/e04_lr_f2.py](eval/e04_lr_f2.py): Runs main S1 nested patient-grouped CV for Tier B L2 LR (`F2`), 500-rep patient bootstrapping, sensitivity-at-fixed-specificity, and saves `oof_predictions.parquet` + `demo/held_out_samples.csv`.
- [eval/e05_ablation.py](eval/e05_ablation.py): Compares feature sets `F0`, `F1`, `F2`, and `F2 @ FRS >= 0.75` on identical S1 folds.
- [eval/e05b_model_comparison.py](eval/e05b_model_comparison.py): Head-to-head comparison of L2 LR, L1 LR, Random Forest, HistGradientBoosting, and Linear SVM on fixed S1 F2 features.
- [eval/e06_cross_source.py](eval/e06_cross_source.py): Runs bidirectional S3 cross-source generalization (`CRyPTIC -> NICD` and `NICD -> CRyPTIC`) with strictly nested inner calibration.
- [eval/e07_loso.py](eval/e07_loso.py): Evaluates Leave-One-Site-Out (`S2`) generalization and records per-site confusion/specificity metrics (`per_site.csv`).
- [eval/e08_novel_variant.py](eval/e08_novel_variant.py): Runs variant-set grouped (`S4`) novel-missense evaluation comparing region+count, PDB 3D distance/burial, and exploratory biophysical properties.
- [eval/e09_hybrid.py](eval/e09_hybrid.py): Exploratory analysis applying nested fold `lo`/`hi` thresholds inside the `661` Tier A `U` samples.
- [eval/e10_prevalence_shift.py](eval/e10_prevalence_shift.py): Audits raw vs. prior-odds-adjusted ECE (10 quantile bins) and Brier score across S3 domain shifts.
- [eval/e11_imbalance.py](eval/e11_imbalance.py): Evaluates `1:1` balanced training degradation and downsampled resistant learning curves (`0.10` to `1.00`).
- [eval/e12_priority.py](eval/e12_priority.py): Ranks unresolved `U` variant combinations by resistance frequency and carrier count, emitting priority tables for `/priority`.
- [eval/e13_ppv_npv.py](eval/e13_ppv_npv.py): Generates Bayes PPV/NPV tables across `1%`, `2%`, `5%`, `10%` prevalences using CV sensitivity and specificity.
- [eval/e14_conformal.py](eval/e14_conformal.py): Audits Bentkus conformal abstention guarantees (`conformal.py`) and computes empirical-Bayes variant shrinkage (`shrinkage.py`).
- [eval/e15_voi.py](eval/e15_voi.py): Simulates retrospective label-reveal sequencing queues comparing Beta-entropy VOI vs simple count ordering.
- [eval/e16_mechanism.py](eval/e16_mechanism.py): Replicates Snobre et al. (2024) `Rv0678 + mmpL5/S5` co-occurrence frequency (`25.04%`) and evaluates the S1 ablation of the epistatic veto.
- [eval/e17_mic_audit.py](eval/e17_mic_audit.py): Audits quantitative MIC distributions (`borderline_by_mic.csv`), exposing `454` isolates (`191 R`) near critical concentration cutoffs.
- [eval/e18_expanded_training.py](eval/e18_expanded_training.py): Tests whether retraining on the expanded `23,441`-isolate CRyPTIC v3.4 table improves S1/KM metrics (negative result).
- [eval/e18b_high_only.py](eval/e18b_high_only.py): Isolates sample quantity from MEDIUM-quality label noise by training on `16,923` HIGH-only CRyPTIC v3.4 samples (still degrades AUC to `0.6668`).
- [eval/e18_head_to_head.py](eval/e18_head_to_head.py): Same-data benchmarks comparing patient-grouped vs shuffled CV (`g7_adlard`) and WHO v2 across `FRS 0.1/0.25/0.75` (`g7_who_v2`).
- [eval/km_cohort.py](eval/km_cohort.py): Shared helper (`load_km_cohort`) loading the held-out `792`-isolate KM validation dataset (`35 R`) into Discern format.
- [eval/ext_km_validate.py](eval/ext_km_validate.py): Evaluates both Tier A catalogue and Tier B L2 model cleanly on the unseen KM external cohort.
- [eval/ext_validate.py](eval/ext_validate.py): Documents accession inventories and honest access blockers for independent Sobol (`744`) and Choi (`49`) datasets.
- [eval/e19_ship_decision.py](eval/e19_ship_decision.py): Evaluates pre-registered ship-gate criteria (`D011`) for shipping the audited epistatic veto using S1 and KM outcomes.
- [eval/e21_tuned_comparison.py](eval/e21_tuned_comparison.py): Nested hyperparameter grid tuning across RF, HGB, L1 LR, L2 LR, and SVM on S1 + KM under D012/D013 ship criteria.
- [eval/e22_biophysics_full_cohort.py](eval/e22_biophysics_full_cohort.py): Tests whether `structure.py` biophysical descriptors improve S1/KM performance across all `11,868` isolates (negative result: `0.90578`).
- [eval/e23_novel_pipeline.py](eval/e23_novel_pipeline.py): Standalone research prototype exploring lineage-filtered rare `mmpL5` mutations, Lineage 2/4 indicators, and continuous `log2(MIC)` regression (`KM AUC 0.89355`).
- [eval/panel_audit.py](eval/panel_audit.py): Computes theoretical sensitivity ceilings and blind spots across all `5` sequencing panel profiles (`e24_panel_coverage`).
- [eval/g8_context.py](eval/g8_context.py): Generates consolidated Wilson confidence intervals across archived evaluation tables (`g8_evaluation_context`).
- [eval/run_all.py](eval/run_all.py): Ordered CLI runner executing `e01_prep.py` through `e13_ppv_npv.py` (or a chosen subset).
- [eval/build_production.py](eval/build_production.py): Trains production Tier A catalogue + Tier B L2 LR + Platt calibrator + thresholds + evidence store on all `11,868` samples into `models/lr_model/`.
- [eval/generate_docs.py](eval/generate_docs.py): Regenerates [docs/MODEL_CARD.md](docs/MODEL_CARD.md) and [docs/RESULTS.md](docs/RESULTS.md) straight from `results/` JSON/CSV files (never hand-edited).
- [eval/verify.py](eval/verify.py): Independent acceptance checker verifying mathematical consistency and tolerance adherence across all `14` acceptance records (`make verify`).

### 10.5 Test Suite: `tests/` (`15` files, `98` tests total)

- [tests/test_archive_checks.py](tests/test_archive_checks.py): Core invariant suite (GARC grammar flags, E01/E02/E03 metrics, strict train/test fold ID disjointness, calibration error bounds, banner presence, and anti-leakage vocabulary check).
- [tests/test_backend_api_wiring.py](tests/test_backend_api_wiring.py): Verifies live JSON `/predict`, `/batch` uploads & replay, SVG domain track context on replay, in-memory report memoization, `/evidence` validation, and content-negotiated errors.
- [tests/test_backend_hardening.py](tests/test_backend_hardening.py): Ensures upload limit enforcement, non-echoing CSV syntax errors, replay reference prevalence exposure, archive-wide ternary PPV/NPV rates, strict legacy JSON validation, and safe `U` deferral on catalogue failures.
- [tests/test_batch_recurrence.py](tests/test_batch_recurrence.py): Unit tests batch surveillance deduplication by `patient_id`, rare-variant WGS escalation, and homopolymer lab contamination alerts.
- [tests/test_e21_tuned_comparison.py](tests/test_e21_tuned_comparison.py): Verifies `D013` promotion guard requiring `>0.005` S1 improvement while preserving external KM AUC.
- [tests/test_explorer.py](tests/test_explorer.py): Tests empty variant profile behavior (`NO_VARIANT_INVISIBLE_RESISTANCE_POSSIBLE`), immutable replay behavior, and rejection of modified replay profiles.
- [tests/test_importers_and_reflex.py](tests/test_importers_and_reflex.py): Validates HGVS→GARC translations, TB-Profiler JSON, Deeplex Myc-TB CSV, AmPORE-TB TSV/CSV parsing, clinical reflex guidance cards, and web review upload.
- [tests/test_indel_canon.py](tests/test_indel_canon.py): Verifies clean H37Rv `Rv0678` ORF translation, training indel canonical consistency, and mapping diverse HGVS indel spellings to matching training mutations.
- [tests/test_label_shift.py](tests/test_label_shift.py): Tests Saerens EM MLLS prevalence estimation accuracy, small-batch fallback gating (`n < 20`), and chronological patient deduplication.
- [tests/test_longitudinal.py](tests/test_longitudinal.py): Tests relative/ISO date parsing, de-novo/rising resistance trajectories, polyclonal competition, 4-col vs 6-col backward compatibility, and `/patient/{id}` rendering.
- [tests/test_novelty_biophysics_and_epistasis.py](tests/test_novelty_biophysics_and_epistasis.py): Unit and regression tests for Grantham/charge calculations, α4 recognition helix flags, subclonal heteroresistance alerts, static high-confidence epistatic veto triggering, and `atpE`/`pepQ` protection.
- [tests/test_panels.py](tests/test_panels.py): Tests all `5` sequencing panel definitions, missed-locus blind spots, homopolymer warnings, below-LOD alerts, and `Predictor` assay integration.
- [tests/test_upload_integration.py](tests/test_upload_integration.py): End-to-end tests covering private cookie-scoped batch links, assay retention across replay, manual prevalence overrides, Deeplex fractional-percent frequency scaling (`<1%`), AmPORE filename IDs, malformed uploads, and longitudinal timelines.
- [tests/test_upload_quality_and_live_scoring.py](tests/test_upload_quality_and_live_scoring.py): Tests live CSV quality completeness scoring, live Tier A+B scoring for novel uploads, and strict separation of quality scores from resistance predictions.
- [tests/test_verify.py](tests/test_verify.py): Meta-tests for `eval/verify.py` covering exact tolerance boundaries, numeric diff calculations, `pass_all: null` handling, and empty/inconsistent index rejection.

### 10.6 Frozen Models (`models/`), Docs (`docs/`), Demo Data (`demo/`), Local Data (`data/`), & Results (`results/`)

#### Production & Experimental Artifacts (`models/` — `7` files)
- `models/lr_model/catalogue_bdq.csv`: Production GARC rule catalogue (`331` lines) built on all `11,868` isolates.
- `models/lr_model/catalogue_version.json`: Catalogue build settings (`NC_000962.3`, `frs=0.1`, `Binomial two`, `manifest_sha256`).
- [models/lr_model/bdq_lr.joblib](models/lr_model/bdq_lr.joblib): Pickled scikit-learn production bundle (`C=0.3` L2 LogisticRegression + Platt calibrator).
- [models/lr_model/metadata.json](models/lr_model/metadata.json): Frozen feature schema (`120` ordered feature names), `pi_tr: 0.072969`, thresholds (`t_bin=0.023, lo=0.0549, hi=0.7273`), and verified E03/E04 hashes.
- `models/lr_model/evidence_store.csv`: Empirical resistance statistics (`698` variant rows + header) with total carriers, `n_resistant`, `pct_resistant`, Wilson 95% CIs, and split `CRyPTIC` vs `NICD` counts.
- `models/lr_mechanism/lof_effects.csv`: Audited pre-computed `mmpL5`/`mmpS5` truncating LoF table (`9.5 KB`) supplying `veto_high_conf` to live `Predictor`.
- `models/lr_shrinkage/variant_shrinkage.csv`: Offline empirical-Bayes Beta-Binomial shrinkage estimates across `1,042` variants (`127 KB`, research artifact from `E14`).

#### Documentation (`docs/` — `6` files + `assets/`)
- [docs/DECISIONS.md](docs/DECISIONS.md): Pre-registered engineering/methodology log (`D001–D013`, locking environment, cohorts, enhanced-mode isolation, and ship criteria).
- [docs/MODEL_CARD.md](docs/MODEL_CARD.md): Auto-generated clinical & technical model card rendered at `/model-card`.
- [docs/RESULTS.md](docs/RESULTS.md): Auto-generated tables (`979` lines) quoting every metric and CSV link across `results/`.
- [docs/PROGRESS.md](docs/PROGRESS.md): Milestone verification logs (`G0` onward) recording exact evaluation outputs and fidelity notes.
- [docs/DEMO_SCRIPT.md](docs/DEMO_SCRIPT.md): Structured live walkthrough script covering upload, unresolved evidence cards, prevalence sensitivity, queue filtering, and verification.
- [docs/GLOSSARY.md](docs/GLOSSARY.md): Plain-language reference table of core terms.
- `docs/assets/discern-logo.svg`: Project brand emblem displayed in [README.md](README.md).

#### Bundled Demos (`demo/` — `3` CSVs)
- [demo/example_batch.csv](demo/example_batch.csv): `19`-row multi-sample example batch showcasing `192_ins_g`, subclonal `G41V`/`R134!`/`R50G`, promoter/indel mixtures, and `mmpL5` background variants.
- [demo/held_out_samples.csv](demo/held_out_samples.csv): `43`-row curated subset (`SAMPLE_ID`, `cat_call`, `p`, `y`) mapped to `oof_predictions.parquet` for zero-leakage UI replay.
- [demo/longitudinal_example.csv](demo/longitudinal_example.csv): `7`-row 6-column dated patient trajectory dataset (`PT_A_V0/V1/V2` showing subclonal emergence `M146T` → `141_ins_c` rising from `0.12` to `0.88`, alongside stable `PATIENT_B`).

#### Data (`data/`) & Results (`results/` — `211` files)
- `data/raw/` & `data/external_cryptic_v34/`: Downloaded tarballs, raw gzip tables, `MANIFEST.json`, `Rv0678_WT.pdb`, WHO v2 CSV, wildcards JSON, and CRyPTIC v3.4 parquets (git-ignored).
- `data/processed/`: Generated `samples.parquet`, `mutations.parquet`, `s1_folds.npy`, `e03_oof_calls.parquet`, `e04_oof.parquet`, and expanded experiment parquets (git-ignored).
- `data/choi_runs.tsv` & `data/sobol_runs.tsv`: External study accession manifests (`49` Korea and `744` SA runs).
- [results/INDEX.md](results/INDEX.md) & `results/index.json`: Master index across **`211` files** partitioned into `results/models/` (`hgb_model`, `hgb_tuned`, `lr_expanded`, `lr_l1`, `lr_l1_tuned`, `lr_l2_tuned`, `lr_mechanism`, `lr_model`, `lr_shrinkage`, `rf_model`, `rf_tuned`, `structure_test`, `svm_model`, `svm_tuned`) and `results/approaches/` (`catalogue`, `conformal`, `data_prep`, `external_choi`, `external_km`, `external_sobol`, `head_to_head`, `hybrid`, `mic_audit`, `panel_audit`, `predictive_values`, `prevalence_correction`, `priority_list`, `research_archive`, `ship_decision`, `voi_queue`, `who_baseline`), plus `results/oof_predictions.parquet`.

---

## 11. Developer Workflows, Routes, CLI & FAQ

### 11.1 Launching, Verifying & Developing Locally

Run from the repository root (`d:\Desktop\project\Discern`) using the project environment (`C:\Users\shaah\miniconda3\envs\discern\python.exe` or your activated venv):

```powershell
# 1. Start the local offline web workbench (with zero request-body logging)
.\serve.ps1
# or explicitly:
python -m uvicorn discern.api.app:app --app-dir src --host 127.0.0.1 --port 8000 --no-access-log

# 2. Verify all 14 frozen acceptance tests + 98 pytest tests (same as `make verify`)
python eval/verify.py
python -m pytest tests/ -q

# 3. Regenerate docs/MODEL_CARD.md and docs/RESULTS.md from results/
python eval/generate_docs.py

# 4. Build & run completely offline inside Docker
docker build -t discern:local .
docker run --rm -p 8000:8000 discern:local
```

### 11.2 CLI Subcommands ([src/discern/cli.py](src/discern/cli.py))

```bash
discern replay --sample-id <id>                  # Print immutable archived OOF prediction
discern batch <input.csv> <output.csv>           # Score a 4-column CSV and rank U referrals
discern eval --all                               # Trigger eval/run_all.py across E01-E13
discern verify                                   # Run the 14/14 acceptance verifier
```

### 11.3 Complete HTTP API Route Reference ([src/discern/api/app.py](src/discern/api/app.py), `21` Routes)

| Method & Path | Response | Purpose |
|---|---|---|
| `GET /health` | JSON | Offline health check (`status`, `banner`). |
| `POST /predict` | JSON | Strictly schema-validated live single-isolate triage (`sample_id`, `prevalence`, `assay`, `variants[]`). |
| `POST /batch` | CSV | Programmatic batch CSV scoring returning calls and descending `probability_at_prevalence` queue rank. |
| `GET /evidence/{gene}/{mutation}` | JSON | Carrier counts, Wilson CI, and separate CRyPTIC vs NICD totals for a validated variant. |
| `GET /replay/{sample_id}` | JSON | Immutable archived out-of-fold score + explanation from `Archive`. |
| `GET /replay-samples.csv`, `/example.csv`, `/longitudinal-example.csv` | CSV | Download bundled demo inputs directly from the workbench. |
| `GET /` | HTML | Primary intake workbench (`home.html`). |
| `POST /review`, `POST /longitudinal-example` | HTML/Redirect | Parse uploaded/native lab exports or bundled longitudinal example into private session store. |
| `GET /result` | HTML | Full single-sample verdict, domain lollipop SVGs, warnings, reflex guidance, and evidence cards (`result.html`). |
| `GET /queue`, `POST /queue` | HTML | Batch review queue (`/queue`) with MLLS prevalence estimation, QC status, recurrence alerts, and CSV upload (`queue.html`). |
| `GET /patient/{patient_id}`, `GET /api/patient/{patient_id}` | HTML / JSON | Interactive serial-isolate clonal-ribbon timeline for dated multi-visit patients (`patient_timeline.html`). |
| `GET /priority` | HTML | Unresolved referral variant combinations priority table from `E12` (`priority.html`). |
| `GET /panels` | HTML | Sequencing assay & tNGS blind-spot coverage matrix from `E24` (`panels.html`). |
| `GET /performance` | HTML | Pre-registered validation metrics and uncertainty intervals from `results/` (`performance.html`). |
| `GET /model-card` | HTML | Generated model card rendered from `docs/MODEL_CARD.md` (`model_card.html`). |
| `GET /about` | HTML | Research-use governance and plain-English glossary (`about.html`). |

### 11.4 Quick Troubleshooting & Beginner FAQ

- **"Why does `results/INDEX.md` say `14/49 tests pass`, while `eval/verify.py` says `14/14 acceptance tests passed`?"**
  `INDEX.md` totals *all* `49` registered experiments including descriptive studies (`pass_all: null`), negative baseline ablations (`e18` degraded expanded training), and ENA-blocked cohorts. `verify.py` checks only the `14` acceptance-gated experiments (`E01–E13`), and **all `14/14` pass**.
- **"Why did `Rv0678@Q131!` get a `30.2%` model probability yet still say `Send for lab test` (`U`) instead of `S` or `R`?"**
  Rule 1: **Tier B never overrides Tier A.** `Q131!` is absent from the curated catalogue (`NOVEL_VARIANT`), so Discern sends it to culture and gives you `30.2%` (`60.7`th queue percentile) so your lab tests high-risk referrals before lower-risk ones.
- **"Why wasn't `model` swapped to Random Forest (`0.9103` AUC) or Gradient Boosting (`0.8892` KM AUC)?"**
  Pre-registered ship rule `D013` required beating L2 LR (`0.9085`) by more than `+0.005` on S1 while keeping KM AUC `>= 0.8682`. No complex model cleared both bars (`E21`), so Discern keeps interpretable, fast, linearly additive L2 Logistic Regression.
- **"What happens if an isolate has `Rv0678@141_ins_c` AND an inactivated `mmpL5` pump?"**
  In live mode, if the `mmpL5`/`mmpS5` mutation is a pre-audited high-confidence truncating variant in `models/lr_mechanism/lof_effects.csv` with `FRS >= Rv0678_FRS - 0.20` (and no `atpE`/`pepQ` `R` is present), the Epistatic Veto overrides the `Rv0678` `R` to **`S`** and logs **`EPISTATIC_PUMP_INACTIVATION`**. If the selected assay is `deeplex_myc_tb`, `ampore_tb`, or `icmr_nirt_minion`, Discern cannot observe `mmpL5`, so it emits **`EPISTATIC_VETO_BLIND`**.
