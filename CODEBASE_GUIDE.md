# Discern: The Complete Beginner-to-Expert Codebase Guide

> **Research use only. Clinical decisions require independent laboratory confirmation.**

---

## Table of Contents
1. [The Big Picture: What Problem Are We Solving?](#1-the-big-picture-what-problem-are-we-solving)
2. [Biological Context](#2-biological-context)
3. [The Core Philosophy: Why Discern is Built Differently](#3-the-core-philosophy-why-discern-is-built-differently)
4. [High-Level Architecture: The Two-Tiered Brain](#4-high-level-architecture-the-two-tiered-brain)
5. [Data Ingestion & Mutation Grammar](#5-data-ingestion--mutation-grammar)
6. [The Core Engine: Module-by-Module Walkthrough](#6-the-core-engine-module-by-module-walkthrough)
7. [Five Research Workflows](#7-five-research-workflows)
8. [The Web Application & API (No CDNs, Strict CSP)](#8-the-web-application--api-no-cdns-strict-csp)
9. [The Evaluation Harness & "Golden Mode" Invariance](#9-the-evaluation-harness--golden-mode-invariance)
10. [Walkthrough: A Sample's Complete Journey Through Discern](#10-walkthrough-a-samples-complete-journey-through-discern)
11. [Complete File-by-File Directory Reference](#11-complete-file-by-file-directory-reference)

---

## 1. The Big Picture: What Problem Are We Solving?

Discern compares reported tuberculosis variants with archived estimates of bedaquiline resistance. A saved study case displays its original held-out call, laboratory phenotype, and model score. A live upload displays a research result for the submitted profile.

- **Resistant (R)** flags a sample for confirmation.
- **Susceptible (S)** is the catalogue's classification for the submitted profile.
- **Send for lab test (U)** identifies an unresolved result. The calibrated model score orders unresolved samples for research review.

## 2. Biological Context

The input is a report of observed variants across bedaquiline-relevant loci. The [model card](docs/MODEL_CARD.md) describes the research cohorts, model limitations, and reliance on laboratory confirmation. Assay coverage and reported allele fractions affect how much evidence is available for interpreting an upload.

---

## 3. The Core Philosophy: Why Discern is Built Differently

In machine learning, people often prioritize training massive deep-learning models. In clinical TB genomics, that approach fails catastrophically for three reasons:

### Hard Principle 1: Never Invent or Round Numbers
Every statistic in Discern's web pages, user interface, and documentation is generated programmatically from archived cross-validation result files (`results/*.json`). No developer or researcher is permitted to hand-edit a metric.

### Hard Principle 2: Patient-Grouped Cross-Validation (No Leakage)
In TB research, a single patient might provide multiple sputum samples over several months. If Sample 1 from Patient $A$ is in the training set and Sample 2 from Patient $A$ is in the test set, the model memorizes Patient $A$'s unique genetic quirks.
Discern enforces strict **patient grouping (`SUBJID`)**: all samples from a single patient must stay together in either training or testing—never split across both.

### Hard Principle 3: "Admit When You Don't Know"
A bad model makes a confident guess on everything. In medicine, a wrong guess can lead to patient death. Discern uses a **ternary decision policy**: $R$, $S$, or $U$ (Unresolved). If an isolate has an unseen mutation or falls into an ambiguous zone, Discern explicitly refuses to guess, flagging it for laboratory confirmation.

### Hard Principle 4: Offline by Design (Zero Internet, Zero CDNs)
Discern uses local artefacts and assets at runtime.
- No Google Fonts.
- No external Bootstrap or Tailwind CDNs.
- No telemetry or cloud tracking.
- The app does not log request bodies. Live uploads are held in a bounded in-memory cache and scoped to the uploading browser.

---

## 4. High-Level Architecture: The Two-Tiered Brain

Discern uses a hybrid, multi-layered architecture:

```
┌────────────────────────────────────────────────────────┐
│                      INPUT DATA                        │
│   List of mutations + Optional Assay/Patient Metadata │
└───────────────────────────┬────────────────────────────┘
                            │
                            ▼
┌────────────────────────────────────────────────────────┐
│           TIER A: Deterministic Catalogue              │
│   Statistical Binomial rules built with catomatic      │
│   Checks known WHO & curated resistance mutations      │
└─────────────┬────────────────────────────┬─────────────┘
              │ Calls R or S               │ Calls U (Uncertain)
              ▼                            ▼
┌───────────────────────────┐ ┌──────────────────────────┐
│   EPISTATIC VETO FILTER   │ │   TIER B: Calibrated ML  │
│ Checks if mmpL5 is broken │ │ L2 Logistic Regression   │
│ If pump broken -> Downgrade│ │ Outputs raw logit score │
└─────────────┬─────────────┘ └────────────┬─────────────┘
              │                            │
              │                            ▼
              │               ┌──────────────────────────┐
              │               │      Platt Scaling       │
              │               │ Inner out-of-fold fit    │
              │               │ Calibrated probability   │
              │               └────────────┬─────────────┘
              │                            │
              ▼                            ▼
┌────────────────────────────────────────────────────────┐
│            BAYESIAN PREVALENCE ADJUSTMENT              │
│    Transforms probability to match local clinic rate   │
│    Calculates Negative & Positive Predictive Values    │
└───────────────────────────┬────────────────────────────┘
                            │
                            ▼
┌────────────────────────────────────────────────────────┐
│                    FINAL CALL & UI                     │
│  Triage Call (R / S / U) + Evidence + Warnings + Queue │
└────────────────────────────────────────────────────────┘
```

### Tier A (The Rule-Based Anchor)
- Evaluates non-synonymous mutations using a statistical binomial resistance catalogue built on 11,868 samples.
- If any variant has proven statistical association with resistance, Tier A calls **Resistant ($R$)**.
- If all variants are known benign variants, Tier A calls **Susceptible ($S$)**.
- If a mutation is rare, unseen, or contested, Tier A calls **Unresolved ($U$)**.

### Tier B (The Machine Learning Calibrator)
- An $L_2$-penalized Logistic Regression model ($C=0.3$) trained on gene aggregates (missense counts, frameshift flags, DNA-binding domain alterations) and recurrent variant indicators.
- Tier B **never overrides Tier A to call resistance on its own**. Instead, Tier B:
  1. Computes the continuous probability of resistance.
  2. Ranks unresolved ($U$) samples so laboratory staff know which samples are most urgent.
  3. Raises disagreement warnings if Tier A and Tier B conflict.

---

## 5. Data Ingestion & Mutation Grammar

How does Discern read mutation strings? Let's look at [`src/discern/grammar.py`](src/discern/grammar.py).

### The Mutation Grammar
Mutations are formatted using the standard GARC format: `GENE@MUTATION` (for example, `Rv0678@141_ins_c`).

| Mutation String | Mutation Type | What It Means in Plain English |
|---|---|---|
| `Rv0678@R50Q` | Missense | Amino acid 50 was Arginine (R), but changed to Glutamine (Q). |
| `Rv0678@Q131!` | Stop Codon (Nonsense) | Amino acid 131 turned into a STOP signal (`!`), truncating the protein. |
| `Rv0678@141_ins_c` | Insertion (Indel) | A cytosine nucleotide was inserted at position 141, causing a frameshift. |
| `Rv0678@193_del_g` | Deletion (Indel) | A guanine nucleotide was deleted at position 193. |
| `Rv0678@c-11a` | Promoter | In the promoter region 11 bases before the gene start, Cytosine became Adenine. |
| `Rv0678@S68S` | Synonymous | Amino acid 68 was Serine and stayed Serine. Harmless silent mutation! |
| `Rv0678@S68Z` | Null / Uncallable | Sequencer could not read position 68. Triggers a data warning! |

### The Fraction of Reads Supporting (FRS)
When sequencing sputum, not every bacterium is identical. If 80% of reads show the mutation and 20% show wild-type, the **Fraction of Reads Supporting (FRS)** is `0.80`.
- Mutations with $FRS \ge 0.75$ are considered **fixed** (dominant clone).
- Mutations with $0.10 \le FRS < 0.75$ are **minor alleles** (subclonal heteroresistance).
- Mutations with $FRS < 0.10$ are treated as sequencing noise and filtered out.

---

## 6. The Core Engine: Module-by-Module Walkthrough

Let's explore the core Python files in [`src/discern/`](src/discern/).

### 1. `prep.py` — Raw Data Harmonization
- Reads raw phenotypic DST measurements from two massive public cohorts:
  - **CRyPTIC** (10,704 samples; 83 resistant): Global multi-country surveillance using microtiter plates. Only `QUALITY == "HIGH"` is kept.
  - **NICD** (1,164 samples; 783 resistant): South African National Institute for Communicable Diseases cohort using MGIT960 liquid culture, heavily enriched for drug resistance.
- Merges duplicates prioritizing resistant calls and validated minimum inhibitory concentrations (MIC).
- Generates clean, fast-loading Parquet files: `data/processed/samples.parquet` and `data/processed/mutations.parquet`.

### 2. `splits.py` — Patient-Safe Cross-Validation
- Defines 4 distinct evaluation split strategies:
  - **S1 (Main 5-Fold Patient CV):** Uses `GroupKFold(n_splits=5)` grouping on `SUBJID`.
  - **S2 (Leave-One-Site-Out / LOSO):** Holds out each of the 23 clinical sites one by one to test geographic generalization.
  - **S3 (Cross-Source Shift):** Train on CRyPTIC, test on NICD (and vice versa).
  - **S4 (Variant-Set Grouping):** Groups by unique combinations of mutations to test performance on unseen variants.

### 3. `catalogue.py` — Tier A Rule Builder
- Interfaces with the `catomatic` and `piezo` bioinformatics libraries.
- Constructs two-tailed Binomial test catalogues with strict unclocking rules and wildcards.
- Outputs `models/lr_model/catalogue_bdq.csv`.

### 4. `features.py` — Numerical Representation
Transforms genomic text into a numerical matrix for machine learning:
1. **Vocabulary features:** Only mutations that appear in $\ge 3$ training samples get their own column.
2. **Structural / Functional Aggregates:**
   - `rv_missense`, `rv_indel`, `rv_stop`, `rv_prom`: Maximum FRS of each mutation class in *Rv0678*.
   - `rv_mis_dbd`: Missense mutation inside the *Rv0678* DNA-binding domain (amino acids 34 to 99).
   - `atpE_any`, `pepQ_any`: Mutations in the target or peptidase genes.
   - `mmpL5_LoF`: Presence of an inactivating indel or stop codon in the pump gene *mmpL5*.
   - `rv_x_L5LoF`: The interaction term ($Rv0678 \times mmpL5 \text{ LoF}$).

### 5. `model.py` — Tier B Classifier & Platt Scaling
- Trains `LogisticRegression(C=0.3, penalty='l2', max_iter=4000)`.
- **Nested Out-of-Fold Calibration:** Raw logistic regression scores are not true probabilities. To fix this, Discern uses Platt Scaling: fitting a secondary logistic regression on the log-odds of inner out-of-fold predictions.
- **Threshold Selection:**
  - $t_{\text{bin}}$: Youden index maximizing $(\text{Sensitivity} + \text{Specificity} - 1)$.
  - $lo$: Highest threshold ensuring Negative Predictive Value ($\text{NPV}) \ge 0.985$.
  - $hi$: Lowest threshold ensuring Positive Predictive Value ($\text{PPV}) \ge 0.800$.

### 6. `prevalence.py` — Bayesian Probability Correction
A probability calculated on a dataset where 7.3% of strains are resistant is invalid in a clinic where only 2% or 20% are resistant!
Discern applies Bayes' rule in odds format:
$$\text{odds}(p) = \frac{p}{1 - p}$$
$$\text{odds}' = \frac{\text{odds}(p)}{\text{odds}(\pi_{\text{train}})} \times \text{odds}(\pi_{\text{clinic}})$$
$$p_{\text{adjusted}} = \frac{\text{odds}'}{1 + \text{odds}'}$$
These predictive values depend on the selected prevalence and the archived call rates.

### 7. `mechanism.py` and `structure.py`

`mechanism.py` implements the configured interpretation rules used in the research audit. Allele-fraction similarity supports a flag and cannot establish physical linkage. `structure.py` reads local structural annotations for exploratory evaluations. See the [model card](docs/MODEL_CARD.md) for the evidence and limitations.

---

## 7. Five Research Workflows

### Addition 1: Sequencing-Assay & tNGS Panel Blind-Spot Auditor
* **File:** [`src/discern/panels.py`](src/discern/panels.py)
* **The Clinical Problem:** Different sequencing kits sequence different genes! Standard **Deeplex Myc-TB** and **ICMR-NIRT MinION** only sequence *Rv0678*—they do not sequence *mmpL5* or *mmpS5*. If a tool doesn't know this, it assumes unsequenced genes are wild-type, silently disabling the Epistatic Veto! Furthermore, Oxford Nanopore sequencers struggle with homopolymers (stretches of identical letters like `GGGGGG` at nt 192–198 of *Rv0678*), frequently reporting false frameshift insertions.
* **The Solution:**
  - Profiles 5 major platforms: WGS, Deeplex Myc-TB, Deeplex XL, AmPORE-TB, and ICMR-NIRT MinION.
  - Warns when pump genes are missing (`EPISTATIC_VETO_BLIND`).
  - Warns when core target genes are unsequenced (`UNSEQUENCED_CORE_LOCI`).
  - Flags low-VAF or Nanopore homopolymer insertions as probable artifacts (`HOMOPOLYMER_ARTIFACT_RISK`).

### Addition 2: Unsupervised Batch Prevalence Estimation (MLLS / Saerens EM)
* **File:** [`src/discern/label_shift.py`](src/discern/label_shift.py)
* **The Clinical Problem:** Doctors using Discern are asked: *"What is the bedaquiline resistance prevalence in your clinic?"* In low-resource settings, doctors have no phenotypic testing data and have to guess!
* **The Solution:**
  - Uses the **Maximum Likelihood Label Shift (Saerens Expectation-Maximization)** algorithm.
  - Estimates prevalence and bootstrap uncertainty from one score per patient, using the earliest dated isolate when dates are available. The method assumes label shift.

### Addition 3: Longitudinal Serial-Isolate Tracking & Clonal Ribbons
* **File:** [`src/discern/longitudinal.py`](src/discern/longitudinal.py)
* **The Solution:**
  - Ingests optional `patient_id` and `collection_date` metadata.
  - Summarizes changes between dated samples. The ribbon displays reported, unphased allele fractions.
  - Renders an interactive **Muller Plot** (clonal ribbon timeline) at `/patient/{patient_id}` showing subclone evolution over time.

### Addition 4: Native Lab-Export Importers & Mechanism-Stratified Reflex Guidance
* **File:** [`src/discern/importers.py`](src/discern/importers.py)
* **The Clinical Problem:** Real hospital sequencers do not output Discern's 4-column CSV. They output TB-Profiler JSON files, Deeplex Myc-TB CSVs, or AmPORE-TB TSV tables using standard HGVS mutation notation (`p.Gln131*` instead of `Q131!`, `c.-11C>A` instead of `c-11a`).
* **The Solution:**
  - Auto-detects and ingests TB-Profiler, Deeplex, and AmPORE files seamlessly.
  - Shows the catalogue rule roles that contributed to the reported call.
  - Maps mutations to a 165-amino-acid structural domain track.
  - Displays research interpretation notes with the result and points readers to laboratory confirmation.

### Addition 5: Lightweight Batch-Recurrence & Contamination Surveillance
* **File:** [`src/discern/surveillance.py`](src/discern/surveillance.py)
* Identifies shared observations across distinct patients in an uploaded batch.
* Surfaces quality-control or WGS follow-up flags. Establishing transmission requires additional evidence.

---

## 8. The Web Application & API (No CDNs, Strict CSP)

Let's look at the web layer in [`src/discern/api/app.py`](src/discern/api/app.py) and [`src/discern/ui/`](src/discern/ui/).

### Technology Stack
- **FastAPI:** High-performance asynchronous Python web framework.
- **Jinja2:** Server-rendered HTML templates.
- **Zero External Dependencies:** All CSS and JavaScript are strictly local. No external web calls are permitted under any circumstances.
- **Strict Content Security Policy (CSP):** `default-src 'self'`. Inline `onclick=` attributes are banned; all event listeners are wired unobtrusively in [`src/discern/ui/static/app.js`](src/discern/ui/static/app.js).

### Public Routes Overview

| Route | HTTP Method | What It Shows / Does |
|---|---|---|
| `/` | `GET` | Home screen: single-isolate variant input, CSV/JSON file upload, and replay selector. |
| `/result` | `GET` | Single isolate triage result: R/S/U call, domain track, evidence breakdown, reflex guidance, and panel coverage strip. |
| `/queue` | `GET`, `POST` | Batch triage queue: priority ranking for laboratory testing, MLLS prevalence curve, and recurrence surveillance alerts. |
| `/panels` | `GET` | Sequencing panel matrix comparing WGS, Deeplex, AmPORE, and MinION coverage and blind spots. |
| `/patient/{patient_id}` | `GET` | Longitudinal timeline: Muller plot (clonal ribbon) and subclone progression table for a specific patient. |
| `/performance` | `GET` | Full scientific audit: every single archived experiment (E01–E24) with complete metrics and confidence intervals. |
| `/model-card` | `GET` | Model card generated from archived report data. |
| `/about` | `GET` | Research context, system architecture, and clinical guidance. |
| `/replay/{sample_id}` | `GET` | Archived held-out result as JSON. |
| `/api/patient/{patient_id}` | `GET` | Timeline JSON for the uploader's batch. |

---

## 9. The Evaluation Harness & "Golden Mode" Invariance

The saved acceptance checks catch drift in the registered evaluation outputs.

### What is the "Golden Mode" (`D003`)?
The core production artifacts:
- Model weights: `models/lr_model/bdq_lr.joblib`
- Production catalogue: `models/lr_model/catalogue_bdq.csv`
- Model metadata: `models/lr_model/metadata.json`

These files are **frozen and cryptographically hashed**. Every experiment in the evaluation harness must reproduce the original numbers within strict scientific tolerance.

### The Acceptance Verifier (`eval/verify.py`)
Run `python eval/verify.py` anytime. It evaluates **14 pre-registered gates**:

| Experiment ID | Experiment Name | Scientific Purpose | Tolerance Check |
|---|---|---|---|
| **E01** | Prep Counts | Verifies exact cohort sample counts (11,868 samples; 866 resistant). | Exact match |
| **E02** | WHO v2 Baseline | Verifies baseline sensitivity of WHO Mutation Catalogue v2. | $\pm 0.005$ |
| **E03** | Tier A CV on S1 | Evaluates 5-fold cross-validated catalogue performance. | $\pm 0.010$ |
| **E04** | Tier B Main (F2) | Evaluates cross-validated ROC-AUC (0.909) and PR-AUC (0.684). | $\pm 0.010$ |
| **E05** | Feature Ablation | Confirms feature additions improve performance over baseline. | $\pm 0.015$ |
| **E06** | Cross-Source Generalization | Evaluates generalization across CRyPTIC and NICD platforms. | $\pm 0.030$ |
| **E07** | Leave-One-Site-Out (LOSO) | Evaluates geographic robustness across 23 international hospitals. | $\pm 0.020$ |
| **E08** | Novel-Variant Test | Evaluates behavior on completely unseen variants (negative result). | Range [0.50, 0.65] |
| **E09** | Hybrid Resolution | Evaluates resolution of Tier A uncertain ($U$) samples by Tier B. | $\pm 15$ samples |
| **E10** | Prevalence Shift | Verifies calibration error drops after Bayesian adjustment. | $\pm 0.010$ |
| **E11** | Imbalance & Learning Curve | Demonstrates class re-weighting does not help and hurts calibration. | $\pm 0.010$ |
| **E12** | Priority Variant Ledger | Audits top variants causing clinical uncertainty. | $\pm 3.0\%$ |
| **E13** | PPV / NPV Tables | Audits positive/negative predictive value across prevalences. | Exact match |

---

## 10. Walkthrough: A Sample's Complete Journey Through Discern

To follow the live path, download [the illustrative longitudinal CSV](demo/longitudinal_example.csv) and submit it from Review queue. The importer validates the file. The upload handler scores the samples and creates a queue scoped to the uploading browser. Open a result from the queue to inspect reported variants and prevalence estimates. Open the visits link to see the dated timeline.

For the archive path, open the saved example from the home page. Its model score and call come from the saved held-out evaluation. See the [demo script](docs/DEMO_SCRIPT.md) for an offline walkthrough.

---

## 11. Complete File-by-File Directory Reference

To help you navigate the codebase like an expert, here is what every single file in the repository does:

### Root Files
* [`AGENTS.md`](AGENTS.md): The master engineering specification and golden numbers.
* [`config.yaml`](config.yaml): Global system configuration (file paths, random seeds, hyper-parameters).
* [`pyproject.toml`](pyproject.toml): Python package build configuration and pytest configuration.
* [`requirements.txt`](requirements.txt): Exact pinned third-party Python dependencies.
* [`Makefile`](Makefile): Quick terminal commands (`make verify`, `make test`, `make app`).
* [`Dockerfile`](Dockerfile): Container build configuration for offline deployment.
* [`CHANGELOG.md`](CHANGELOG.md): Chronicle of all updates and scientific enhancements.

### Core Application (`src/discern/`)
* [`grammar.py`](src/discern/grammar.py): Mutation string parsing and classification.
* [`prep.py`](src/discern/prep.py): Data cleaning, filtering, and merging of CRyPTIC and NICD cohorts.
* [`splits.py`](src/discern/splits.py): Patient-grouped cross-validation splitting algorithms.
* [`catalogue.py`](src/discern/catalogue.py): Tier A deterministic resistance catalogue builder (`catomatic`/`piezo`).
* [`features.py`](src/discern/features.py): Feature matrix extraction (aggregates, vocabulary, interaction terms).
* [`model.py`](src/discern/model.py): Tier B $L_2$ Logistic Regression, nested Platt scaling, and threshold search.
* [`prevalence.py`](src/discern/prevalence.py): Bayesian prevalence adjustment and PPV/NPV calculation.
* [`mechanism.py`](src/discern/mechanism.py): Epistatic veto engine and allele-frequency linkage verification.
* [`structure.py`](src/discern/structure.py): 3D PDB structure parser, Grantham matrix, and $\alpha 4$ helix biophysics.
* [`panels.py`](src/discern/panels.py): Sequencing panel profiles and blind-spot auditor.
* [`label_shift.py`](src/discern/label_shift.py): MLLS (Saerens EM) unsupervised batch prevalence estimation.
* [`longitudinal.py`](src/discern/longitudinal.py): Serial-isolate tracking, clonal dynamics, and Muller plot generator.
* [`importers.py`](src/discern/importers.py): TB-Profiler, Deeplex, and AmPORE file importers and HGVS normalizer.
* [`surveillance.py`](src/discern/surveillance.py): Batch recurrence and artifact alerts.
* [`evidence.py`](src/discern/evidence.py): Historical mutation carrier and resistance database query store.
* [`quality.py`](src/discern/quality.py): Live uploaded CSV data completeness and syntax quality auditor.
* [`decision.py`](src/discern/decision.py): Main runtime prediction pipeline coordinating Tier A, Tier B, and vetoes.
* [`research.py`](src/discern/research.py): Loader for pre-calculated golden cross-validation results.
* [`cli.py`](src/discern/cli.py): Command-line interface (`discern replay`, `discern batch`, `discern verify`).

### Web Interface & API (`src/discern/api/` & `src/discern/ui/`)
* [`src/discern/api/app.py`](src/discern/api/app.py): FastAPI application serving all UI routes and API endpoints.
* [`src/discern/api/schemas.py`](src/discern/api/schemas.py): Pydantic input/output schemas for API requests.
* [`src/discern/ui/templates/base.html`](src/discern/ui/templates/base.html): Base HTML layout with mandatory research-use banner and navigation.
* [`src/discern/ui/templates/home.html`](src/discern/ui/templates/home.html): Input interface (single sample entry, file upload, replay dropdown).
* [`src/discern/ui/templates/result.html`](src/discern/ui/templates/result.html): Detailed triage result page with domain track and reflex guidance.
* [`src/discern/ui/templates/queue.html`](src/discern/ui/templates/queue.html): Prioritized laboratory testing queue with MLLS curve.
* [`src/discern/ui/templates/panels.html`](src/discern/ui/templates/panels.html): 5-gene sequencing panel comparison table.
* [`src/discern/ui/templates/patient_timeline.html`](src/discern/ui/templates/patient_timeline.html): Muller plot clonal ribbon visualization.
* [`src/discern/ui/templates/performance.html`](src/discern/ui/templates/performance.html): Scientific audit page rendering all E-tables.
* [`src/discern/ui/static/app.css`](src/discern/ui/static/app.css): Lightweight, fully local, accessible CSS styling.
* [`src/discern/ui/static/app.js`](src/discern/ui/static/app.js): Vanilla JavaScript for tabs, timeline scrubber, and interactive controls.

### Evaluation & Verification (`eval/`)
* [`eval/verify.py`](eval/verify.py): Acceptance test validator running 14 golden checks.
* [`eval/build_production.py`](eval/build_production.py): Script that builds the production catalogue and model artifacts.
* [`eval/panel_audit.py`](eval/panel_audit.py): Audits theoretical sensitivity of sequencing panels against the 866-resistant cohort.
* [`eval/e01_prep.py`](eval/e01_prep.py) to [`eval/e23_novel_pipeline.py`](eval/e23_novel_pipeline.py): The 23 individual scientific experiments validating every metric.

### Test Suite (`tests/`)
* [`tests/test_archive_checks.py`](tests/test_archive_checks.py): Tests grammar parsing, fold disjointness, and archived metrics.
* [`tests/test_backend_hardening.py`](tests/test_backend_hardening.py): Tests upload security, error handling, and predictive values.
* [`tests/test_panels.py`](tests/test_panels.py): Tests panel blind spots, LOD warnings, and homopolymer checks.
* [`tests/test_label_shift.py`](tests/test_label_shift.py): Tests MLLS prevalence estimation across simulated shifts.
* [`tests/test_longitudinal.py`](tests/test_longitudinal.py): Tests date parsing, rising subclones, and polyclonal competition.
* [`tests/test_importers_and_reflex.py`](tests/test_importers_and_reflex.py): Tests TB-Profiler, Deeplex, and AmPORE ingestion and reflex codes.
* [`tests/test_batch_recurrence.py`](tests/test_batch_recurrence.py): Tests WGS follow-up and quality-control alerts.
* [`tests/test_novelty_biophysics_and_epistasis.py`](tests/test_novelty_biophysics_and_epistasis.py): Tests epistatic veto and Grantham biophysics.
* [`tests/test_verify.py`](tests/test_verify.py): Meta-tests ensuring the verification harness itself never allows broken acceptance criteria.

---

### Quick Start: Running Discern Locally

```powershell
# 1. Activate conda environment
conda activate discern

# 2. Run the 14 Golden Acceptance Tests
python eval/verify.py

# 3. Run the automated tests
python -m pytest tests/ -q

# 4. Start the offline web application
python -m uvicorn discern.api.app:app --app-dir src --host 127.0.0.1 --port 8000 --no-access-log
```

Open `http://127.0.0.1:8000/` in your browser.
