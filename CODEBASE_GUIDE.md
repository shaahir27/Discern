# Discern: The Complete Beginner-to-Expert Codebase Guide

> **Research use only. For evaluation and education. Clinical decisions require independent laboratory confirmation.**

---

## Welcome to Discern!

If you have never worked with bacterial genomics, machine learning, or clinical tuberculosis research—take a deep breath. You are in the right place.

This guide was written to be read like a personal conversation between two colleagues. No dense jargon without an immediate, clear explanation. No hand-waving or magical black boxes. By the time you reach the end of this document, you will understand:
1. **The biological battle** between tuberculosis bacteria and the world's most critical new antibiotic.
2. **The machine learning and statistical foundations** that allow us to evaluate resistance risk with mathematical rigor.
3. **Every single metric**—what it means, how it is calculated, what was actually tested, and why it matters in research.
4. **The entire codebase architecture**, from raw sequencing ingestion to the final web interface.

Let's dive in!

---

## Table of Contents
1. [The Big Picture: What Problem Are We Solving?](#1-the-big-picture-what-problem-are-we-solving)
2. [The Biology Made Simple: The TB Cell and Bedaquiline](#2-the-biology-made-simple-the-tb-cell-and-bedaquiline)
3. [The Core Philosophy: Why Discern is Built Differently](#3-the-core-philosophy-why-discern-is-built-differently)
4. [High-Level Architecture: The Two-Tiered Brain](#4-high-level-architecture-the-two-tiered-brain)
5. [The Metrics Masterclass: Every Number Explained and Validated](#5-the-metrics-masterclass-every-number-explained-and-validated)
6. [Data Ingestion, Grammar, and Cohort Harmonization](#6-data-ingestion-grammar-and-cohort-harmonization)
7. [The Core Engine: Module-by-Module Walkthrough](#7-the-core-engine-module-by-module-walkthrough)
8. [The Five Cutting-Edge Clinical Research Frontiers](#8-the-five-cutting-edge-clinical-research-frontiers)
9. [Production Hardening, Privacy, and Transparent UI](#9-production-hardening-privacy-and-transparent-ui)
10. [The Evaluation Harness & "Golden Mode" Invariance](#10-the-evaluation-harness--golden-mode-invariance)
11. [Walkthrough: A Sample's Complete Journey Through Discern](#11-walkthrough-a-samples-complete-journey-through-discern)
12. [Complete File-by-File Directory Reference](#12-complete-file-by-file-directory-reference)
13. [How to Run, Test, and Verify](#13-how-to-run-test-and-verify)

---

## 1. The Big Picture: What Problem Are We Solving?

### What is Tuberculosis (TB)?
Tuberculosis is a bacterial infection caused by a rod-shaped organism called *Mycobacterium tuberculosis* ($M.tb$). It primarily attacks the lungs. When an infected individual coughs or speaks, tiny microscopic droplets spread the bacteria through the air. TB remains one of the deadliest infectious diseases on Earth, claiming over 1.3 million lives every year.

### What is Bedaquiline (BDQ)?
For decades, treating drug-resistant TB was a grueling ordeal: 18 to 24 months of taking up to 20 toxic pills daily, alongside painful daily injections that frequently caused permanent hearing loss and psychiatric illness.

Then came **Bedaquiline**. Approved as the first new class of TB medicine in over four decades, bedaquiline fundamentally transformed care. Today, it serves as the essential backbone of the modern **BPaLM / BPaL regimens** (Bedaquiline, Pretomanid, Linezolid, and Moxifloxacin) recommended by the World Health Organization (WHO). Instead of two years of painful therapy, many patients can now be treated in **6 months of all-oral therapy**.

### The Diagnostic Challenge: Does the Drug Work for this Strain?
Bedaquiline only cures a patient if the specific strain of bacteria infecting their lungs is susceptible to it. If the bacteria are resistant, the drug will fail. 

### The Traditional Method: Phenotypic Testing (pDST)
Traditionally, clinics send sputum (phlegm) to a biosafety laboratory. The lab attempts to grow the bacteria in liquid growth media or on agar plates containing bedaquiline:
- If the bacteria grow, the sample is classified as **Resistant**.
- If the bacteria are killed, the sample is classified as **Susceptible**.

**The Catch:** *Mycobacterium tuberculosis* grows exceptionally slowly. It can take **2 to 8 weeks** to obtain a phenotypic result! In many high-burden settings, central reference laboratories are far away, creating prolonged backlogs during which patients wait for definitive confirmation.

### The Modern Alternative: Targeted DNA Sequencing (tNGS)
Instead of waiting weeks to grow bacteria, modern labs can extract bacterial DNA and read its genetic letters in **24 to 48 hours** using benchtop sequencers. 

Sequencing yields a list of "spelling changes" (genetic mutations) in the bacterial genome. For example:
> *"Position 141 of gene Rv0678 has an extra Cytosine inserted."*

**The Dilemma:** What does that spelling change actually mean? Is it a dangerous mutation that breaks the drug's mechanism, or is it a harmless natural variation?

Existing reference tools (such as the WHO Mutation Catalogue v2) classify well-known mutations. But for rare, newly observed, or contested mutations, they output **"Uncertain ($U$)"**. In clinical practice, up to **15% to 25% of patient isolates** carry mutations that fall into this uncertain category.

### What Discern Does
Discern is an **open-source, offline, clinically calibrated triage engine for research and evaluation**:
1. It ingests the raw list of mutations found in an isolate.
2. It evaluates the sample through an explainable, **two-tiered system**:
   - **Tier A (Rule-Based Catalogue):** Provides a conservative genetic classification:
     - **Resistant ($R$)**: Flags an isolate with established genetic evidence of resistance, referring the case for laboratory confirmation.
     - **Susceptible ($S$)**: Indicates that no catalogued resistance variants were identified in interrogated genes. (Clinical decisions still require laboratory confirmation, as resistance from unsequenced loci or unknown mechanisms cannot be ruled out).
     - **Send for lab test ($U$)**: Identifies variants whose phenotypic consequence is uncertain or disputed, routing the isolate for prioritized physical laboratory testing.
   - **Tier B (Statistical Machine Learning):** Computes a calibrated probability of resistance.
3. Crucially, **Discern never leaves unresolved samples in an unordered pile**. It ranks the unresolved ($U$) queue by calibrated risk so laboratories can prioritize which physical cultures to confirm first.

```
       Patient Sputum Sample
                 │
                 ▼
       DNA Sequencing (tNGS / WGS)
                 │
                 ▼
   List of Mutations (e.g., Rv0678@141_ins_c)
                 │
                 ▼
      ┌─────────────────────┐
      │       DISCERN       │
      │  Triage Engine & ML │
      └──────────┬──────────┘
                 │
      ┌──────────┼──────────┐
      ▼          ▼          ▼
  Resistant  Susceptible  Unresolved
   (Flag)    (No Var.)   (Ranked for pDST)
```

---

## 2. The Biology Made Simple: The TB Cell and Bedaquiline

To understand how Discern's code works, you only need to understand five key genes in the TB genome. Think of them as mechanical components inside the bacterial cell:

```
┌────────────────────────────────────────────────────────┐
│                   THE BACTERIAL CELL                   │
│                                                        │
│   ┌───────────────────────────┐                        │
│   │   atpE: The Power Turbine │ ◄─── BEDAQUILINE       │
│   │   (ATP Synthase Motor)    │      (Jams the rotor)  │
│   └───────────────────────────┘                        │
│                                                        │
│   ┌───────────────────────────┐                        │
│   │   Rv0678: Security Guard  │                        │
│   │   (Transcriptional Switch)│                        │
│   └─────────────┬─────────────┘                        │
│                 │ (If broken, guard falls off switch)  │
│                 ▼                                      │
│   ┌───────────────────────────┐                        │
│   │   mmpL5 / mmpS5: Sump Pump│ ───► Pumps bedaquiline │
│   │   (Efflux Pump Machinery) │      OUT of the cell!  │
│   └───────────────────────────┘                        │
│                                                        │
│   ┌───────────────────────────┐                        │
│   │   pepQ: Cytoplasmic Peptidase                      │
│   │   (Amino Acid Recycler)   │                        │
│   └───────────────────────────┘                        │
└────────────────────────────────────────────────────────┘
```

### 1. `atpE` — The Power Turbine (Direct Drug Target)
- **What it does:** Produces subunit $c$ of the ATP synthase enzyme. ATP synthase is a rotating molecular motor that manufactures ATP (cellular energy). Without ATP, the bacteria starve and die.
- **How Bedaquiline Works:** Bedaquiline physically binds into the rotor of the ATP synthase turbine, jamming the gears.
- **How Resistance Occurs:** If the DNA of `atpE` mutates at the exact drug-binding pocket (such as codon changes `A63P` or `I66M`), bedaquiline can no longer bind. This causes **high-level, efflux-independent resistance**.

### 2. `Rv0678` — The Security Guard (Transcriptional Repressor)
- **What it does:** `Rv0678` produces a homodimeric protein that binds directly to the bacterial chromosome. It acts like a security guard holding an electrical switch in the **OFF** position—specifically repressing the production of efflux pumps.
- **How Resistance Occurs:** If `Rv0678` is damaged—by a premature stop codon (`!`), an insertion/deletion that disrupts the reading frame (`indel`), or an amino acid substitution in its DNA-binding domain—the guard falls off the DNA.

### 3. `mmpL5` and `mmpS5` — The Sump Pump (Efflux System)
- **What it does:** A transmembrane transporter that pumps lipids across the bacterial cell wall.
- **How Resistance Occurs:** When the security guard (`Rv0678`) is broken, the cell over-expresses `mmpL5/mmpS5` pumps. These pumps catch bedaquiline inside the cytoplasm and pump it **out** into the extracellular space before it can reach the ATP synthase turbine. This causes low-to-moderate minimum inhibitory concentration (MIC) elevation.

### 4. `pepQ` — The Cytoplasmic Peptidase
- **What it does:** A cytoplasmic metallopeptidase that hydrolyzes proline dipeptides.
- **How Resistance Occurs:** Mutations in `pepQ` cause low-level resistance through cellular stress pathways, often acting as a secondary contributor alongside other mechanisms.

---

### The Biological Interaction: Epistasis (The Pump Veto)

Here is where naive software tools fail, but Discern's mechanistic reasoning succeeds:
> **Question:** What happens if the security guard (`Rv0678`) is broken (signaling the cell to build pumps), BUT the pump gene (`mmpL5`) is ALSO mutated and inactivated?
> 
> **Answer:** Even though the signal says "turn on the pumps," the physical pump motor is broken! The cell cannot pump out the drug. Bedaquiline enters unimpeded, and the bacteria remains **susceptible**!

This biological phenomenon—where a mutation in one gene alters or cancels the effect of a mutation in another—is called **negative epistasis**.

Standard software inspects mutations in isolation. It detects a broken `Rv0678`, flags the sample as "Resistant", and raises an alarm. Discern inspects both loci together. When an inactivating loss-of-function mutation in *mmpL5* is detected, Discern executes an **Epistatic Veto**, downgrading the call to an unresolved referral ($U$) to prevent false resistance alarms.

#### The Allele-Frequency Linkage Condition
Under the **WHO Mutation Catalogue 3rd Edition (2026)** and **TB-Profiler v6.6.5**, there is a crucial requirement:
If the `Rv0678` mutation is present in **95%** of sequencing reads, but the broken pump mutation (`mmpL5`) is present in only **10%**, the pump mutation cannot protect the dominant population. They likely exist in different subclonal populations!

Discern checks this allele-frequency relationship:
$$\text{FRS}(mmpL5) \ge \text{FRS}(Rv0678) - 0.20$$
The pump mutation must be present at an allele fraction consistent with co-occurring in the same bacterial lineage. (Note that short-read sequencing cannot definitively prove physical linkage on the same single chromosome; this remains an important biological caveat).

---

## 3. The Core Philosophy: Why Discern is Built Differently

Discern was engineered around four strict scientific and ethical principles:

### Principle 1: Never Invent or Round Numbers
Every single metric, sensitivity percentage, specificity bound, and confidence interval displayed in the web UI, model card, and documentation is **read directly from programmatic evaluation files** (`results/*.json`). No developer is permitted to hand-type or round up an evaluation number.

### Principle 2: Patient-Grouped Cross-Validation (No Data Leakage)
A single tuberculosis patient might provide multiple sputum samples across several hospital visits. If Sample 1 from Patient $A$ is placed in the training set and Sample 2 from Patient $A$ is in the test set, an algorithm can easily memorize Patient $A$'s unique genetic quirks. 
Discern enforces strict **patient-grouped splitting (`SUBJID`)**: all samples from a given human being are kept together in either the training set or the test set—never split across both!

### Principle 3: "Admit When You Don't Know"
A dangerous algorithm provides an overconfident guess on everything. In medicine, a wrong guess can lead to patient harm. Discern uses a **ternary triage policy**:
- It outputs $R$ or $S$ only when supported by statistical catalogue evidence.
- Otherwise, it assigns **Unresolved ($U$)** and refers the sample for physical laboratory confirmation.

### Principle 4: 100% Offline by Design (Zero CDNs, Zero Internet)
Discern is built to operate within regional hospital laboratories, often inside air-gapped hospital intranets with no external connectivity:
- Zero Google Fonts or external web assets.
- Zero CDNs (no remote Bootstrap or Tailwind).
- Zero tracking scripts or analytics telemetry.
- **Request bodies are never logged** to protect patient confidentiality.

---

## 4. High-Level Architecture: The Two-Tiered Brain

Discern uses a two-tiered decision architecture where deterministic rules and statistical machine learning work together transparently:

```
┌────────────────────────────────────────────────────────┐
│                      INPUT DATA                        │
│   Mutations (GARC syntax) + Sequencing Assay Info      │
└───────────────────────────┬────────────────────────────┘
                            │
                            ▼
┌────────────────────────────────────────────────────────┐
│           TIER A: Deterministic Catalogue              │
│   Two-tailed Binomial test rules built with catomatic  │
│   Evaluates established WHO & curated mutations        │
└─────────────┬────────────────────────────┬─────────────┘
              │ Flags R or S               │ Flags U (Unresolved)
              ▼                            ▼
┌───────────────────────────┐ ┌──────────────────────────┐
│   EPISTATIC VETO FILTER   │ │   TIER B: Calibrated ML  │
│ Checks if mmpL5 pump died │ │ L2 Logistic Regression   │
│ If pump broken -> Referral│ │ Continuous risk score    │
└─────────────┬─────────────┘ └────────────┬─────────────┘
              │                            │
              │                            ▼
              │               ┌──────────────────────────┐
              │               │      Platt Scaling       │
              │               │ Inner out-of-fold fit    │
              │               └────────────┬─────────────┘
              │                            │
              ▼                            ▼
┌────────────────────────────────────────────────────────┐
│            BAYESIAN PREVALENCE ADJUSTMENT              │
│    Transforms probability to match user clinic rate    │
│    Computes True Positive (PPV) & Negative (NPV) values│
└───────────────────────────┬────────────────────────────┘
                            │
                            ▼
┌────────────────────────────────────────────────────────┐
│                    FINAL CLINICAL OUTPUT               │
│  Triage Call (R/S/U) + Evidence Cards + Priority Queue │
└────────────────────────────────────────────────────────┘
```

### Tier A: The Deterministic Catalogue Baseline
- Acts as a conservative, rule-based algorithmic baseline.
- Evaluates non-synonymous mutations using a statistical binomial catalogue built on 11,868 isolates using the `catomatic` and `piezo` libraries.
- If a sample carries a variant with statistically verified resistance association, Tier A flags **Resistant ($R$)**.
- If a sample carries only benign variants, Tier A classifies it as **Susceptible ($S$)**.
- If a mutation is rare, unstudied, or ambiguous, Tier A outputs **Unresolved ($U$)**.

### Tier B: The Calibrated Machine Learning Model
- An $L_2$-regularized Logistic Regression model ($C=0.3$) trained on gene-level aggregate features (frameshift indicators, stop codons, DNA-binding domain alterations) and recurrent variant counts.
- **The Core Constraint:** Tier B **never overrides Tier A to call resistance on its own**.
- Instead, Tier B:
  1. Computes a continuous, calibrated probability of resistance.
  2. Ranks unresolved ($U$) cases so laboratories can prioritize which physical cultures to test first.
  3. Raises disagreement warnings if Tier A and Tier B conflict.

---

## 5. The Metrics Masterclass: Every Number Explained and Validated

In genomic machine learning, high-level metrics are often quoted without explaining what was actually measured. Below is an honest breakdown of Discern's validated metrics:

### 1. Sensitivity: Firm vs. All-Case
- **What is Sensitivity?** *Out of 100 truly resistant patient isolates confirmed by laboratory testing, how many did the system successfully flag as resistant?*
- **Firm Sensitivity: 79.24% (95% CI: 75.96%–82.43%)**  
  This number (from experiment [E03](results/approaches/catalogue/e03_s1_catalogue_cv/results.json)) measures performance **only among firm calls**—that is, the 94.43% of cases where the catalogue reached a definite ($R$ or $S$) conclusion.
- **All-Case Sensitivity (U as S): 60.39% (95% CI: 57.06%–63.63%)**  
  If you look across *all* 866 laboratory-confirmed resistant isolates in the reference cohort and count any unresolved ($U$) call as not flagged ($S$), the catalogue flagged **60.39%** as resistant. (Conversely, if all $U$ calls were aggressively treated as resistant, sensitivity would be 84.18%, but specificity would drop to 94.28%).
- **Why isn't sensitivity 100%? (The "Invisible 16%"):** In the reference cohort, **15.8% of phenotypically resistant strains carry zero mutations** in any of the known bedaquiline-relevant genes (*Rv0678*, *atpE*, *pepQ*). Because their resistance is mediated by unmapped biological mechanisms, no targeted DNA tool can achieve 100% sensitivity on this data.

### 2. Specificity: 98.35% (95% CI: 98.10%–98.58%)
- **What is Specificity?** *Out of 100 drug-susceptible patient isolates, how often did the catalogue correctly avoid raising a false resistance alarm?*
- **Our Number:** The catalogue achieves **98.35%** specificity among firm calls.

### 3. Discrimination: Internal Cross-Validation vs. External Validation
- **What is ROC-AUC?** The Area Under the Receiver Operating Characteristic curve measures ranking quality. If you pick one random resistant strain and one random susceptible strain, the AUC represents the probability that the model assigns a higher resistance score to the resistant strain.
- **Internal Cross-Validation (S1): 0.9085 (95% CI: 0.8944–0.9230)**  
  Evaluated across 5 patient-grouped folds on the reference cohort (experiment [E04](results/models/lr_model/e04_s1_patient_cv/results.json)).
- **External Validation (KM Cohort): 0.8682**  
  Evaluated on an independent, non-overlapping external cohort of 792 isolates (35 resistant) from Zenodo dataset 15680920 (experiment [ext_km_validate.py](results/approaches/external_km/e_km_validation/REPORT.md)). The drop from 0.9085 to 0.8682 reflects real-world domain transfer across different laboratories and lineage distributions.
- **AUPRC / Average Precision: 0.6836 (95% CI: 0.6471–0.7194)**  
  Because true resistance is relatively rare (7.3% baseline prevalence), Area Under the Precision-Recall Curve provides a more stringent evaluation than ROC-AUC by focusing on positive predictive performance.

### 4. Positive & Negative Predictive Values (PPV & NPV)
- **What are PPV and NPV?**  
  - **PPV:** *When the tool outputs a Resistant flag, what is the probability the isolate is truly resistant?*  
  - **NPV:** *When the tool outputs a Susceptible classification, what is the probability the isolate is truly susceptible?*
- **Prevalence-Dependent Estimates:** PPV and NPV are **not universal constants**. They are model-derived estimates calculated via Bayes' theorem for specific assumed prevalence scenarios:
  - In a clinic with **1% resistance prevalence**: Estimated PPV is **14%** (Model) / **33%** (Catalogue). (Because resistance is rare, most flags will be false alarms).
  - In a clinic with **5% resistance prevalence**: Estimated PPV is **47%** (Model) / **72%** (Catalogue).
  - In a clinic with **10% resistance prevalence**: Estimated PPV is **65%** (Model) / **85%** (Catalogue).
- Under these model estimates, the Negative Predictive Value (NPV) remains high across evaluated scenarios (typically exceeding 98% to 99%), reflecting the high specificity of the catalogue. However, these figures are model estimates, not empirical guarantees for all clinical settings.

### 5. Calibration Metrics: Brier Score & Expected Calibration Error (ECE)
- **Platt Scaling:** Raw logistic regression scores are not calibrated probabilities. Nested Platt scaling fits a sigmoid curve on inner out-of-fold log-odds to output probabilities that reflect empirical frequencies.
- In experiment [E10](results/approaches/prevalence_correction/e10_prevalence_shift/REPORT.md), when adjusting across divergent cohorts, Bayesian prevalence correction reduced the Brier error from **0.102 down to 0.017**, and lowered the Expected Calibration Error (ECE) from **0.288 down to 0.013**.

---

## 6. Data Ingestion, Grammar, and Cohort Harmonization

### The Mutation Grammar (`src/discern/grammar.py`)
In raw sequencing files, genetic variants are reported in diverse formats. Discern standardizes variants into GARC format: `GENE@MUTATION`.

| Mutation String | Classification | Biological Interpretation |
|---|---|---|
| `Rv0678@R50Q` | Missense | Amino acid 50 Arginine (R) mutated to Glutamine (Q). |
| `Rv0678@Q131!` | Nonsense (Stop) | Amino acid 131 turned into a premature STOP (`!`), truncating the protein. |
| `Rv0678@141_ins_c` | Insertion (Indel) | A Cytosine inserted at nucleotide 141, shifting the reading frame. |
| `Rv0678@193_del_g` | Deletion (Indel) | A Guanine deleted at nucleotide 193. |
| `Rv0678@c-11a` | Promoter | In the promoter region 11 bases before the start codon, Cytosine became Adenine. |
| `Rv0678@S68S` | Synonymous | Amino acid 68 remained Serine. Silent mutation. |
| `Rv0678@S68Z` | Null / Uncallable | Position 68 could not be resolved by the sequencer. Triggers a warning. |

### The Reference Cohorts (`src/discern/prep.py`)
Discern was trained and evaluated on **11,868 clinical isolates** pooled from two major open-access public studies:
1. **CRyPTIC Consortium:** 10,704 high-quality isolates collected across 23 international hospital laboratories using 14-day microtiter plates (UKMYC). It represents broad global surveillance (83 resistant strains).
2. **NICD (National Institute for Communicable Diseases, South Africa):** 1,164 isolates tested using liquid MGIT960 culture. This cohort was heavily enriched for drug-resistant cases (783 resistant strains).

Pooling these cohorts yields **866 resistant strains** and **11,002 susceptible strains**, providing a robust, diverse benchmark for bedaquiline resistance evaluation.

---

## 7. The Core Engine: Module-by-Module Walkthrough

Let's inspect the core engine inside [`src/discern/`](src/discern/):

### [`prep.py`](src/discern/prep.py) — Raw Data Harmonization
- Filters CRyPTIC to retain only records with `QUALITY == "HIGH"`.
- Merges CRyPTIC and NICD, resolving duplicate sample IDs by prioritizing resistant records and validated MIC measurements.
- Exports compressed Parquet tables: `data/processed/samples.parquet` and `data/processed/mutations.parquet`.

### [`splits.py`](src/discern/splits.py) — Patient-Grouped Cross-Validation
- Generates 4 evaluation split schemes:
  - **S1:** 5-fold patient-grouped (`GroupKFold` on `SUBJID`).
  - **S2:** Leave-One-Site-Out (LOSO) across all 23 international hospital sites.
  - **S3:** Cross-Source shift (train on CRyPTIC, test on NICD, and vice versa).
  - **S4:** Variant-set grouping (clustering by unique mutation combinations to evaluate novel variants).

### [`catalogue.py`](src/discern/catalogue.py) — Tier A Rule Builder
- Interfaces with `catomatic` and `piezo` bioinformatics tools.
- Builds two-tailed Binomial test catalogues with strict unlocking rules and wildcards, generating `models/lr_model/catalogue_bdq.csv`.

### [`features.py`](src/discern/features.py) — Numerical Representation
Transforms genomic text into a numerical feature matrix:
1. **Recurrent Vocabulary:** Any mutation observed in $\ge 3$ training samples receives its own indicator column.
2. **Biological Aggregates:** Maximum allele fraction for *Rv0678* missense, indels, stop codons, promoter changes, DNA-binding domain alterations (codons 34–99), *atpE* target mutations, and *mmpL5* loss-of-function flags.

### [`model.py`](src/discern/model.py) — Tier B Classifier & Calibration
- Fits $L_2$-regularized Logistic Regression ($C=0.3$).
- Executes **nested out-of-fold Platt scaling**: fits an inner logistic regression on out-of-fold log-odds to output calibrated probabilities.
- Identifies Youden threshold $t_{\text{bin}}$, low threshold $lo$ (NPV $\ge 0.985$), and high threshold $hi$ (PPV $\ge 0.800$).

### [`prevalence.py`](src/discern/prevalence.py) — Bayesian Prevalence Adjustment
Adjusts predicted probabilities from the reference training prevalence ($\pi_{\text{train}} \approx 7.3\%$) to a user-specified clinic prevalence ($\pi_{\text{clinic}}$) using Bayes' rule in odds formulation:
$$\text{odds}(p) = \frac{p}{1 - p}$$
$$\text{odds}' = \frac{\text{odds}(p)}{\text{odds}(\pi_{\text{train}})} \times \text{odds}(\pi_{\text{clinic}})$$
$$p_{\text{adjusted}} = \frac{\text{odds}'}{1 + \text{odds}'}$$

### [`mechanism.py`](src/discern/mechanism.py) — The Epistatic Veto Engine
- Parses the *Rv0678* coding sequence (165 amino acids).
- Classifies variants into DNA-binding domain ($\le 99$), mid-region (100–160), and late-region ($> 160$).
- Enforces allele-frequency linkage: $\text{FRS}(pump) \ge \text{FRS}(Rv0678) - 0.20$.
- Downgrades eligible *Rv0678* false-resistant calls to referral ($U$) when an inactivating pump mutation co-occurs.

### [`structure.py`](src/discern/structure.py) — 3D Biophysical Feature Engine
- Parses 3D atomic coordinates from `data/raw/Rv0678_WT.pdb`.
- Calculates minimum atomic distance from each amino acid to the DNA double helix.
- Computes residue burial (count of alpha-carbon atoms within a 10-Ångström sphere).
- Computes Grantham physicochemical distance values to quantify amino acid radicality.

---

## 8. The Five Cutting-Edge Clinical Research Frontiers

Discern incorporates 5 specialized modules designed to address real-world laboratory workflows:

```
┌────────────────────────────────────────────────────────────────────────┐
│                   DISCERN'S 5 CLINICAL RESEARCH FRONTIERS              │
├────────────────────────────────────────────────────────────────────────┤
│ 1. PANEL AUDITOR       ──► Audits sequencing kit blind spots & poly-G  │
│ 2. MLLS PREVALENCE     ──► Estimates batch prevalence under shift model│
│ 3. LONGITUDINAL RIBBON ──► Tracks within-patient subclone trajectories │
│ 4. UNIVERSAL IMPORTERS ──► Ingests native TB-Profiler / Deeplex / AmPORE│
│ 5. SURVEILLANCE AUDIT  ──► Flags batch recurrence vs technical artifacts│
└────────────────────────────────────────────────────────────────────────┘
```

### Frontier 1: Sequencing-Assay & tNGS Panel Blind-Spot Auditor
* **Source:** [`src/discern/panels.py`](src/discern/panels.py) & [`eval/panel_audit.py`](eval/panel_audit.py)
* **The Clinical Issue:** Different sequencing kits sequence different genomic targets. Standard **Deeplex Myc-TB** and **ICMR-NIRT MinION** kits sequence *Rv0678*, but omit *mmpL5* and *mmpS5*. If software assumes unsequenced loci are wild-type, the epistatic veto is silently blinded. Furthermore, Nanopore platforms can encounter homopolymer errors (runs of identical bases like `GGGGGG` at positions 192–198 of *Rv0678*), occasionally producing false indel calls.
* **What Discern Does:**
  - Profiles 5 platforms: WGS, Deeplex Myc-TB, Deeplex Myc-TB XL, ONT AmPORE-TB, and ICMR-NIRT MinION.
  - Flags `EPISTATIC_VETO_BLIND` when pump loci are unsequenced.
  - Flags `UNSEQUENCED_CORE_LOCI` when target loci (*atpE* / *pepQ*) are omitted.
  - Flags `HOMOPOLYMER_ARTIFACT_RISK` when low-VAF indels occur in homopolymer tracts on Nanopore platforms.
  - Provides a 5-gene coverage matrix at `/panels` and an interactive visual strip on `/result`.

### Frontier 2: Unsupervised Batch Prevalence Estimation (MLLS / Saerens EM)
* **Source:** [`src/discern/label_shift.py`](src/discern/label_shift.py)
* **The Clinical Issue:** Discern's prevalence control asks for the local resistance prevalence. Regional labs without routine phenotypic culture rarely know this parameter.
* **What Discern Does:**
  - Implements **Maximum Likelihood Label Shift (Saerens Expectation-Maximization)**.
  - When a lab uploads an unlabelled batch of sequenced isolates ($N \ge 30$), Discern **estimates** the local prevalence $\hat{\pi}$ from the distribution of model scores, assuming conditional score distributions $P(X|Y)$ transfer from the reference cohort.
  - Computes a 100-replicate patient-bootstrap 95% confidence interval and automatically deduplicates serial isolates from the same patient.
  - Displays an interactive, inline SVG curve of Positive Predictive Value (PPV) across 0.5%–25% prevalence on `/queue`.

### Frontier 3: Longitudinal Serial-Isolate Tracking & Clonal Ribbons
* **Source:** [`src/discern/longitudinal.py`](src/discern/longitudinal.py)
* **The Clinical Issue:** TB therapy lasts months. Resistance often emerges gradually as a minority subclone (e.g. 12% frequency) that outcompetes sensitive bacteria under drug selection until becoming dominant (e.g. 95% frequency). Static single-timepoint tools miss these dynamics.
* **What Discern Does:**
  - Ingests optional `patient_id` and `collection_date` columns.
  - Classifies clonal dynamics: `DE_NOVO_EMERGENCE`, `RISING_SUBCLONE`, `FIXED_RESISTANCE`, `POLYCLONAL_COMPETITION`, and `TRANSIENT_SUBCLONE`.
  - Serves an interactive **Muller Plot** (clonal ribbon timeline) at `/patient/{patient_id}` with a time scrubber and visit inspector. (Demonstrated with [`demo/longitudinal_example.csv`](demo/longitudinal_example.csv)).

### Frontier 4: Universal Lab Importers & Mechanism-Stratified Reflex Guidance
* **Source:** [`src/discern/importers.py`](src/discern/importers.py) & [`src/discern/decision.py`](src/discern/decision.py)
* **The Clinical Issue:** Real sequencers output varied native formats: TB-Profiler JSON files, Deeplex Myc-TB CSV exports, or AmPORE-TB TSV spreadsheets using HGVS mutation notation (`p.Gln131*` instead of `Q131!`, `c.-11C>A` instead of `c-11a`).
* **What Discern Does:**
  - Ingests TB-Profiler, Deeplex, and AmPORE files automatically.
  - Normalizes HGVS syntax into GARC mutation notation.
  - Generates standardized **Reason Codes** (`RC_RV0678_LOF_STOP`, `RC_ATPE_TARGET_MUTATION`, etc.).
  - Maps mutations to a 165-amino-acid structural domain track on `/result`.
  - Displays mechanism-stratified clinical reflex notes: highlights clofazimine (CFZ) cross-resistance risk and prioritizes samples for phenotypic confirmation.

### Frontier 5: Lightweight Batch-Recurrence & Contamination Surveillance
* **Source:** [`src/discern/surveillance.py`](src/discern/surveillance.py)
* **The Clinical Issue:** If the same rare bedaquiline-resistant mutation appears across multiple patients in a single sequencing run:
  1. *Scenario A:* It could reflect an active transmission cluster in the community.
  2. *Scenario B:* It could reflect cross-contamination or an amplicon sequencing artifact.
* **What Discern Does:**
  - Deduplicates patient records.
  - If an ultra-rare mutation ($\le 5$ historical cohort occurrences) appears across multiple distinct patients, it triggers a `WGS_SURVEILLANCE_ESCALATION` alert.
  - If the recurring mutation is a known homopolymer indel, it flags `CONTAMINATION_OR_ARTIFACT_RISK` rather than alarming clinicians prematurely.

---

## 9. Production Hardening, Privacy, and Transparent UI

Discern incorporates robust production safeguards:

### Bounded In-Memory Session Storage
- Live uploaded batches are maintained in a bounded in-memory cache (`UploadedBatch` in [`src/discern/api/app.py`](src/discern/api/app.py)).
- Scoped to the user's browser via a secure `discern_upload` cookie.
- Live links carry a unique batch identifier (`?batch=...`). A live upload from one user **never leaks to another session** and **never shadows pre-computed replay records**.
- Users can export their batch as a private CSV via `/queue?download=true`.

### Transparent Evidence Display
On [`src/discern/ui/templates/result.html`](src/discern/ui/templates/result.html):
1. **Explicit Referral Reasons:** When an isolate is referred to the laboratory ($U$), Discern explicitly identifies the reason:
   - Was it an unresolved mutation?
   - Did the epistatic veto downgrade an *Rv0678* finding?
   - Did the catalogue parser encounter an unsupported format?
   - Or is it an archived record where the historical reason was not saved?
2. **Cohort Evidence First:** Each variant card displays historical carrier counts broken down by study source (`CRyPTIC` vs `NICD`) *before* showing the pooled summary, because NICD was enriched for multidrug-resistant cases.
3. **The Phasing Notice:** When an isolate carries multiple mutations, Discern displays an explicit note:
   > *"This file does not show whether these variants are on the same genome."*
   Short-read bulk sequencing cannot establish whether two mutations co-occur in the same bacterial genome (*cis*) or reside in separate sub-populations (*trans*).
4. **Comprehension Protocol:** The interface's clarity is governed by a formal protocol ([`docs/RESULT_REVIEW_PROTOCOL.md`](docs/RESULT_REVIEW_PROTOCOL.md)) designed to test whether users correctly distinguish model scores from laboratory results.

---

## 10. The Evaluation Harness & "Golden Mode" Invariance

In medical software, scientific reproducibility is non-negotiable.

### What is the Golden Mode (`D003`)?
The production model artifacts:
- `models/lr_model/bdq_lr.joblib`
- `models/lr_model/catalogue_bdq.csv`
- `models/lr_model/metadata.json`

These files are **cryptographically hashed and frozen**. No code modification is ever permitted to alter or retrain these files.

### The Acceptance Verifier (`eval/verify.py`)
Run `python eval/verify.py` anytime. It tests **14 pre-registered scientific gates**:

| Gate ID | Experiment | What It Proves | Acceptance Tolerance |
|---|---|---|---|
| **E01** | Data Prep | Validates exact cohort sample counts (11,868 samples, 866 resistant). | Exact match |
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
| **E13** | PPV / NPV Tables | Audits predictive values across prevalence rates. | Exact match |

---

## 11. Walkthrough: A Sample's Complete Journey Through Discern

Below is an illustrative walkthrough of how Discern processes a sample:

```
Step 1: Ingestion
A user uploads an Oxford Nanopore AmPORE-TB file:
  Sample ID: PATIENT_42
  Gene: Rv0678
  Mutation: p.Gln131*
  Frequency: 92% (FRS = 0.92)
  Assay: ampore_tb

Step 2: Universal Import & Normalization (src/discern/importers.py)
Discern detects AmPORE format. Converts HGVS `p.Gln131*` into GARC `Rv0678@Q131!`.

Step 3: Assay Blind-Spot Audit (src/discern/panels.py)
Panel detected: AmPORE-TB (covers Rv0678 and atpE, but omits mmpL5/mmpS5).
Warning raised: `EPISTATIC_VETO_BLIND` (Pump loci were not sequenced).

Step 4: Tier A Rule Evaluation (src/discern/catalogue.py)
`Q131!` is a known loss-of-function stop mutation in Rv0678.
Catalogue call: RESISTANT (R).

Step 5: Epistatic Veto Check (src/discern/mechanism.py)
Did mmpL5 carry an inactivating mutation?
-> No mmpL5 mutation was sequenced. Veto cannot fire. Call remains R.

Step 6: Tier B Statistical Scoring & Calibration (src/discern/model.py)
Aggregates computed: `rv_stop = 0.92`, `rv_any = 0.92`.
Raw logistic regression predicted probability: 0.861.
Platt-calibrated probability at reference training prevalence (pi_tr = 0.073): 0.892.

Step 7: Bayesian Prevalence Adjustment (src/discern/prevalence.py)
Assumed clinic prevalence selected: 5.0% (pi = 0.05).
Adjusted probability of resistance:
  odds(0.892) = 0.892 / (1 - 0.892) = 8.259
  odds(0.073) = 0.073 / (1 - 0.073) = 0.0787
  odds(0.050) = 0.050 / (1 - 0.050) = 0.0526
  odds' = 8.259 * (0.0526 / 0.0787) = 5.522
  p_adj = 5.522 / (1 + 5.522) = 0.8467 (84.7%)
Estimated Positive Predictive Value (PPV): 72.0%.

Step 8: Output Delivery (result.html)
- Triage Badge: RESISTANT (Red).
- Call Driver: `Rv0678@Q131!` supports Resistant.
- Domain Track: Pins residue 131 in the Rv0678 dimerization domain.
- Clinical Reflex Advice: Warns of potential cross-resistance to Clofazimine (CFZ); recommends phenotypic confirmation.
- Panel Strip: Visual indicator showing *Rv0678* and *atpE* sequenced, *mmpL5/S5* not interrogated.
```

---

## 12. Complete File-by-File Directory Reference

### Root Directory
* [`AGENTS.md`](AGENTS.md): The master technical specification and golden numbers.
* [`config.yaml`](config.yaml): Global system configuration (file paths, random seeds, hyperparameters).
* [`pyproject.toml`](pyproject.toml): Python package build metadata and pytest configuration.
* [`requirements.txt`](requirements.txt): Exact pinned Python dependencies.
* [`Makefile`](Makefile): Developer shortcuts (`make verify`, `make test`, `make app`).
* [`Dockerfile`](Dockerfile): Offline container build configuration.
* [`CHANGELOG.md`](CHANGELOG.md): Detailed chronicle of all features and architectural decisions.

### Core Source Code (`src/discern/`)
* [`src/discern/grammar.py`](src/discern/grammar.py): Mutation string parsing and classification.
* [`src/discern/prep.py`](src/discern/prep.py): Data cleaning, filtering, and merging of CRyPTIC and NICD cohorts.
* [`src/discern/splits.py`](src/discern/splits.py): Patient-grouped cross-validation splitting algorithms.
* [`src/discern/catalogue.py`](src/discern/catalogue.py): Tier A deterministic resistance catalogue builder (`catomatic`/`piezo`).
* [`src/discern/features.py`](src/discern/features.py): Feature extraction (aggregates, vocabulary, interaction terms).
* [`src/discern/model.py`](src/discern/model.py): Tier B $L_2$ Logistic Regression, nested Platt scaling, and threshold search.
* [`src/discern/prevalence.py`](src/discern/prevalence.py): Bayesian prevalence adjustment and predictive value calculation.
* [`src/discern/mechanism.py`](src/discern/mechanism.py): Epistatic veto engine and allele-frequency linkage verification.
* [`src/discern/structure.py`](src/discern/structure.py): 3D PDB structure parser, Grantham matrix, and $\alpha 4$ helix biophysics.
* [`src/discern/panels.py`](src/discern/panels.py): Sequencing panel profiles and blind-spot auditor.
* [`src/discern/label_shift.py`](src/discern/label_shift.py): MLLS (Saerens EM) unsupervised batch prevalence estimation.
* [`src/discern/longitudinal.py`](src/discern/longitudinal.py): Serial-isolate tracking, clonal dynamics, and Muller plot generator.
* [`src/discern/importers.py`](src/discern/importers.py): TB-Profiler, Deeplex, and AmPORE file importers and HGVS normalizer.
* [`src/discern/surveillance.py`](src/discern/surveillance.py): Batch recurrence, outbreak escalation, and artifact alerts.
* [`src/discern/evidence.py`](src/discern/evidence.py): Historical mutation carrier and resistance query store.
* [`src/discern/quality.py`](src/discern/quality.py): Live uploaded CSV data completeness and syntax auditor.
* [`src/discern/decision.py`](src/discern/decision.py): Main runtime prediction pipeline coordinating Tier A, Tier B, and vetoes.
* [`src/discern/research.py`](src/discern/research.py): Loader for pre-calculated golden cross-validation results.
* [`src/discern/cli.py`](src/discern/cli.py): Command-line interface (`discern predict`, `discern batch`, `discern verify`).

### Web Interface & API (`src/discern/api/` & `src/discern/ui/`)
* [`src/discern/api/app.py`](src/discern/api/app.py): FastAPI web application serving all UI routes and API endpoints.
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
* [`tests/test_batch_recurrence.py`](tests/test_batch_recurrence.py): Tests outbreak surveillance escalation vs. homopolymer artifacts.
* [`tests/test_novelty_biophysics_and_epistasis.py`](tests/test_novelty_biophysics_and_epistasis.py): Tests epistatic veto and Grantham biophysics.
* [`tests/test_upload_integration.py`](tests/test_upload_integration.py): Full end-to-end integration tests for web uploads, queues, and timeline views.
* [`tests/test_upload_quality_and_live_scoring.py`](tests/test_upload_quality_and_live_scoring.py): Tests CSV syntax quality grading and live prediction scoring.
* [`tests/test_verify.py`](tests/test_verify.py): Meta-tests ensuring the verification harness itself never allows broken acceptance criteria.

---

## 13. How to Run, Test, and Verify

You can run the entire system locally in under 30 seconds:

```powershell
# 1. Activate the conda environment
conda activate discern

# 2. Run the 14 Golden Acceptance Tests
python eval/verify.py

# 3. Run all 89 automated tests
python -m pytest tests/ -q

# 4. Start the offline web application
python -m uvicorn discern.api.app:app --app-dir src --host 127.0.0.1 --port 8000 --no-access-log
```

Open your web browser and navigate to:
**`http://127.0.0.1:8000/`**

You can now explore single isolates, upload batch runs, inspect sequencing kit blind spots, and review the full scientific performance archive. Congratulations—you now understand every single dimension of Discern!
