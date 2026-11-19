# AGENTS.md: BDQ-Triage build specification

You are the coding agent for **BDQ-Triage**. Read this whole file before writing code. It is the single source of truth.
The folder `reference/` holds working scripts that produced every number in this spec. Treat them as an **executable specification**: port them into the package, keep the outputs identical, and do not "improve" the method.

---

## 0. Mission

Build an open, offline, explainable tool that reads the DNA variants found in a TB germ's bedaquiline-relevant genes and returns **Resistant / Susceptible / Send for lab test**, with evidence, plus an evaluation harness that regenerates every number in section 11.

**What the tool is:** research-use triage for TB labs. Input: a variant table. Output: a call, evidence per variant, a calibrated probability at the lab's resistance prevalence, warnings, and a ranked queue of unresolved samples.

**Non-goals (do not build, do not claim):** a diagnostic test, treatment advice, deep learning, prediction of never-seen variants, cost-saving claims, Indian sensitivity claims, anything needing the internet at runtime.

## 1. Working rules (hard)

1. **Never invent or round up numbers.** Every number in the UI, README and model card must be read from a file written by the evaluation harness (`results/*.json`).
2. **Splits are patient-grouped** (`SUBJID`) unless a section says otherwise. Never split by sample.
3. **Never fit anything on test rows**: not the variant vocabulary, not the calibrator, not the thresholds.
4. **No `class_weight`, no resampling** in the calibrated model. (Section 11, E11 shows balancing does not help and it breaks calibration.)
5. If any acceptance number in section 11 misses its tolerance, **stop** and write `results/ACCEPTANCE_FAILURES.md` with the diff before continuing. Do not tune the model to hit numbers.
6. All randomness uses `numpy.random.default_rng(seed)`; seeds live in `config.yaml`. `GroupKFold` is deterministic.
7. Offline at runtime. No network calls, no telemetry, **no logging of request bodies**.
8. A research-use banner appears on every page, every report and every API response (`banner` field).
9. Commit small. Keep `CHANGELOG.md`. After each milestone run `make verify` and paste the output into `docs/PROGRESS.md`.
10. Do not commit raw data. `data/raw/` and `data/processed/` are in `.gitignore`.
11. Python 3.12, type hints, `ruff`-clean, tests with `pytest`.
12. When something in this file is ambiguous, pick the simplest option, write the decision into `docs/DECISIONS.md`, and continue.

## 2. Repository layout (create exactly)

```
bdq-triage/
  AGENTS.md                 (this file)
  README.md  CHANGELOG.md  Makefile  config.yaml  requirements.txt  Dockerfile  .gitignore
  reference/                (given scripts; do not edit; used as golden implementation)
  data/raw/                 (downloaded; ignored)
  data/processed/           (parquet; ignored)
  src/bdq_triage/
    __init__.py  config.py  grammar.py  prep.py  splits.py
    catalogue.py  features.py  model.py  prevalence.py  decision.py
    evidence.py  structure.py  report.py  cli.py
    api/ app.py  schemas.py
    ui/ templates/*.html  static/(css, js; no external CDNs)
  eval/
    run_all.py  e01_prep.py ... e13_ppv_npv.py  metrics.py  tables.py
  results/                  (JSON + CSV + MD written by eval; committed)
  models/                   (joblib + metadata.json; committed)
  tests/
  docs/  PROGRESS.md  DECISIONS.md  MODEL_CARD.md  RESULTS.md  DEMO_SCRIPT.md  GLOSSARY.md
  demo/  held_out_samples.csv  example_batch.csv
```

## 3. Environment

`requirements.txt` (pin exactly):
```
pandas>=2.2
numpy>=1.26
scipy>=1.11
scikit-learn>=1.4
pyarrow>=15
joblib>=1.3
matplotlib>=3.8
catomatic==0.2.3
piezo==0.9.2
fastapi>=0.110
uvicorn>=0.29
pydantic>=2.6
jinja2>=3.1
pyyaml>=6
pytest>=8
httpx>=0.27
ruff>=0.4
```
Create a venv (`python3.12 -m venv .venv`). `matplotlib-venn` is needed **only** to import the upstream repo's own `protocols/` code. Do not import upstream code in production; reimplement.

## 4. Data acquisition

Source: public repository `fowler-lab/tb-bdq-cat` (MIT licence). It holds processed genotype and phenotype tables, catalogues, wildcard rules and a structure file.

```bash
mkdir -p data/raw
curl -L -o data/raw/tb-bdq-cat.tar.gz https://codeload.github.com/fowler-lab/tb-bdq-cat/tar.gz/refs/heads/main
tar -xzf data/raw/tb-bdq-cat.tar.gz -C data/raw        # creates data/raw/tb-bdq-cat-main/
```
**Do not** use the GitHub REST API (rate-limited) or a shallow `git clone` (it failed in testing). Use the tarball.

After download, write `data/raw/MANIFEST.json` with the SHA-256 and byte size of every file under `tb-bdq-cat-main/data`, `catalogues`, `wildcards`. The model metadata must record this manifest hash.

Files used (under `data/raw/tb-bdq-cat-main/`):

| File | Rows (observed) | Required columns (assert they exist) |
|---|---|---|
| `data/GENOMES.csv.gz` | 44,139 | `UNIQUEID, SITEID, SUBJID` |
| `data/DST_MEASUREMENTS_BDQ_CRYPTIC.csv.gz` | 21,261 | `UNIQUEID, DRUG, SOURCE, METHOD_3, METHOD_MIC, PHENOTYPE, QUALITY` |
| `data/DST_MEASUREMENTS_BDQ_NICD.csv.gz` | 1,821 (1,169 distinct ids) | same columns; all `METHOD_3 == MGIT960` |
| `data/MUTATIONS_BDQ.csv.gz` | large (all genes) | `UNIQUEID, GENE, MUTATION, IS_MINOR_ALLELE, MINOR_MUTATION, FRS, CODES_PROTEIN, AMINO_ACID_NUMBER` |
| `data/Rv0678_WT.pdb` | n/a | ATOM records; chains A,B protein (residues 15 to 163) and C,D DNA |
| `catalogues/WHOv2/NC_000962.3_WHO-UCN-TB-2023.5_v2.0_GARC1_RFUS.csv` | 1 file | WHO v2 catalogue (GARC format) |
| `wildcards/bdq_wildcards_no_mmpL5.json` | 1 file | rule file for the catalogue builder |

## 5. Data preparation (`prep.py`, golden copy: `reference/prep.py`)

Steps, in order:

1. Read the two DST tables. CRyPTIC: **keep only `QUALITY == "HIGH"`**. NICD: keep all (all one quality grade). Add `SRC` = `"CRyPTIC"` or `"NICD"`.
2. Stack them. If a `UNIQUEID` has several rows, keep one row: prefer `PHENOTYPE == "R"`, then rows with a non-null `METHOD_MIC`. Implement with `sort_values` plus `drop_duplicates(keep="first")`. **Do not use `groupby().apply()`**: in current pandas it drops the grouping column.
3. Inner-join with `GENOMES` on `UNIQUEID` (drop duplicate genome rows first) to add `SITEID`, `SUBJID`.
4. `y = (PHENOTYPE == "R")`.
5. Mutations: restrict to `UNIQUEID` in the sample table. Build the name `MUT = GENE + "@" + (MINOR_MUTATION if IS_MINOR_ALLELE else MUTATION)`. Set `FRSf = FRS if IS_MINOR_ALLELE else 1.0`.
6. Let `suffix = MUT.split("@")[-1]`. Flags:
   - `IS_SYN` = `suffix[0] == suffix[-1]` (synonymous: same first and last character, e.g. `S68S`).
   - `IS_NULL` = `suffix[-1] == "Z"` (null call, position not callable).
7. Parsed flags for non-synonymous, non-null rows: `is_stop` (`suffix` ends with `!`), `is_indel` (`suffix` contains `_ins_`, `_del_`, `indel` or `del_`), `is_prom` (matches `^[acgt]-\d+`), `is_mis` (matches `^[A-Z]\d+[A-Z]$`), `pos` (digits of a missense).
8. Save `data/processed/samples.parquet` and `data/processed/mutations.parquet`.

**Mutation grammar (GARC-style):**

| Example | Meaning |
|---|---|
| `R50Q` | missense, residue 50, R to Q |
| `141_ins_c` | insertion of c at nucleotide 141 (frameshift) |
| `193_del_g` | deletion at nucleotide 193 |
| `135_minorindel` | minor-allele indel (appears in `MINOR_MUTATION`) |
| `Q131!` | stop codon at residue 131 |
| `c-11a` | promoter change 11 bases upstream |
| `S68S` | synonymous, ignore |
| `S68Z` | null call, ignore |

**Asserted counts (exact):** 11,868 samples; 866 resistant; CRyPTIC 10,704 samples with 83 resistant; NICD 1,164 samples with 783 resistant; 11,611 unique `SUBJID`; 25,554 mutation rows for these samples (before dropping synonymous/null); 1,115 minor-allele rows (median `FRSf` about 0.38); 11 null rows in 4 samples.

**Gene roles:** `CORE = [Rv0678, atpE, pepQ]` provide variant features. `mmpL5` and `mmpS5` are **excluded from the catalogue** and used only as an `mmpL5` loss-of-function flag (indel or stop in `mmpL5`).

## 6. Splits (`splits.py`)

| ID | Definition |
|---|---|
| S1 | `GroupKFold(n_splits=5)`, groups = `SUBJID`, on samples in table order |
| S2 | leave-one-site-out: hold out each distinct `SITEID` (23 values). NICD samples have `SITEID == 10`, the same site as CRyPTIC South Africa |
| S3 | cross-source: train on one `SRC`, test on the other, both directions |
| S4 | variant-set grouping for the novel-variant test: group key = sorted `"|".join(MUTATION strings)` |

## 7. Tier A: catalogue backbone (`catalogue.py`, golden copy: `reference/catbase.py`)

```python
from catomatic.BinaryCatalogue import BinaryBuilder
import piezo, json, tempfile

samples_cols = ["UNIQUEID","DRUG","SOURCE","METHOD_3","PHENOTYPE","METHOD_MIC"]
muts = mutations[(~IS_SYN) & (~IS_NULL) & (~GENE.isin(["mmpL5","mmpS5"]))]   # columns renamed: MUTATION=MUT, FRS=FRSf
cat = BinaryBuilder(samples=train_samples[samples_cols], mutations=muts_train, frs=0.1).build(
        test="Binomial", background=0.1, tails="two", strict_unlock=True)
wild = json.load(open("wildcards/bdq_wildcards_no_mmpL5.json"))
df = cat.build_piezo("NC_000962.3", "BDQ-cv", "1.0", "BDQ", wild, json_dumps=True)   # json_dumps=True is REQUIRED
df.to_csv(path, index=False); rc = piezo.ResistanceCatalogue(path)
```
Library facts (catomatic 0.2.3, piezo 0.9.2): the constructor argument is lowercase **`frs`**; `build_piezo(..., json_dumps=True)` is needed or `piezo` cannot read the CSV; `rc.predict(variant)` returns a dict like `{"BDQ": "R"}` (also `"S"` or `"U"`). A full catalogue build takes about 17 seconds, so cache per fold.

**Applying a catalogue to a sample:** take that sample's mutations with `FRSf >= 0.1` (same filters as above), call `rc.predict("GENE@MUT")` for each, and aggregate: `R` if any variant is R; else `U` if any is U; else `S` (including samples with no variants).

**WHO v2 as-is (reference baseline):** load the WHO file with `piezo.ResistanceCatalogue`. For each sample, use non-synonymous, non-null variants with `FRSf >= thr`. Wrap `predict` in `try/except` and treat an exception as `S`. Run for `thr = 0.75` and `thr = 0.1`.

## 8. Features for Tier B (`features.py`, golden copy: `reference/ml.py`)

Per sample, using non-synonymous, non-null mutation rows with `FRSf >= frs_min` (default 0.1; ablation uses 0.75):

1. **Variant table**: for rows in `CORE` genes, one value per `(UNIQUEID, MUT)` = max `FRSf`.
2. **Aggregates** (each is the **max `FRSf`** over rows in the class, 0 if none):
   - `rv_missense`, `rv_mis_dbd` (missense with `34 <= pos <= 99`), `rv_indel`, `rv_stop`, `rv_prom` (all Rv0678)
   - `atpE_any`, `pepQ_any`
   - `mmpL5_LoF` (mmpL5 indel or stop)
   - `rv_any = max(rv_missense, rv_indel, rv_stop, rv_prom)`
   - `rv_x_L5LoF = rv_any * mmpL5_LoF`
3. **Vocabulary**: variants that appear in at least **3 training samples** (count over the training ids only). One column per variant, value = max `FRSf`.
4. **Feature sets:** `F0` = vocabulary only. `F1` = F0 + `rv_missense, rv_mis_dbd, rv_indel, rv_stop, rv_prom, atpE_any, pepQ_any`. `F2` = F1 + `mmpL5_LoF, rv_x_L5LoF`. **Production = F2 with `frs_min = 0.1`.**

## 9. Model, calibration, thresholds (`model.py`)

- Model: `LogisticRegression(C=0.3, penalty="l2", max_iter=4000)` with **no** class weights.
- **Nested calibration** inside every training set: inner `GroupKFold(min(5, n_groups))` over the training rows with groups = `SUBJID` gives out-of-fold raw scores. Fit Platt scaling: `LogisticRegression(C=1e6)` on `logit(raw)` (clip raw to `[1e-6, 1-1e-6]`). The final model is fit on all training rows; test scores go raw, then Platt.
- **Thresholds** from the calibrated inner out-of-fold scores `poof` and labels `ytr` (use unique values of `round(poof, 4)` as candidates):
  - `t_bin`: threshold maximising Youden `sens + spec - 1`.
  - `lo`: the **largest** `t` such that at least 20 rows have `poof <= t` and their NPV (`1 - mean(ytr)`) is `>= 0.985`.
  - `hi`: the **smallest** `t` such that at least 20 rows have `poof >= t` and their PPV (`mean(ytr)`) is `>= 0.80`.
  - enforce `hi >= lo + 1e-6`.
- Ternary rule: `S` if `p <= lo`, `R` if `p >= hi`, else `U`.
- Production artefact: train on **all** 11,868 samples with the same inner procedure; save `models/bdq_lr.joblib` and `models/metadata.json` (versions, manifest hash, training prevalence `866/11868`, `lo`, `hi`, `t_bin`, feature names, metrics file hashes).

## 10. Prevalence correction (`prevalence.py`)

Probabilities are valid only at the training prevalence `pi_tr`. For a lab prevalence `pi`:
```
odds(p)   = p / (1 - p)
p_adj     = odds' / (1 + odds'),   odds' = odds(p) / odds(pi_tr) * odds(pi)
```
For threshold selection at a target prevalence, **reweight** the inner out-of-fold rows: weight `pi/pi_tr` for resistant rows and `(1-pi)/(1-pi_tr)` for susceptible rows, compute weighted NPV/PPV, and apply the same `0.985` / `0.80` rules on `p_adj`. Show PPV and NPV at the user's prevalence using cross-validated sensitivity and specificity from `results/`. Default prevalence in the UI is **5%, labelled as a placeholder the lab must replace**.

## 11. Experiments and acceptance numbers (`eval/`)

Write each result as JSON (+ CSV or MD table) in `results/`. Tolerances are absolute. If outside tolerance, follow rule 5.

| ID | What | Golden script | Expected (tolerance) |
|---|---|---|---|
| E01 | Prep counts | `prep.py` | exact numbers in section 5 |
| E02 | WHO v2 as-is, pooled | `catbase.py` | `thr 0.75`: firm 0.968, sens among firm 0.561, spec 0.987, **U counted as S: sens 0.497 (+-0.005)**. `thr 0.1`: firm 0.961, sens 0.801, spec 0.981, **U as S: sens 0.679 (+-0.005)** |
| E03 | Catalogue CV on S1 (build per fold, predict test) | `stage2.py` | firm 0.944 (+-0.01); sens among firm 0.792 (+-0.015); spec 0.984 (+-0.004); U as S: sens 0.604, spec 0.984; U as R: sens 0.842, spec 0.943; **661 U samples, 206 resistant** (+-15). In-sample (built on all data, optimistic, label it so): firm 0.951, sens 0.816, spec 0.984 |
| E04 | Tier B main (F2) on S1, nested | `run_ml.py` | AUC 0.909 (+-0.01); AUPRC 0.684 (+-0.02); binary sens 0.835 (+-0.02), spec 0.948 (+-0.01); ternary firm 0.947 (+-0.01), sens 0.735, spec 0.989; sensitivity at specificity 0.95/0.97/0.98/0.985 = 0.830/0.782/0.652/0.604 (+-0.03). Patient-bootstrap 500 reps for 95% CIs |
| E05 | Ablation on S1 | `run_ml.py` | AUC / AUPRC: F0 0.791/0.503; F1 0.905/0.672; F2 0.909/0.684; F2 with `frs_min=0.75`: 0.799/0.548 (+-0.015) |
| E06 | Cross-source (S3), nested inside training source | `run_ml.py` | CRyPTIC to NICD: AUC 0.91, sens 0.844, spec 0.963. NICD to CRyPTIC: AUC 0.83, sens 0.675, spec 0.953 (+-0.03) |
| E07 | Leave-one-site-out (S2) | `run_ml.py` | pooled binary sens 0.788 (+-0.02), spec 0.954 (+-0.01). Per-site specificity for sites with n >= 300 lies in 0.944 to 0.979. **Do not report a pooled LOSO AUC** (mixes 23 calibrations). Report per-site sensitivity only when the site has >= 25 resistant samples |
| E08 | Novel-variant test (S4) | section 11.1 | AUC about 0.57 (accept 0.50 to 0.65). Report as a negative result. Do not tune |
| E09 | Hybrid on catalogue-U samples (fold `lo/hi`) | `stage4.py` | ML AUC inside U = 0.744 (+-0.03). Resolves 173 as S (17 resistant) and 126 as R (51 susceptible); 362 remain (114 resistant) (+-15 each). Report as exploratory; **never auto-call** |
| E10 | Prevalence-shift test | `stage4.py` | NICD to CRyPTIC: ECE raw 0.288 to 0.013 (+-0.01), Brier 0.102 to 0.017 (+-0.005). CRyPTIC to NICD: ECE 0.568 to 0.101 (+-0.03), Brier 0.501 to 0.124 (+-0.02). ECE uses 10 quantile bins |
| E11 | Imbalance + learning curve | `stage5.py` | trained 1:1 balanced: AUC 0.904, AUPRC 0.653, sens@spec98 0.592. Fraction of resistant training rows 0.10/0.25/0.50/0.75/1.00: AUC 0.898/0.905/0.908/0.906/0.909 (+-0.01) |
| E12 | Variant priority list + evidence cards | `stage3.py` | 332 distinct variant sets among the 661 U samples; top-5 sets cover 28%, top-10 34%, top-20 40%, top-50 50% of U samples (+-3 points). Largest: `Rv0678@c-11a` (102 samples, 9% resistant), then `pepQ@G197R` (32, 0 resistant) |
| E13 | PPV/NPV table | `stage3.py` | model (sens .83, spec .95) at prevalence 1/2/5/10%: PPV 14/25/47/65%, NPV 99.8/99.6/99.1/98.1%. Catalogue (.79, .984): PPV 33/50/72/85%, NPV 99.8/99.6/98.9/97.7% |

Also record: model predict latency (about 1 ms per sample, 135 features) as `results/latency.json`.

### 11.1 Novel-variant test (E08) and structure features (`structure.py`)

Parse `Rv0678_WT.pdb` (plain text, no extra libraries): chain A residues 15 to 163 with atoms; chains C and D are DNA. For each chain-A residue compute (a) the minimum distance from any atom to any DNA atom, (b) burial = number of CA atoms (chains A and B) within 10 angstrom of the residue's CA. Assumption: the file's residue numbers equal the gene's amino-acid numbers; write this into `docs/DECISIONS.md` and the model card.
Subset: samples whose Rv0678 variants are all missense and that carry no atpE/pepQ variant (expect about 328 samples, 170 distinct variant sets). Features per sample: min distance over its missense positions (default 30.0 if none in structure), max burial (default 0), DBD flag, count of positions, outside-structure flag. Model: `StandardScaler` + `LogisticRegression(C=0.5, class_weight="balanced")`. CV: `GroupKFold(5)` grouped by variant-set key (S4). Compare "region flag + count" against "+ structure". Both are expected near AUC 0.56 to 0.57. **Report it as the reason the tool abstains on new variants.**

## 12. Product (`src/bdq_triage`)

### 12.1 Decision policy (fixed)
- **Final call = Tier A catalogue aggregate** (R / S / U). The model **never changes the call**.
- The model supplies: `probability_at_prevalence`, a **rank** for unresolved samples, and disagreement warnings.
- Production catalogue: built on all 11,868 samples with the section 7 settings; stored in `models/catalogue_bdq.csv` (piezo format) plus `models/catalogue_version.json`.

### 12.2 Warnings (codes)
| Code | Trigger |
|---|---|
| `MINOR_ALLELE_ONLY` | all contributing variants have `FRS < 0.9` |
| `NOVEL_VARIANT` | variant not in the evidence store |
| `LOW_EVIDENCE` | variant carriers < 10 |
| `NULL_CALL` | any input variant ends in `Z` |
| `PUMP_DAMAGE_MODIFIER` | `mmpL5_LoF > 0` together with an Rv0678 variant |
| `NO_VARIANT_INVISIBLE_RESISTANCE_POSSIBLE` | call S with no variants (about 1 in 6 resistant germs carry no known variant) |
| `MODEL_DISAGREES` | catalogue S but model `p >= hi`, or catalogue R but `p <= lo` |

### 12.3 API (FastAPI; `uvicorn bdq_triage.api.app:app`)
- `POST /predict` body: `{"sample_id": "S1" (optional, regex ^[A-Za-z0-9_-]{1,40}$), "prevalence": 0.05 (optional, 0.001 to 0.5), "variants": [{"gene": "Rv0678", "mutation": "141_ins_c", "frs": 1.0}]}`. Validate `gene` in `{Rv0678, atpE, pepQ, mmpL5, mmpS5}`; validate mutation against the grammar in section 5; reject free text.
- Response fields: `call` (`R|S|U`), `call_label`, `tier_a` (`variant_calls[]`), `model` (`score_raw`, `probability_train_prevalence`, `probability_at_prevalence`, `queue_rank_percentile`), `predictive_values` (`prevalence`, `ppv_if_R`, `npv_if_S`), `evidence[]`, `warnings[]`, `versions` (package, catalogue, model, data manifest hash), `banner`.
- `evidence[]` item: `variant`, `n_carriers`, `n_resistant`, `pct_resistant`, `ci95_low`, `ci95_high` (Wilson), `by_source` (`CRyPTIC`, `NICD` with `n`, `r`), `note` ("NICD is enriched for resistance; read per-source numbers").
- `GET /evidence/{gene}/{mutation}`; `POST /batch` (long-format CSV: `sample_id,gene,mutation,frs`) returns a CSV with the call and, for U samples, a rank by `probability_at_prevalence` (descending); `GET /model-card` returns `docs/MODEL_CARD.md` rendered; `GET /health`.
- PPV and NPV use cross-validated sensitivity and specificity read from `results/` (catalogue values for the shown call), via Bayes' rule at the user's prevalence.

### 12.4 CLI
`bdq-triage predict --variants file.csv --prevalence 0.05`, `bdq-triage batch in.csv out.csv`, `bdq-triage eval --all`, `bdq-triage verify` (runs acceptance checks).

### 12.5 UI (server-rendered Jinja; no external CDN, no external fonts)
Pages: `/` input (paste or CSV upload, prevalence slider), `/result` (big call, evidence cards, warnings, PPV/NPV, "what this means" text), `/queue` (batch table sorted by priority), `/priority` (variant priority list from E12), `/performance` (every E-table from `results/` with intervals), `/model-card`, `/about`.
**Replay mode** on `/`: a dropdown of held-out germs from `demo/held_out_samples.csv`. These must show the **out-of-fold** catalogue call and model probability stored in `results/oof_predictions.parquet` (not the production model, which has seen them). Label it "held-out in cross-validation". Live prediction is labelled "model trained on all public data".
Structure viewer: P2. Minimum: a table of Rv0678 positions with distance to DNA and DNA-binding-region flag. Only add a 3D viewer if the library is bundled locally.

### 12.6 Priorities
P0: sections 5 to 11, `/predict`, CLI, result page, performance page, tests, README, model card. P1: batch, queue, priority list, replay mode, Docker. P2: structure viewer, drug-plugin skeleton (`DrugPlugin` protocol: genes, tables, catalogue, features, model), PDF export.

## 13. Tests (`tests/`, must pass in CI)

1. Grammar: every row of the section 5 table parses to the expected flags; `Z` is null; first==last is synonymous.
2. Prep: the exact counts in E01.
3. Golden: WHO v2 gives 0.497 and 0.679 within +-0.005.
4. No leakage: for each fold, the vocabulary and the calibrator use only training ids (assert id sets are disjoint).
5. Patient grouping: no `SUBJID` appears in both train and test of any S1 fold.
6. Calibration: on the pooled out-of-fold scores, binned reliability error below 0.03.
7. Prevalence formula: `p_adj` at `pi == pi_tr` equals `p`.
8. API: `/predict` for `Rv0678@141_ins_c` returns call `R`; for `pepQ@G197R` returns a non-`R` call; invalid gene returns HTTP 422; request bodies are not logged.
9. Determinism: running E04 twice gives identical JSON.
10. Banner present in every HTML page and API response.

## 14. Documents to generate (from `results/`, never hand-typed numbers)

- `README.md`: what it is, how to run in 5 commands, what it does not do.
- `docs/MODEL_CARD.md`: intended use, data (sources, counts, known issues), method, metrics with intervals, limits (no independent validation, about 16% invisible resistance, lab dependence, India sensitivity untestable, new variants unresolved), ethical notes, versions.
- `docs/RESULTS.md`: all E-tables.
- `docs/DEMO_SCRIPT.md`: the 5-minute demo (clear R case `141_ins_c`; unresolved case `c-11a` shown as Send-for-lab-test; prevalence slider 7% to 2% showing PPV falling; performance page; reproduction table).
- `docs/GLOSSARY.md`: plain-language definitions.

## 15. Milestones and gates

| Milestone | Work | Gate |
|---|---|---|
| M0 | repo, venv, config, Makefile, tarball download, manifest | tarball verified |
| M1 | `grammar.py`, `prep.py`, E01 | **exact counts** |
| M2 | Tier A catalogue, E02, E03 | **G1: E02 and E03 within tolerance** |
| M3 | features, model, nested calibration, E04, E05 | **G2: E04 within tolerance** |
| M4 | E06 to E13, `results/` complete, `oof_predictions.parquet` | all section 11 rows pass |
| M5 | `decision.py`, `evidence.py`, production artefacts in `models/` | `/predict` unit tests pass |
| M6 | API + CLI | tests 8 to 10 pass |
| M7 | UI pages (P0 then P1) | click-through works offline |
| M8 | docs + demo CSV + Dockerfile | `docker run` serves the app |
| M9 | hardening, `make verify` green | all tests green |

**Fallback rule:** if G1 or G2 fails and cannot be fixed within a day, ship Tier A only (catalogue tool) with the honest evaluation, and say so in the model card.

`Makefile` targets: `setup`, `data`, `prep`, `eval`, `verify`, `test`, `app`, `demo`, `docs`, `docker`.

## 16. Pitfalls already hit (do not repeat)

1. GitHub API rate limit and failed shallow clone: use the tarball.
2. `catomatic` argument is `frs`, not `FRS`; `build_piezo` needs `json_dumps=True`.
3. `groupby().apply()` drops the grouping column in current pandas.
4. Importing the upstream `protocols/` code needs `matplotlib-venn`; do not import it in production.
5. `Z` suffix means null call, not synonymous; synonymous means first and last character are equal.
6. Minor alleles are named in `MINOR_MUTATION`; use it when `IS_MINOR_ALLELE`.
7. `class_weight="balanced"` distorts probabilities; calibrated model must not use it.
8. Probabilities do not transfer across prevalence; always correct (section 10).
9. Do not report pooled leave-one-site-out AUC.
10. NICD ids carry site 10, shared with CRyPTIC South Africa (817 resistant in total at that site).
11. Numeric `SITEID` has no public country key. Site 5 (4,030 genomes) is **inferred** to be India; mark it "unverified" everywhere and never present Indian sensitivity (only about 5 resistant).
12. Catalogue build takes about 17 s; cache per fold in `data/processed/cache/`.

## 17. Claims policy (UI, README, model card)

**Allowed statements:** "reproduces the published catalogue performance"; "about 79 to 84% sensitivity and 95 to 98% specificity by cross-validation"; "cross-lab sensitivity ranged from about 63% to 84%"; "about 1 in 6 resistant germs carry no known variant"; "R calls are flags for confirmation at realistic prevalence"; "research use only".
**Forbidden statements:** "more accurate than WHO/Oxford"; "validated in hospitals"; "works on Indian patients" (beyond false-alarm audit); "saves money" or any rupee figure; "diagnoses" or "recommends treatment"; "predicts new variants".

## 18. Definition of done

- [ ] `make verify` passes: E01 to E13 within tolerance, all tests green.
- [ ] `docker run` serves the app offline.
- [ ] Every number on `/performance` and in `docs/` traces to a file in `results/`.
- [ ] Banner and model card present; forbidden statements absent (grep test).
- [ ] `docs/PROGRESS.md` shows each gate result.
- [ ] Demo runs from `docs/DEMO_SCRIPT.md` with no internet.