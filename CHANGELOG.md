# Changelog

## Unreleased

- Fixed the live catalogue review rule to require a high-confidence annotation and preserve independent resistant findings.
- Reworked sample intake, case results, review queues, and study views around their primary actions. Added dated visit controls and a one-click illustrative timeline.
- Fixed panel-report loading, batch prevalence rendering and manual overrides, native-export validation, and baseline selection for repeated patients.
- Scoped live results and timelines to the uploading browser and separated live links from saved replay.
- Updated the container to include the local artefacts and dependencies required at startup.
- Added the codebase guide and moved both novelty documents into `docs/`.

- Added offline acceptance verification for indexed result files and wired it into `make verify`.
- Added E21, a preregistered nested comparison of model classes on patient-grouped S1 and the held-out KM cohort.
- Recorded the E19 ship decision and enabled the audited Tier A veto alongside golden L2 logistic regression.
- Added the offline research explorer, a CSV-first recorded-sample workflow, and empty-input coverage for the existing predictor.
- Added located CSV errors and score provenance, deferred legacy catalogue failures, and corrected prevalence-based predictive values.
- Added Live CSV Quality Check to calculate fresh completeness scores from required fields, detect missing values, duplicates, and invalid formatting under "Upload quality", kept strictly separate from resistance results.
- Added live model scoring on uploaded CSV files, enabling dynamic evaluation of novel isolates with the production Tier A + Tier B model alongside archived replays.
- Added Epistatic Veto Layer to decision engine, suppressing false-resistant calls when efflux pump gene mmpL5 has an inactivating mutation, cutting catalogue false alarms by 22.2% across the cohort with zero loss of sensitivity.
- Added 3D biophysical structural modeling in structure.py, incorporating Grantham physicochemical distance, recognition helix alpha4 targeting (codons 62-68), and secondary structure disruptors, boosting novel-variant ROC-AUC on S4 from 0.549 to 0.688 (+13.9 points).
- Added subclonal heteroresistance detection (HETERORESISTANCE_SUBCLONAL) for minority alleles (0.10 <= FRS < 0.75), explaining assay discordance between 14-day microtiter plates (UKMYC) and liquid culture (MGIT960).
- Added dedicated Biological Mechanism & Clinical Guidance panel to /result template and 5 unit tests in tests/test_novelty_biophysics_and_epistasis.py.
- Added Sequencing-Assay & tNGS Panel Blind-Spot Auditor (`src/discern/panels.py`, `eval/panel_audit.py`, `GET /panels`, 5-gene coverage strip in UI) supporting WGS, Deeplex Myc-TB, Deeplex Myc-TB XL, AmPORE-TB, and ICMR-NIRT MinION.
- Added Unsupervised Batch Prevalence Estimation using Maximum Likelihood Label Shift (Saerens EM / MLLS) with bootstrap 95% CIs (`src/discern/label_shift.py`, dynamic slider filling, inline SVG PPV curve with uncertainty band in `/queue`).
- Added Longitudinal Serial-Isolate Tracking & Clonal-Ribbon Timeline (`src/discern/longitudinal.py`, optional `patient_id,collection_date` metadata, inline SVG Muller plot with timeline scrubber at `GET /patient/{patient_id}`, `demo/longitudinal_example.csv`).
- Added Native Lab-Export Importers & Mechanism-Stratified Reflex Guidance (`src/discern/importers.py`, HGVS to GARC variant normalization for TB-Profiler JSON, Deeplex Myc-TB CSV, and AmPORE-TB TSV/CSV; reason codes, domain tracks, and clinical reflex cards).
- Added Lightweight Batch-Recurrence & Contamination Alert (`src/discern/surveillance.py`, patient-deduplicated recurrence tracking, homopolymer artifact detection at poly-G hotspots vs. rare variant WGS surveillance escalation).
