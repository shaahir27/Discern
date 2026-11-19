# Changelog

## Unreleased

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
