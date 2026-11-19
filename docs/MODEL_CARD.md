# Discern model card

Research use only. For evaluation and education. Clinical decisions require laboratory confirmation.

## Intended use

Discern presents archived research on bedaquiline resistance in tuberculosis. It supports reviewing published evidence, patient-grouped evaluation, and example records. Clinical decisions require independent laboratory confirmation.

## Training data

The reference cohort has 11,868 isolates, including 866 resistant. CRyPTIC contributes 10,704 high-quality isolates (83 resistant). NICD contributes 1,164 isolates (783 resistant) and is enriched for resistance.

Cohort composition and susceptibility-testing methods affect how estimates transfer.

## Evaluation

Patient-grouped cross-validation gave the golden catalogue 79.24% sensitivity and 98.35% specificity among firm calls (E03). The calibrated F2 logistic regression had ROC AUC 0.9085 and average precision 0.6836 (E04).

The audited catalogue downgrade changed 41 S1 calls to lab referral, including 2 phenotypically resistant isolates (E16).

KM includes 792 non-overlapping isolates with 35 resistant. The archived model AUC was 0.8682. The cohort is a sensitivity analysis with medium-quality labels and undocumented assay and critical concentration.

The performance archive presents the uncertainty intervals that accompany the saved reports. Evidence cards give Wilson intervals for the observed carrier proportions.

## Limits and interpretation

In S1, 15.8% of resistant isolates lacked a retained core-gene variant. An archived susceptible call can therefore coexist with a resistant laboratory phenotype.

The available external cohort does not establish sensitivity in India. Results can vary by laboratory, lineage mix, and local prevalence.

The conformal S-zone missed its target. Its best achievable bound was 1.5516% at delta 0.05. Certificate bounds failed to transfer across labs in the archived experiments.

## Uncertainty intervals

S1 catalogue sensitivity 0.7924 (95% interval 0.7596–0.8243).

S1 catalogue specificity 0.9835 (95% interval 0.9810–0.9858).

S1 model AUC 0.9085 (95% interval 0.8944–0.9230).

S1 model AUPRC 0.6836 (95% interval 0.6471–0.7194).

KM catalogue sensitivity 0.7917 (95% interval 0.5953–0.9076).

KM catalogue specificity 0.9881 (95% interval 0.9757–0.9942).

## Replay provenance

Replay retrieves an immutable prediction from the fold that held out the displayed isolate. The recorded catalogue call and the model probability describe different parts of the study. The evidence cards summarize training-cohort carriers, with source counts shown separately.
