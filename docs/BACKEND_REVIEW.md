# Backend review

Research use only.

## Input and output

The explorer accepts a CSV with recorded sample IDs or a complete recorded profile with the columns `sample_id,gene,mutation,frs`. The server checks the headers and cells, checks each ID against the bundled collection, and matches full profiles by their rows and fractions. Rows can be reordered. A successful match retrieves the saved catalogue call, laboratory finding, and model probability. A multi-sample upload opens the review queue. Changed or unknown profiles receive an actionable validation message. Uploaded content is never included in the error page.

The browser and CLI use the same import and matching code. Replay reads from the out-of-fold prediction file and checks identifiers, outcomes, fold labels, scores, cohort totals, and agreement with the demo fixture at startup. The catalogue result and model score stay separate. The prevalence control changes aggregate catalogue estimates while preserving the saved result.

## What the result screen needs from the backend

| User question | Data needed | Coverage today |
|---|---|---|
| Which sample did I open? | Stable sample ID, result origin, study source, and evaluation fold | Available for saved results. |
| What did the catalogue call? | Call, label, and structured reason codes tied to the evidence that affected the call | The saved call is available. The replay record lacks call-level reason codes and cannot say which observed variant drove the call. |
| What does the score mean? | A nullable score, its definition, reference population and prevalence, and an explanation when unavailable | Replay includes the saved probability, source experiment, and prevalence and counts for the training fold. Individual score uncertainty is marked as unavailable. |
| What did the lab find? | A separate lab outcome with provenance and an explicit unavailable state | The replay includes the recorded phenotype and study source. |
| Why did the upload fail? | Error category, optional CSV column and row, and a correction the user can make | The batch API returns a stable error code and the CSV line and column where available. The browser and CLI show the location and recovery message. Oversized inputs return HTTP 413. |

The screen presents the available values and labels the variant table as study
evidence. A variant listed there may have had no effect on the saved call. Supply
reason codes and evidence links before adding a section that claims to explain a call.

## Engine hardening

- Empty inputs used an empty, object-typed mutation table and caused column loss in pandas. Feature assembly now returns zero-filled aggregate rows for this case. The regression test checks the resulting research warning.
- A legacy catalogue exception or unsupported return value now defers the call for lab testing and adds a processing warning. A successful catalogue call still follows the frozen catalogue and veto policy.
- The legacy schema validates the complete mutation field, enforces matching gene prefixes, and rejects extra fields. The explorer continues to use exact archive matching.
- Scores depend on the study cohort and prevalence. KM provided a limited sensitivity analysis. The model card links the evaluation and its uncertainty.

## Predictive-value calculation

For an assumed prevalence across all samples, reweight separately the rates of
resistant calls and susceptible calls across the full cross-validation cohort.
The archive uses E03 call rates. The legacy engine uses the matching E03 or E16
policy rates. Unresolved calls contribute to the population rates and stay out of
the denominator for either firm call. An independent test calculates the same
estimates by reweighting the archived rows directly.

This corrects the estimates shown at a chosen prevalence. The archived calls,
model scores, acceptance metrics and training artifacts use the same frozen inputs.

[Engine](../src/discern/decision.py) · [Legacy validator](../src/discern/api/schemas.py) · [Archive data layer](../src/discern/research.py) · [Model card](MODEL_CARD.md)
