# Discern demo

Research use only. Use the bundled recorded examples throughout the demo.

- Open **Check a sample**, download the example CSV, select it in the file picker, and choose **Open result**. Compare the catalogue call, lab finding, and saved model score. Open the sample details to find the study method.
- Choose a **Send for lab test** example. Open its carrier evidence and read the CRyPTIC and NICD counts separately. Explain the Wilson interval and the difference between the model score and the call.
- Move the prevalence slider to a lower value and update the estimates. The PPV moves with the assumed prevalence while the archived prediction stays fixed.
- Open **Review queue**, select the unresolved view, and inspect the order. Download the roster and import a subset of IDs to create a smaller view.
- Open **Study results**, select the catalogue and model experiments, then open the evaluation context report for the uncertainty intervals. Show the preregistered ship decision and the model card.

Reproduce the saved acceptance checks with `make verify`. Read [the result registry](RESULTS.md) for sources and [the decision log](DECISIONS.md) for the rationale.
