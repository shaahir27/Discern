"""Model-card content rendered from the registered evaluation reports."""

from .research import REPORTS, Reports


def model_card(reports: Reports) -> list[tuple[str, list[str]]]:
    prep, cat, model, mechanism, km, conformal = (
        reports.observed(name) for name in ("prep", "catalogue", "model", "mechanism", "km", "conformal")
    )
    context = reports.get(REPORTS["intervals"])
    coverage = context["observed"]["interval_level"]
    invisible = context["observed"]["fraction_without_retained_core_variant"]
    intervals = [
        f"{row['metric']} {row['estimate']:.4f} ({coverage:.0%} interval "
        f"{row['lower']:.4f}–{row['upper']:.4f})."
        for row in context["extra"]["intervals"]
    ]
    return [
        ("Intended use", [
            ("Discern presents archived research on bedaquiline resistance in tuberculosis. "
            "It supports reviewing published evidence, patient-grouped evaluation, and example records. "
            "Clinical decisions require independent laboratory confirmation."),
        ]),
        ("Training data", [
            (f"The reference cohort has {prep['n_samples']:,} isolates, including {prep['n_resistant']:,} resistant. "
            f"CRyPTIC contributes {prep['cryptic_n']:,} high-quality isolates "
            f"({prep['cryptic_R']:,} resistant). NICD contributes {prep['nicd_n']:,} isolates "
            f"({prep['nicd_R']:,} resistant) and is enriched for resistance."),
            "Cohort composition and susceptibility-testing methods affect how estimates transfer.",
        ]),
        ("Evaluation", [
            (f"Patient-grouped cross-validation gave the golden catalogue {cat['sens_firm']:.2%} sensitivity "
            f"and {cat['spec_firm']:.2%} specificity among firm calls (E03). "
            f"The calibrated F2 logistic regression had ROC AUC {model['auc']:.4f} "
            f"and average precision {model['auprc']:.4f} (E04)."),
            (f"The audited catalogue downgrade changed {mechanism['s1_n_flips']} S1 calls to lab referral, "
            f"including {mechanism['s1_n_flips_resistant']} phenotypically resistant isolates (E16)."),
            (f"KM includes {km['n_isolates']} non-overlapping isolates with {km['n_resistant']} resistant. "
            f"The archived model AUC was {km['model_auc']:.4f}. The cohort is a sensitivity analysis with "
            "medium-quality labels and undocumented assay and critical concentration."),
            ("The performance archive presents the uncertainty intervals that accompany the saved reports. "
            "Evidence cards give Wilson intervals for the observed carrier proportions."),
        ]),
        ("Limits and interpretation", [
            (f"In S1, {invisible:.1%} of resistant isolates lacked a retained core-gene variant. "
             "An archived susceptible call can therefore coexist with a resistant laboratory phenotype."),
            ("The available external cohort does not establish sensitivity in India. Results can vary by "
            "laboratory, lineage mix, and local prevalence."),
            (f"The conformal S-zone missed its target. Its best achievable bound was "
            f"{conformal['s_best_achievable_bound']:.4%} at delta {conformal['delta']:.2f}. "
            "Certificate bounds failed to transfer across labs in the archived experiments."),
        ]),
        ("Uncertainty intervals", intervals),
        ("Replay provenance", [
            ("Replay retrieves an immutable prediction from the fold that held out the displayed isolate. "
            "The recorded catalogue call and the model probability describe different parts of the study. "
            "The evidence cards summarize training-cohort carriers, with source counts shown separately."),
        ]),
    ]
