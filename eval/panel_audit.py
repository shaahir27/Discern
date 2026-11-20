"""Cohort-level panel coverage audit across commercial and research tNGS panels."""
from __future__ import annotations

import json
from pathlib import Path
import pandas as pd

from discern.panels import PANELS, compute_cohort_panel_coverage

RESULTS_DIR = Path("results/approaches/panel_audit/e24_panel_coverage")


def run_cohort_panel_audit() -> dict:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    results_dict = {
        "title": "Sequencing-Assay & tNGS Panel Blind-Spot Audit",
        "description": "Coverage upper bounds and blind spots across tNGS panels on the 866-isolate BDQ-R reference cohort",
        "reference_counts": {
            "total_resistant": 866,
            "known_locus_resistant": 729,
            "invisible_resistant_no_core_variant": 137,
        },
        "panels": {},
        "caveat": "Coverage is an upper bound on detection; tNGS depth, low VAF, and homopolymer calling errors can further reduce real-world sensitivity (Hall et al., eLife 2024; Ismail et al., AAC 2025).",
    }

    for key, p in PANELS.items():
        cov = compute_cohort_panel_coverage(p.covered_genes)
        row = {
            "panel_key": key,
            "panel_label": p.label,
            "platform": p.platform,
            "vaf_floor": p.vaf_floor,
            "covered_genes": ",".join(sorted(list(p.covered_genes))),
            "known_locus_sensitivity": cov["known_locus_sensitivity"],
            "overall_sensitivity_upper_bound": cov["overall_sensitivity_upper_bound"],
            "missed_known_pct": cov["missed_known_pct"],
            "missed_loci": ",".join(cov["missed_loci"]),
        }
        summary_rows.append(row)
        results_dict["panels"][key] = {**row, "notes": p.notes}

    # Negative controls
    ctrl_full = compute_cohort_panel_coverage(frozenset({"Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5"}))
    ctrl_empty = compute_cohort_panel_coverage(frozenset())
    results_dict["controls"] = {
        "full_coverage_control": ctrl_full,
        "empty_coverage_control": ctrl_empty,
    }

    df = pd.DataFrame(summary_rows)
    df.to_csv(RESULTS_DIR / "panel_summary.csv", index=False)

    with open(RESULTS_DIR / "results.json", "w", encoding="utf-8") as f:
        json.dump(results_dict, f, indent=2)

    print(f"Panel audit written to {RESULTS_DIR}")
    return results_dict


if __name__ == "__main__":
    run_cohort_panel_audit()
