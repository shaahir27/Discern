"""Unit tests for Sequencing-Assay & tNGS Panel Blind-Spot Auditor."""
from __future__ import annotations

from discern.decision import Predictor
from discern.panels import (
    PANELS,
    PanelProfile,
    audit_sample_panel,
    compute_cohort_panel_coverage,
)


def test_panel_profiles_exist():
    expected_keys = {"wgs", "deeplex_myc_tb", "deeplex_xl", "ampore_tb", "icmr_nirt_minion"}
    assert expected_keys.issubset(set(PANELS.keys()))
    for p in PANELS.values():
        assert isinstance(p, PanelProfile)
        assert p.vaf_floor > 0
        assert len(p.covered_genes) > 0


def test_panel_blind_spot_deeplex_and_ampore():
    # Sample with Rv0678 mutation
    variants = [{"gene": "Rv0678", "mutation": "M146T", "frs": 0.85}]

    # 1. Deeplex Myc-TB: covers only Rv0678
    audit_deeplex = audit_sample_panel(variants, "deeplex_myc_tb")
    status_d = audit_deeplex["gene_status"]
    assert status_d["Rv0678"] == "mutated"
    assert status_d["atpE"] == "not_interrogated"
    assert status_d["pepQ"] == "not_interrogated"
    assert status_d["mmpL5"] == "not_interrogated"
    assert status_d["mmpS5"] == "not_interrogated"
    assert audit_deeplex["epistatic_veto_status"] == "blind_missing_pump_genes"
    assert "EPISTATIC_VETO_BLIND" in audit_deeplex["warnings"]
    assert "UNSEQUENCED_CORE_LOCI" in audit_deeplex["warnings"]

    # 2. AmPORE-TB: covers Rv0678 and atpE (missing pepQ, mmpL5, mmpS5)
    audit_ampore = audit_sample_panel(variants, "ampore_tb")
    status_a = audit_ampore["gene_status"]
    assert status_a["Rv0678"] == "mutated"
    assert status_a["atpE"] == "wild_type"
    assert status_a["pepQ"] == "not_interrogated"
    assert status_a["mmpL5"] == "not_interrogated"
    assert status_a["mmpS5"] == "not_interrogated"
    assert audit_ampore["epistatic_veto_status"] == "blind_missing_pump_genes"
    assert "EPISTATIC_VETO_BLIND" in audit_ampore["warnings"]
    assert "UNSEQUENCED_CORE_LOCI" in audit_ampore["warnings"]

    # 3. WGS: all active, no veto blind warning
    audit_wgs = audit_sample_panel(variants, "wgs")
    assert audit_wgs["epistatic_veto_status"] == "active"
    assert "EPISTATIC_VETO_BLIND" not in audit_wgs["warnings"]
    assert "UNSEQUENCED_CORE_LOCI" not in audit_wgs["warnings"]


def test_homopolymer_guardrail():
    # Rv0678 192_ins_g is in 7-G homopolymer tract (nt 192-198)
    # Low FRS (0.18) on WGS (mixed platform) -> should warn HOMOPOLYMER_ARTIFACT_RISK
    audit_wgs_low = audit_sample_panel([{"gene": "Rv0678", "mutation": "192_ins_g", "frs": 0.18}], "wgs")
    assert "HOMOPOLYMER_ARTIFACT_RISK" in audit_wgs_low["warnings"]

    # High FRS (0.80) on Nanopore (icmr_nirt_minion) -> should warn HOMOPOLYMER_ARTIFACT_RISK because nanopore platform
    audit_nano_high = audit_sample_panel([{"gene": "Rv0678", "mutation": "192_ins_g", "frs": 0.80}], "icmr_nirt_minion")
    assert "HOMOPOLYMER_ARTIFACT_RISK" in audit_nano_high["warnings"]

    # High FRS (0.80) on Illumina WGS -> should NOT warn homopolymer risk
    audit_wgs_high = audit_sample_panel([{"gene": "Rv0678", "mutation": "192_ins_g", "frs": 0.80}], "deeplex_xl")
    assert "HOMOPOLYMER_ARTIFACT_RISK" not in audit_wgs_high["warnings"]

    # Non-homopolymer mutation (M146T) -> should NEVER warn homopolymer risk
    audit_clean = audit_sample_panel([{"gene": "Rv0678", "mutation": "M146T", "frs": 0.15}], "wgs")
    assert "HOMOPOLYMER_ARTIFACT_RISK" not in audit_clean["warnings"]


def test_below_lod_warning():
    # AmPORE has a 15% (0.15) LOD floor. A variant at 12% (0.12) should trigger BELOW_PANEL_LOD
    audit_lod = audit_sample_panel([{"gene": "Rv0678", "mutation": "M146T", "frs": 0.12}], "ampore_tb")
    assert "BELOW_PANEL_LOD" in audit_lod["warnings"]


def test_predictor_integration():
    pred = Predictor()
    res = pred.predict([{"gene": "Rv0678", "mutation": "141_ins_c", "frs": 0.95}], assay="deeplex_myc_tb")
    assert "panel_audit" in res
    assert res["panel_audit"]["assay_key"] == "deeplex_myc_tb"
    assert "EPISTATIC_VETO_BLIND" in res["warnings"]


def test_panel_audit_negative_controls():
    # Negative controls test for compute_cohort_panel_coverage
    ctrl_full = compute_cohort_panel_coverage(frozenset({"Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5"}))
    assert ctrl_full["known_locus_sensitivity"] == 1.0

    ctrl_empty = compute_cohort_panel_coverage(frozenset())
    assert ctrl_empty["known_locus_sensitivity"] == 0.0
