"""Unit tests for Longitudinal Serial-Isolate Tracking & Clonal-Ribbon Timeline."""
from __future__ import annotations

from fastapi.testclient import TestClient

from discern.api.app import app
from discern.longitudinal import analyze_patient_timelines, parse_timepoint_date


def test_parse_timepoint_date():
    # ISO date string
    d1, s1 = parse_timepoint_date("2025-01-15")
    assert s1 == "2025-01-15"
    assert d1 is not None

    # Day offset as integer or string
    d2, s2 = parse_timepoint_date("28")
    assert s2 == "Day 28"
    assert d2 == 28.0

    d3, s3 = parse_timepoint_date("Day 84")
    assert s3 == "Day 84"
    assert d3 == 84.0

    # Unparseable date returns None
    d_bad, _ = parse_timepoint_date("not-a-date")
    assert d_bad is None


def test_de_novo_and_rising_emergence():
    # Patient with 3 timepoints
    # t0: Day 0, WT (no core variants)
    # t1: Day 28, Rv0678@M146T frs=0.20
    # t2: Day 84, Rv0678@M146T frs=0.85 (rising from 0.20 -> 0.85, crossing 0.75)
    records = [
        {"sample_id": "S1", "patient_id": "PT_01", "collection_date": "Day 0", "variants": []},
        {"sample_id": "S2", "patient_id": "PT_01", "collection_date": "Day 28", "variants": [
            {"gene": "Rv0678", "mutation": "M146T", "frs": 0.20}
        ]},
        {"sample_id": "S3", "patient_id": "PT_01", "collection_date": "Day 84", "variants": [
            {"gene": "Rv0678", "mutation": "M146T", "frs": 0.85}
        ]},
    ]
    timelines = analyze_patient_timelines(records)
    assert "PT_01" in timelines
    pt = timelines["PT_01"]
    assert len(pt["timepoints"]) == 3
    alert_codes = [a["code"] for a in pt["alerts"]]
    assert "DE_NOVO_VARIANT_EMERGENCE" in alert_codes
    assert "RISING_SUBCLONE_TRAJECTORY" in alert_codes


def test_polyclonal_competition_and_replacement():
    # Patient with 2 competing Rv0678 variants:
    # Visit 1: Rv0678@M146T frs=0.45, Rv0678@141_ins_c frs=0.15
    # Visit 2: Rv0678@M146T frs=0.05 (drops by 0.40 >= 0.15), Rv0678@141_ins_c frs=0.80 (rises by 0.65 >= 0.15)
    records = [
        {"sample_id": "SMP1", "patient_id": "PT_POLY", "collection_date": "2025-01-01", "variants": [
            {"gene": "Rv0678", "mutation": "M146T", "frs": 0.45},
            {"gene": "Rv0678", "mutation": "141_ins_c", "frs": 0.15},
        ]},
        {"sample_id": "SMP2", "patient_id": "PT_POLY", "collection_date": "2025-03-01", "variants": [
            {"gene": "Rv0678", "mutation": "M146T", "frs": 0.05},
            {"gene": "Rv0678", "mutation": "141_ins_c", "frs": 0.80},
        ]},
    ]
    timelines = analyze_patient_timelines(records)
    assert "PT_POLY" in timelines
    pt = timelines["PT_POLY"]
    alert_codes = [a["code"] for a in pt["alerts"]]
    assert "POLYCLONAL_RV0678_COMPETITION" in alert_codes
    assert "CLONAL_REPLACEMENT" in alert_codes


def test_backward_compatible_4col_csv():
    # Records without patient_id or single timepoint
    records = [
        {"sample_id": "S1", "variants": [{"gene": "Rv0678", "mutation": "M146T", "frs": 0.95}]},
        {"sample_id": "S2", "patient_id": "SOLO", "collection_date": "Day 0", "variants": []},
    ]
    timelines = analyze_patient_timelines(records)
    assert timelines == {}


def test_patient_timeline_route_renders():
    # Upload multi-visit CSV through app and check /patient/{patient_id}
    csv_content = (
        "sample_id,gene,mutation,frs,patient_id,collection_date\n"
        "V1,Rv0678,M146T,0.20,PAT_DEMO_01,Day 0\n"
        "V2,Rv0678,M146T,0.85,PAT_DEMO_01,Day 28\n"
    )
    with TestClient(app) as client:
        resp = client.post("/review", data={"csv": csv_content, "assay": "wgs"}, follow_redirects=False)
        assert resp.status_code in (200, 303)
        # Check GET /patient/PAT_DEMO_01
        res_pt = client.get("/patient/PAT_DEMO_01")
        assert res_pt.status_code == 200
        assert "PAT_DEMO_01" in res_pt.text
        assert "timeline-data" in res_pt.text
