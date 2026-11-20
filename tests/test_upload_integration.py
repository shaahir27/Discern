import json
import re

import pytest
from fastapi.testclient import TestClient

from discern.api.app import app
from discern.importers import detect_and_parse_upload
from discern.label_shift import deduplicate_patient_probabilities
from discern.longitudinal import analyze_patient_timelines


@pytest.mark.parametrize("route", ["/review", "/queue"])
def test_batch_prevalence_renders_and_manual_default_can_override(route: str) -> None:
    payload = "sample_id\n" + "".join(f"TEST_{i}\n" for i in range(30))
    with TestClient(app) as client:
        queue = client.post(route, data={"csv": payload})
        assert queue.status_code == 200
        assert "Positive Predictive Value (PPV) vs Resistance Prevalence Curve" in queue.text
        link = re.search(r'href="(/result\?sample_id=TEST_0[^"]*)"', queue.text)
        assert link
        destination = link[1].replace("&amp;", "&")
        assert "Pre-filled from batch" in client.get(destination).text
        manual = client.get(destination + "&prevalence=0.05")
        assert manual.status_code == 200
        assert "Pre-filled from batch" not in manual.text


def test_panel_page_loads_saved_comparison() -> None:
    with TestClient(app) as client:
        result = client.get("/panels")
        assert result.status_code == 200
        assert "Deeplex" in result.text


def test_page_links_to_versioned_local_assets() -> None:
    with TestClient(app) as client:
        result = client.get("/")
        for asset in ("app.css", "app.js"):
            url = re.search(rf'(/static/{re.escape(asset)}\?v=\d+)"', result.text)
            assert url
            assert client.get(url[1]).status_code == 200


def test_live_batch_links_are_private_and_do_not_shadow_replay() -> None:
    with TestClient(app) as owner:
        archived_csv = owner.get("/example.csv").text
        archived_id = archived_csv.splitlines()[1].split(",")[0]
        payload = f"sample_id\n{archived_id}\nPRIVATE_SAMPLE\n"
        queue = owner.post("/queue", data={"csv": payload})
        assert queue.status_code == 200
        link = re.search(r'href="(/result\?sample_id=PRIVATE_SAMPLE[^"]*)"', queue.text)
        assert link
        destination = link[1].replace("&amp;", "&")
        assert owner.get(destination).status_code == 200
        assert "Saved model score" in owner.get("/result", params={"sample_id": archived_id}).text
        other = TestClient(app)
        assert other.get(destination).status_code == 404


@pytest.mark.parametrize("route", ["/review", "/queue"])
def test_invalid_optional_date_is_rejected_before_scoring(route: str) -> None:
    payload = (
        "sample_id,gene,mutation,frs,patient_id,collection_date\n"
        "NEW_SAMPLE,Rv0678,S68S,1,PATIENT_1,2026-02-30\n"
    )
    with TestClient(app) as client:
        result = client.post(route, data={"csv": payload})
        assert result.status_code == 422
        assert "collection_date" in result.text


def test_deeplex_header_selects_percent_units_even_below_one_percent() -> None:
    payload = b"Sample ID,Gene name,AA change,Freq (%)\nTEST,Rv0678,S68S,0.5\n"
    assert detect_and_parse_upload(payload)["rows"][0]["frs"] == 0.005


def test_ampore_csv_can_be_imported_without_a_filename() -> None:
    payload = b"sample,gene,variant,vaf\nTEST,Rv0678,S68S,0.5\n"
    assert detect_and_parse_upload(payload)["inferred_assay"] == "ampore_tb"


@pytest.mark.parametrize("variants", [[None], [{"gene": "Rv0678", "change": "S68S", "freq": "invalid"}]])
def test_malformed_native_export_returns_a_validation_error(variants: list) -> None:
    payload = json.dumps({"id": "TEST", "dr_variants": variants})
    with TestClient(app) as client:
        result = client.post("/review", data={"csv": payload})
        assert result.status_code == 422


def test_batch_prevalence_deduplication_uses_chronological_baseline() -> None:
    records = [
        {"sample_id": "LATER", "patient_id": "P", "collection_date": "Day 28", "probability_train": 0.8},
        {"sample_id": "FIRST", "patient_id": "P", "collection_date": "Day 0", "probability_train": 0.2},
    ]
    assert deduplicate_patient_probabilities(records) == [0.2]


def test_timeline_uses_dated_visits_and_excludes_synonymous_rows() -> None:
    records = [
        {"sample_id": "UNDATED", "patient_id": "P", "variants": []},
        {"sample_id": "FIRST", "patient_id": "P", "collection_date": "Day 0",
         "variants": [{"gene": "Rv0678", "mutation": "S68S", "frs": 1.0}]},
        {"sample_id": "SECOND", "patient_id": "P", "collection_date": "Day 28",
         "variants": []},
    ]
    timeline = analyze_patient_timelines(records)["P"]
    assert [visit["sample_id"] for visit in timeline["timepoints"]] == ["FIRST", "SECOND"]
    assert timeline["variants"] == []
    assert timeline["wt_series"] == [1.0, 1.0]


def test_example_timeline_opens_with_selectable_dated_visits() -> None:
    with TestClient(app) as client:
        result = client.post("/longitudinal-example")
        assert result.status_code == 200
        assert "Illustrative example" in result.text
        assert "Day 0" in result.text
        assert "Day 28" in result.text
        assert "Day 84" in result.text
        assert 'data-visit-index="0"' in result.text
