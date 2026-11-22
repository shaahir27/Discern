"""Verification tests for Backend API & Wiring improvements (B1-B6)."""
import csv
import io
from pathlib import Path

from fastapi.testclient import TestClient

from discern.api.app import app
from discern.research import BANNER, Archive, Reports


def test_b1_predict_endpoint_success_and_validation():
    with TestClient(app) as client:
        # Empty profile -> S call, invisible resistance warning
        r0 = client.post("/predict", json={"variants": []})
        assert r0.status_code == 200
        d0 = r0.json()
        assert d0["call"] == "S"
        assert d0["banner"] == BANNER
        assert "NO_VARIANT_INVISIBLE_RESISTANCE_POSSIBLE" in d0["warnings"]
        assert "model" in d0
        assert "tier_a" in d0
        assert "predictive_values" in d0
        assert "evidence" in d0
        assert "versions" in d0

        # Known resistant variant -> R call
        r1 = client.post("/predict", json={
            "variants": [{"gene": "Rv0678", "mutation": "141_ins_c", "frs": 1.0}],
        })
        assert r1.status_code == 200
        d1 = r1.json()
        assert d1["call"] == "R"
        assert d1["tier_a"]["variant_calls"][0]["call"] == "R"

        # Unresolved variant -> U call
        r2 = client.post("/predict", json={
            "variants": [{"gene": "pepQ", "mutation": "G197R", "frs": 1.0}],
        })
        assert r2.status_code == 200
        d2 = r2.json()
        assert d2["call"] != "R"

        # Validation: invalid gene -> 422
        r_gene = client.post("/predict", json={"variants": [{"gene": "unknown_gene", "mutation": "R50Q"}]})
        assert r_gene.status_code == 422

        # Validation: invalid mutation notation -> 422
        r_mut = client.post("/predict", json={"variants": [{"gene": "Rv0678", "mutation": "free_text_bad"}]})
        assert r_mut.status_code == 422

        # Validation: invalid prevalence -> 422
        r_prev = client.post("/predict", json={"prevalence": 0.99, "variants": []})
        assert r_prev.status_code == 422


def test_b2_batch_live_upload_and_replay():
    with TestClient(app) as client:
        # Live multi-sample upload
        csv_payload = (
            "sample_id,gene,mutation,frs\n"
            "LIVE_S1,Rv0678,141_ins_c,1.0\n"
            "LIVE_S2,pepQ,G197R,1.0\n"
        )
        resp = client.post("/batch", content=csv_payload, headers={"Content-Type": "text/csv"})
        assert resp.status_code == 200
        assert "text/csv" in resp.headers["content-type"]
        reader = list(csv.DictReader(io.StringIO(resp.text)))
        assert len(reader) == 2
        calls = {row["sample_id"]: row["call"] for row in reader}
        assert calls["LIVE_S1"] == "R"
        assert BANNER in resp.text

        # Replay batch with recorded sample IDs
        data = Archive(Path("."))
        replay_payload = f"sample_id\n{data.example_id}\n"
        resp_replay = client.post("/batch", content=replay_payload, headers={"Content-Type": "text/csv"})
        assert resp_replay.status_code == 200
        reader_replay = list(csv.DictReader(io.StringIO(resp_replay.text)))
        assert len(reader_replay) == 1
        assert reader_replay[0]["sample_id"] == data.example_id

        # Unknown sample ID in replay-style batch -> 422
        bad_replay = "sample_id\nUNKNOWN_REPLAY_SAMPLE_ID\n"
        resp_bad = client.post("/batch", content=bad_replay, headers={"Content-Type": "text/csv"})
        assert resp_bad.status_code == 422

        # Unsupported media type -> 415
        resp_bad_type = client.post("/batch", content=csv_payload, headers={"Content-Type": "application/json"})
        assert resp_bad_type.status_code == 415


def test_b3_replay_domain_map_context():
    data = Archive(Path("."))
    replay = data.replay(data.example_id)
    assert any("pos" in item and item["pos"] is not None for item in replay["evidence"])
    assert any("lof_annotation" in item and item["lof_annotation"] is not None for item in replay["evidence"])

    with TestClient(app) as client:
        res = client.get(f"/result?sample_id={data.example_id}")
        assert res.status_code == 200
        # Check that the protein domain track or SVG domain elements render
        assert "Rv0678 protein domain map" in res.text or "DNA-Binding" in res.text or "domain-track" in res.text or "dimerization" in res.text.lower()


def test_b4_reports_in_memory_memoization():
    reports = Reports(Path("."))
    assert len(reports._cache) == 0
    cat1 = reports.get("catalogue/e03_s1_catalogue_cv")
    assert "catalogue/e03_s1_catalogue_cv" in reports._cache
    cat2 = reports.get("catalogue/e03_s1_catalogue_cv")
    assert cat1 is cat2


def test_b5_evidence_endpoint_validation():
    with TestClient(app) as client:
        # Existing variant
        r_ok = client.get("/evidence/pepQ/G197R")
        assert r_ok.status_code == 200
        assert r_ok.json()["variant"] == "pepQ@G197R"

        # Missing / unobserved variant -> 404
        r_missing = client.get("/evidence/Rv0678/A999A")
        assert r_missing.status_code == 404

        # Invalid gene -> 422
        r_bad_gene = client.get("/evidence/bad_gene/R50Q")
        assert r_bad_gene.status_code == 422

        # Invalid mutation syntax -> 422
        r_bad_mut = client.get("/evidence/Rv0678/invalid!!")
        assert r_bad_mut.status_code == 422

        # Gene prefix mismatch -> 422
        r_mismatch = client.get("/evidence/Rv0678/atpE@A999A")
        assert r_mismatch.status_code == 422


def test_b6_content_negotiated_errors():
    with TestClient(app) as client:
        # HTML error for browser navigation
        r_html = client.get("/nonexistent-endpoint", headers={"Accept": "text/html"})
        assert r_html.status_code == 404
        assert "text/html" in r_html.headers["content-type"]
        assert "<section class=\"error-state\"" in r_html.text
        assert BANNER in r_html.text

        # JSON error for API client
        r_json = client.get("/nonexistent-endpoint", headers={"Accept": "application/json"})
        assert r_json.status_code == 404
        assert "application/json" in r_json.headers["content-type"]
        assert r_json.json()["detail"] == "Not Found"
        assert r_json.json()["banner"] == BANNER

        # Default accept header (*/*) returns JSON
        r_default = client.get("/nonexistent-endpoint")
        assert r_default.status_code == 404
        assert "application/json" in r_default.headers["content-type"]
