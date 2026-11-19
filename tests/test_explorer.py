import csv
import io

import pytest
from fastapi.testclient import TestClient


def test_empty_profile_keeps_the_research_warning() -> None:
    from discern.decision import Predictor

    result = Predictor().predict([], sample_id="EMPTY")
    assert result["call"] == "S"
    assert "NO_VARIANT_INVISIBLE_RESISTANCE_POSSIBLE" in result["warnings"]


def test_replay_returns_the_saved_result_and_batch_accepts_recorded_ids_only(caplog: pytest.LogCaptureFixture) -> None:
    from discern.api.app import app

    with TestClient(app) as client:
        roster = list(csv.DictReader(io.StringIO(client.get("/replay-samples.csv").text)))
        sample_id = roster[0]["sample_id"]
        replay = client.get(f"/replay/{sample_id}").json()
        assert replay["sample_id"] == sample_id
        assert replay["mode"] == "held-out in cross-validation"
        assert replay["banner"]
        batch = client.post("/batch", content=f"sample_id\n{sample_id}\n", headers={"Content-Type": "text/csv"})
        assert batch.status_code == 200
        assert next(iter(csv.DictReader(io.StringIO(batch.text))))["call"] == replay["call"]
        rejected = client.post(
            "/batch",
            content="sample_id,mutation\nprivate-id,private-payload\n",
            headers={"Content-Type": "text/csv"},
        )
        assert rejected.status_code == 422
        assert "private-payload" not in rejected.text
        assert "private-payload" not in caplog.text


def test_csv_upload_opens_the_saved_score_and_rejects_an_edited_profile() -> None:
    from discern.api.app import app

    with TestClient(app) as client:
        template = client.get("/example.csv")
        assert template.status_code == 200
        rows = list(csv.DictReader(io.StringIO(template.text)))
        replay = client.get(f"/replay/{rows[0]['sample_id']}").json()
        result = client.post("/review", data={"csv": template.text})
        assert result.status_code == 200
        assert "Saved model score" in result.text
        assert f"{replay['probability']:.1%}" in result.text
        assert "Bedaquiline result" in result.text
        assert "Recorded lab result" in result.text
        assert "Observed gene changes" in result.text
        assert "The archive did not record which variants drove this call." in result.text
        assert result.text.count('class="panel variant-card"') == len(replay["evidence"])
        for item in replay["evidence"]:
            gene, mutation = item["variant"].split("@")
            assert f'<span class="variant-gene">{gene}</span>' in result.text
            assert f'<span class="variant-change">{mutation}</span>' in result.text
            if item["n_carriers"]:
                assert f'{item["n_resistant"]} / {item["n_carriers"]}' in result.text
        rows[0]["frs"] = "0.004321"
        output = io.StringIO()
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
        rejected = client.post("/review", data={"csv": output.getvalue()})
        assert rejected.status_code == 422
        assert "could not match" in rejected.text
        assert "0.004321" not in rejected.text
