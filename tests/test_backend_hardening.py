import csv
import io
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
from fastapi.testclient import TestClient
from pydantic import ValidationError

from discern.api.app import app
from discern.api.schemas import PredictRequest, Variant
from discern.decision import Predictor
from discern.research import MAX_CSV_BYTES, Archive


def test_csv_reports_location_without_echoing_a_submitted_value() -> None:
    with TestClient(app) as client:
        rows = list(csv.DictReader(io.StringIO(client.get("/example.csv").text)))
        rows[0]["frs"] = "private-value"
        output = io.StringIO()
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
        response = client.post("/batch", content=output.getvalue(), headers={"Content-Type": "text/csv"})
        assert response.status_code == 422
        assert response.json()["error"] == {"code": "invalid_fraction", "line": 2, "column": "frs"}
        assert "private-value" not in response.text
        form = client.post("/review", data={"csv": output.getvalue()})
        assert form.status_code == 422
        assert "CSV line 2" in form.text
        assert "private-value" not in form.text


def test_upload_limits_preserve_status_and_error_category() -> None:
    with TestClient(app) as client:
        response = client.post("/review", data={"csv": "+" * (MAX_CSV_BYTES + 1)})
        assert response.status_code == 413
        response = client.post("/batch", content=b"a" * (MAX_CSV_BYTES + 1), headers={"Content-Type": "text/csv"})
        assert response.status_code == 413
        assert response.json()["error"]["code"] == "file_too_large"


def test_replay_exposes_reference_prevalence_and_marks_missing_explanations() -> None:
    data = Archive(Path("."))
    replay = data.replay(data.example_id)
    oof = pd.read_parquet("results/oof_predictions.parquet")
    training = oof.loc[oof.fold != replay["fold"]]
    assert replay["score_context"]["reference_prevalence"] == pytest.approx(training.y.mean())
    assert replay["score_context"]["uncertainty"] is None
    assert replay["call_explanation"]["contributing_variants"] is None


def test_predictive_values_use_call_rates_across_all_archived_cases() -> None:
    data = Archive(Path("."))
    oof = pd.read_parquet("results/oof_predictions.parquet")
    prevalence = 0.05
    weights = oof.y.map({1: prevalence / oof.y.mean(), 0: (1 - prevalence) / (1 - oof.y.mean())})
    expected_ppv = weights[oof.cat_call.eq("R") & oof.y.eq(1)].sum() / weights[oof.cat_call.eq("R")].sum()
    expected_npv = weights[oof.cat_call.eq("S") & oof.y.eq(0)].sum() / weights[oof.cat_call.eq("S")].sum()
    with TestClient(app) as client:
        page = client.get("/result", params={"sample_id": data.example_id, "prevalence": prevalence})
        assert f"{expected_ppv:.1%}" in page.text
        assert f"{expected_npv:.1%}" in page.text


def test_legacy_predictive_values_use_the_archived_e16_policy() -> None:
    rates = pd.read_csv("results/models/lr_mechanism/e16_lof_ablation/s1_veto_metrics.csv").set_index("state").loc["ON"]
    prevalence = 0.05
    true_r = prevalence * rates.sens_U_as_S
    false_r = (1 - prevalence) * (1 - rates.spec_U_as_S)
    true_s = (1 - prevalence) * rates.spec_U_as_R
    false_s = prevalence * (1 - rates.sens_U_as_R)
    values = Predictor().predict([], prevalence=prevalence)["predictive_values"]
    assert values["ppv_if_R"] == pytest.approx(true_r / (true_r + false_r))
    assert values["npv_if_S"] == pytest.approx(true_s / (true_s + false_s))


def test_legacy_schema_rejects_free_text_and_mismatched_prefix() -> None:
    for entry in ("random text r", "garbage!", "pepQ@R50Q", "Rv0678@Rv0678@R50Q"):
        with pytest.raises(ValidationError):
            Variant(gene="Rv0678", mutation=entry)
    with pytest.raises(ValidationError):
        PredictRequest(sample_id="valid\n")
    assert Variant(gene="Rv0678", mutation="Rv0678@R50Q").mutation == "R50Q"


@pytest.mark.parametrize("returned", [None, {}])
def test_catalogue_failure_defers_the_call_and_flags_the_error(
    monkeypatch: pytest.MonkeyPatch, returned: dict | None,
) -> None:
    data = Archive(Path("."))
    variants = data.profile(data.example_id)
    predictor = Predictor()

    def unavailable(_variant: str) -> dict:
        if returned is None:
            raise RuntimeError("private detail")
        return returned

    monkeypatch.setattr(predictor, "rc", SimpleNamespace(predict=unavailable))
    result = predictor.predict(variants)
    assert result["call"] == "U"
    assert "CATALOGUE_EVALUATION_FAILED" in result["warnings"]
    assert "private detail" not in str(result)
