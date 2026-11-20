from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from fastapi.testclient import TestClient
from sklearn.metrics import average_precision_score, roc_auc_score

from discern.api.app import app
from discern.features import vocab_from
from discern.grammar import classify
from discern.prevalence import adjust, ece
from discern.research import BANNER, Archive, Reports


@pytest.mark.parametrize(("mutation", "flag"), [
    ("R50Q", "is_mis"), ("141_ins_c", "is_indel"), ("193_del_g", "is_indel"),
    ("135_minorindel", "is_indel"), ("Q131!", "is_stop"), ("c-11a", "is_prom"),
    ("S68S", "IS_SYN"), ("S68Z", "IS_NULL"),
])
def test_documented_grammar_flags(mutation: str, flag: str) -> None:
    assert classify(mutation)[flag] is True


def test_saved_cohort_and_baseline_metrics_agree_with_acceptance() -> None:
    reports = Reports(Path("."))
    prep = reports.observed("prep")
    samples = pd.read_parquet("data/processed/samples.parquet")
    mutations = pd.read_parquet("data/processed/mutations.parquet")
    assert (len(samples), samples.y.sum(), samples.SUBJID.nunique(), len(mutations)) == (
        prep["n_samples"], prep["n_resistant"], prep["unique_subjid"], prep["mutation_rows"],
    )
    for label, target in (("0.75", 0.497), ("0.1", 0.679)):
        record = reports.get(f"who_baseline/e02_who_as_is_frs{label}")
        assert record["observed"]["sens_U_as_S"] == pytest.approx(target, abs=0.005)


def test_patient_folds_are_disjoint_and_archived_scores_are_calibrated() -> None:
    patients = pd.read_parquet("data/processed/samples.parquet", columns=["UNIQUEID", "SUBJID"])
    predictions = pd.read_parquet("results/oof_predictions.parquet").merge(
        patients, on="UNIQUEID", validate="one_to_one",
    )
    assert predictions.groupby("SUBJID").fold.nunique().eq(1).all()
    assert ece(predictions.y, predictions.p) < 0.03
    np.testing.assert_allclose(adjust(np.array([0.1, 0.35, 0.85]), 0.073, 0.073), [0.1, 0.35, 0.85])
    saved = Reports(Path(".")).observed("model")
    assert roc_auc_score(predictions.y, predictions.p) == pytest.approx(saved["auc"], abs=1e-12)
    assert average_precision_score(predictions.y, predictions.p) == pytest.approx(saved["auprc"], abs=1e-12)


def test_queue_evaluation_accounts_for_revealed_and_unrevealed_samples() -> None:
    report = Path("results/approaches/voi_queue/e15_label_reveal_sim/resolution_curves.csv")
    curves = pd.read_csv(report)
    assert curves.n_unrevealed_counted_resolved.ge(0).all()
    assert (curves.n_resolved == curves.n_samples_revealed + curves.n_unrevealed_counted_resolved).all()


def test_vocabulary_ignores_variants_present_only_in_held_out_rows() -> None:
    core = pd.DataFrame({"seen": [1, 1, 1, 0, 0, 0], "held_out_only": [0, 0, 0, 1, 1, 1]})
    assert vocab_from(core, np.array([True, True, True, False, False, False]), 3) == ["seen"]


def test_public_pages_and_reports_have_the_banner_and_registered_values() -> None:
    source = Archive(Path("."))
    with TestClient(app) as client:
        routes = ["/", "/queue", "/priority", "/performance", "/model-card", "/about", "/health",
                  f"/result?sample_id={source.example_id}"]
        for path in routes:
            response = client.get(path)
            assert response.status_code == 200, path
            assert BANNER in response.text
        for report in source.reports.paths:
            response = client.get("/performance", params={"report": report})
            assert response.status_code == 200, report
            assert response.headers["X-Research-Use"] == BANNER
        auc = source.reports.observed("model")["auc"]
        assert str(auc) in client.get("/performance", params={"report": "lr_model/e04_s1_patient_cv"}).text
        shift = client.get("/performance", params={"report": "prevalence_correction/e10_prevalence_shift"})
        assert "Adjusting scores with known target prevalence" in shift.text
        assert "The target prevalence came from recorded lab labels" in shift.text
        queue = client.get("/performance", params={
            "report": "voi_queue/e15_label_reveal_sim", "table": "resolution_curves.csv",
        })
        assert "Retrospective comparison of variant sets" in queue.text
        assert "unrevealed counted as resolved" in queue.text.lower()
        assert client.get("/evidence/pepQ/G197R").json()["banner"] == BANNER
        assert client.post("/predict", json={"variants": []}).status_code == 404
