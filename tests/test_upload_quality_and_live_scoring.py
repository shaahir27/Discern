"""Tests for Live CSV Quality Check and live model scoring on uploaded CSVs."""

import csv
import io

import pytest
from fastapi.testclient import TestClient

from discern.api.app import app
from discern.quality import check_csv_quality


def test_quality_check_calculates_fresh_completeness_and_identifies_issues() -> None:
    # 1. Clean valid CSV (100% complete)
    clean_csv = """sample_id,gene,mutation,frs
NEW_01,Rv0678,141_ins_c,1.0
NEW_01,atpE,A63V,0.85
"""
    rep = check_csv_quality(clean_csv)
    assert rep["is_valid"] is True
    assert rep["completeness_score"] == 100.0
    assert rep["total_rows"] == 2
    assert rep["valid_rows"] == 2
    assert rep["missing_count"] == 0
    assert rep["duplicate_count"] == 0
    assert rep["formatting_count"] == 0

    # 2. CSV with missing values, duplicate row, and invalid formatting
    problem_csv = """sample_id,gene,mutation,frs
NEW_02,,141_ins_c,1.0
NEW_02,Rv0678,141_ins_c,1.0
NEW_02,Rv0678,141_ins_c,1.0
NEW_03,fake_gene,invalid_syntax_text,1.5
"""
    bad_rep = check_csv_quality(problem_csv)
    assert bad_rep["is_valid"] is False
    assert bad_rep["completeness_score"] < 100.0
    assert bad_rep["missing_count"] == 1
    assert bad_rep["duplicate_count"] == 1
    # fake_gene, invalid_syntax_text, and frs 1.5 are formatting issues
    assert bad_rep["formatting_count"] >= 3

    # Check that rows with issues are flagged by line
    types = [issue["type"] for issue in bad_rep["issues"]]
    assert "missing" in types
    assert "duplicate" in types
    assert "format" in types


def test_live_scoring_runs_on_new_uploaded_csv() -> None:
    with TestClient(app) as client:
        # A completely new sample ID not in the pre-trained replay roster
        new_sample_csv = """sample_id,gene,mutation,frs
DEMO_LIVE_ISOLATE,Rv0678,141_ins_c,1.0
"""
        response = client.post("/review", data={"csv": new_sample_csv})
        assert response.status_code == 200

        # Check that Upload quality is present and labeled
        assert "Upload quality" in response.text
        assert "Completeness" in response.text
        assert "100.0%" in response.text

        # Check that the model ran and produced a resistance score
        assert "Bedaquiline result" in response.text
        assert "Catalogue call" in response.text
        assert "Resistant" in response.text
        assert "Live model evaluation" in response.text
        assert "DEMO_LIVE_ISOLATE" in response.text

        # Verify separation: Upload quality vs resistance results
        assert "Kept separate from archived resistance results" in response.text


def test_live_scoring_batch_queue_displays_upload_quality() -> None:
    with TestClient(app) as client:
        batch_csv = """sample_id,gene,mutation,frs
ISOLATE_ALPHA,Rv0678,141_ins_c,1.0
ISOLATE_BETA,Rv0678,c-11a,0.95
"""
        response = client.post("/review", data={"csv": batch_csv})
        assert response.status_code == 200
        assert "Review queue" in response.text
        assert "Upload quality" in response.text
        assert "ISOLATE_ALPHA" in response.text
        assert "ISOLATE_BETA" in response.text
