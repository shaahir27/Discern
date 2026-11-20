"""Unit tests for Batch-Recurrence & Contamination Surveillance Audit."""
from __future__ import annotations

from pathlib import Path

import pytest

from discern.evidence import EvidenceStore
from discern.surveillance import batch_recurrence_audit


@pytest.fixture
def ev_store():
    return EvidenceStore(Path("models/lr_model/evidence_store.csv"))


def test_same_patient_serial_isolates_do_not_trigger(ev_store):
    # 3 serial samples from the same patient carrying rare Rv0678@L32S
    rows_by_sample = {
        "SMP1": [{"gene": "Rv0678", "mutation": "L32S", "frs": 0.85}],
        "SMP2": [{"gene": "Rv0678", "mutation": "L32S", "frs": 0.90}],
        "SMP3": [{"gene": "Rv0678", "mutation": "L32S", "frs": 0.92}],
    }
    patient_map = {"SMP1": "PT_SAME", "SMP2": "PT_SAME", "SMP3": "PT_SAME"}

    alerts = batch_recurrence_audit(rows_by_sample, ev_store, patient_map=patient_map)
    assert len(alerts) == 0


def test_two_distinct_patients_rare_variant_triggers_wgs_escalation(ev_store):
    # 2 distinct patients sharing rare Rv0678@L32S
    rows_by_sample = {
        "SMP_A": [{"gene": "Rv0678", "mutation": "L32S", "frs": 0.85}],
        "SMP_B": [{"gene": "Rv0678", "mutation": "L32S", "frs": 0.90}],
    }
    patient_map = {"SMP_A": "PT_A", "SMP_B": "PT_B"}

    alerts = batch_recurrence_audit(rows_by_sample, ev_store, patient_map=patient_map)
    assert len(alerts) == 1
    alert = alerts[0]
    assert alert["category"] == "RARE_VARIANT_WGS_ESCALATION"
    assert alert["gene"] == "Rv0678"
    assert alert["mutation"] == "L32S"
    assert alert["n_patients"] == 2
    assert "SMP_A" in alert["sample_ids"]
    assert "SMP_B" in alert["sample_ids"]
    # Radical honesty guardrail: must NOT claim confirmed transmission
    assert "confirmed transmission cluster" not in alert["message"].lower()
    assert "escalate to wgs" in alert["message"].lower()


def test_two_patients_homopolymer_triggers_qc_warning(ev_store):
    # 2 distinct patients sharing homopolymer indel 192_ins_g
    rows_by_sample = {
        "SMP_X": [{"gene": "Rv0678", "mutation": "192_ins_g", "frs": 0.22}],
        "SMP_Y": [{"gene": "Rv0678", "mutation": "192_ins_g", "frs": 0.28}],
    }
    patient_map = {"SMP_X": "PT_X", "SMP_Y": "PT_Y"}

    alerts = batch_recurrence_audit(rows_by_sample, ev_store, patient_map=patient_map)
    assert len(alerts) == 1
    alert = alerts[0]
    assert alert["category"] == "QC_CONTAMINATION_OR_ARTIFACT_RISK"
    assert "homopolymer" in alert["message"].lower() or "cross-contamination" in alert["message"].lower()
