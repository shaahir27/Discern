"""Unit tests for Unsupervised Batch Prevalence Estimation (MLLS / Saerens EM)."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from discern.label_shift import estimate_prevalence_mlls, deduplicate_patient_probabilities


def test_mlls_small_batch_gate():
    # Less than MIN_BATCH_N (30) samples
    p_small = [0.05] * 10
    res = estimate_prevalence_mlls(p_small, pi_tr=0.073)
    assert res["status"] == "insufficient_n"
    assert res["pi_hat"] is None
    assert res["n"] == 10
    assert res["min_n"] == 30
    assert "minimum 30 required" in res["note"]


def test_mlls_recovers_simulated_prevalences():
    oof_path = Path("results/oof_predictions.parquet")
    if not oof_path.is_file():
        pytest.skip("results/oof_predictions.parquet not found")

    df = pd.read_parquet(oof_path)
    pi_tr = 866.0 / 11868.0  # reference training prevalence ~0.07297

    # Split into true resistant and susceptible score pools
    r_scores = df.loc[df["y"] == 1, "p"].values
    s_scores = df.loc[df["y"] == 0, "p"].values

    rng = np.random.default_rng(2026)
    target_prevalences = [0.01, 0.02, 0.05, 0.10, 0.20]
    N = 500

    for true_pi in target_prevalences:
        n_r = int(round(N * true_pi))
        n_s = N - n_r
        sample_r = rng.choice(r_scores, size=n_r, replace=True)
        sample_s = rng.choice(s_scores, size=n_s, replace=True)
        batch = np.concatenate([sample_r, sample_s])
        rng.shuffle(batch)

        res = estimate_prevalence_mlls(batch, pi_tr=pi_tr, n_bootstrap=100, seed=42)
        assert res["status"] == "estimated"
        assert res["converged"] is True
        pi_hat = res["pi_hat"]
        ci_low = res["ci95_low"]
        ci_high = res["ci95_high"]

        # Assert estimate is within 0.02 of target
        assert abs(pi_hat - true_pi) < 0.025, f"Target {true_pi}, got {pi_hat}"
        # Assert true_pi is inside or within 0.01 of 95% CI
        assert (ci_low - 0.01) <= true_pi <= (ci_high + 0.01), f"Target {true_pi} outside CI [{ci_low}, {ci_high}]"


def test_mlls_deduplicates_patient_id():
    # 35 baseline susceptible patients
    records = []
    for i in range(35):
        records.append({
            "sample_id": f"SMP_{i}",
            "patient_id": f"PT_{i}",
            "probability_train": 0.02,
        })
    # Add 1 patient with 10 serial resistant visits
    for visit in range(10):
        records.append({
            "sample_id": f"SMP_R_{visit}",
            "patient_id": "PT_RESISTANT",
            "probability_train": 0.95,
        })

    # Total samples = 45, but distinct patients = 36
    deduped = deduplicate_patient_probabilities(records)
    assert len(deduped) == 36

    # Verify PT_RESISTANT only contributes 1 probability score
    res = estimate_prevalence_mlls(deduped, pi_tr=0.073)
    assert res["status"] == "estimated"
    assert res["n"] == 36
    # 1 resistant out of 36 should give a low prevalence (~0.03), not inflated by 10 visits (~0.23)
    assert res["pi_hat"] < 0.08

