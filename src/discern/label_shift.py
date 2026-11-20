"""Unsupervised local prevalence estimation using Maximum Likelihood Label Shift (MLLS).

Applies the Saerens-Latinne-Decaestecker / Alexandari et al. EM fixed-point solver
on Platt-calibrated Tier B probabilities to recover true local bedaquiline resistance
prevalence from unlabeled batch sequences with bootstrap confidence intervals.
"""
from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

MIN_BATCH_N = 30
PI_MIN = 1e-4
PI_MAX = 0.99
MAX_ITER = 200
TOL = 1e-7


def _mlls_em(p: np.ndarray, pi_tr: float) -> tuple[float, bool, int]:
    """Fixed-point EM solver for label shift likelihood maximization."""
    pi_hat = float(np.clip(pi_tr, PI_MIN, PI_MAX))
    for it in range(1, MAX_ITER + 1):
        w1 = pi_hat / pi_tr
        w0 = (1.0 - pi_hat) / (1.0 - pi_tr)
        denom = w1 * p + w0 * (1.0 - p)
        q = (w1 * p) / np.maximum(denom, 1e-12)
        pi_next = float(np.clip(q.mean(), PI_MIN, PI_MAX))
        if abs(pi_next - pi_hat) < TOL:
            return pi_next, True, it
        pi_hat = pi_next
    return pi_hat, False, MAX_ITER


def estimate_prevalence_mlls(
    probs_train_prior: Sequence[float],
    pi_tr: float,
    *,
    n_bootstrap: int = 200,
    seed: int = 42,
) -> dict[str, Any]:
    """Estimate true local prevalence from Platt-calibrated probabilities under label shift."""
    arr = np.asarray(probs_train_prior, dtype=float)
    arr = np.clip(arr, 1e-6, 1.0 - 1e-6)
    n = len(arr)

    if n < MIN_BATCH_N:
        return {
            "status": "insufficient_n",
            "n": n,
            "min_n": MIN_BATCH_N,
            "pi_hat": None,
            "ci95_low": None,
            "ci95_high": None,
            "converged": False,
            "n_iter": 0,
            "note": (
                f"Batch has {n} distinct patients (minimum {MIN_BATCH_N} required for reliable "
                "prevalence estimation). Using manual/default prevalence."
            ),
        }

    pi_hat, converged, n_iter = _mlls_em(arr, pi_tr)

    # Bootstrap 95% confidence intervals
    rng = np.random.default_rng(seed)
    boot_estimates = np.empty(n_bootstrap, dtype=float)
    for b in range(n_bootstrap):
        resample = rng.choice(arr, size=n, replace=True)
        b_pi, _, _ = _mlls_em(resample, pi_tr)
        boot_estimates[b] = b_pi

    ci_low = float(np.percentile(boot_estimates, 2.5))
    ci_high = float(np.percentile(boot_estimates, 97.5))

    return {
        "status": "estimated",
        "n": n,
        "min_n": MIN_BATCH_N,
        "pi_hat": round(pi_hat, 4),
        "ci95_low": round(ci_low, 4),
        "ci95_high": round(ci_high, 4),
        "converged": converged,
        "n_iter": n_iter,
        "method": "MLLS (Saerens/Alexandari EM on Platt-calibrated Tier B scores)",
    }


def deduplicate_patient_probabilities(records: Sequence[dict[str, Any]]) -> list[float]:
    """Extract exactly one baseline probability per distinct patient to avoid longitudinal bias."""
    seen_patients: set[str] = set()
    probs: list[float] = []
    for r in records:
        pid = str(r.get("patient_id") or r.get("sample_id") or "").strip()
        if not pid or pid in seen_patients:
            continue
        seen_patients.add(pid)
        p = r.get("probability_train")
        if p is None:
            p = r.get("probability")
        if p is not None:
            probs.append(float(p))
    return probs
