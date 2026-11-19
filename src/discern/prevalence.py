"""Prevalence correction + calibration metrics (AGENTS.md section 10)."""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np


def odds(p):
    p = np.clip(p, 1e-9, 1 - 1e-9)
    return p / (1 - p)


def adjust(p, pi_tr, pi):
    o = odds(p) / odds(pi_tr) * odds(pi)
    return o / (1 + o)


def brier(y, p):
    return float(np.mean((np.asarray(p) - np.asarray(y)) ** 2))


def ece(y, p, bins=10):
    """Expected calibration error, standard quantile binning (sklearn calibration_curve convention)."""
    y = np.asarray(y, float)
    p = np.clip(np.asarray(p, float), 0, 1)
    edges = np.percentile(p, np.linspace(0, 100, bins + 1))
    edges[0] = 0.0
    edges[-1] = 1.0
    idx = np.searchsorted(edges[1:-1], p, side="right")
    err = 0.0
    n = len(p)
    for b in range(bins):
        sel = idx == b
        if sel.sum() == 0:
            continue
        err += sel.sum() / n * abs(y[sel].mean() - p[sel].mean())
    return float(err)


def prevalence_weights(y, pi_tr, pi):
    y = np.asarray(y)
    return np.where(y == 1, pi / pi_tr, (1 - pi) / (1 - pi_tr))


def ppv_npv(sens, spec, pi):
    ppv = sens * pi / (sens * pi + (1 - spec) * (1 - pi))
    npv = spec * (1 - pi) / (spec * (1 - pi) + (1 - sens) * pi)
    return float(ppv), float(npv)


def archived_predictive_values(metrics: Mapping[str, float], prevalence: float) -> tuple[float, float]:
    """Reweight actual resistant and susceptible call rates across the entire cohort."""
    true_r = prevalence * metrics["sens_U_as_S"]
    false_r = (1 - prevalence) * (1 - metrics["spec_U_as_S"])
    true_s = (1 - prevalence) * metrics["spec_U_as_R"]
    false_s = prevalence * (1 - metrics["sens_U_as_R"])
    return true_r / (true_r + false_r), true_s / (true_s + false_s)
