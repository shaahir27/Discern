"""G3 / P-E14: Learn-then-Test / Bentkus risk-controlled abstention certificates.

S-zone risk = fraction truly-resistant among S-calls (target alpha=0.015, NPV>=98.5%).
R-zone risk = fraction truly-susceptible among R-calls (target beta=0.20, PPV>=80%).
Bentkus p-value tests H0: risk >= target; certify at level delta. The certified UPPER BOUND
(inverted Bentkus) is the guarantee — never the point target. Patient grouping => effective
n ~ patients; marginal split-conformal validity only, not i.i.d. (state in model card).
Prevalence weighting is class-constant (Tibshirani 2019); at an estimated user prevalence the
certificate is approximate (Barber 2022 coverage-gap bound).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq
from scipy.stats import binom

E = np.e


def bentkus_pvalue(k: int, n: int, target: float) -> float:
    """P-value for H0: error-rate >= target, given k errors in n trials (Bentkus bound)."""
    if n <= 0:
        return 1.0
    return float(min(1.0, E * binom.cdf(k, n, target)))


def certified_upper_bound(k: int, n: int, delta: float) -> float:
    """Smallest risk level r rejectable at delta => certified error-rate <= r (conf 1-delta)."""
    if n <= 0 or k >= n:
        return 1.0
    f = lambda r: E * binom.cdf(k, n, r) - delta
    if f(1.0) > 0:
        return 1.0
    return float(brentq(f, 0.0, 1.0, xtol=1e-6))


def zone_counts(y, p, t, zone: str):
    y = np.asarray(y)
    p = np.asarray(p)
    sel = p <= t if zone == "S" else p >= t
    n = int(sel.sum())
    k = int((y[sel] == 1).sum()) if zone == "S" else int((y[sel] == 0).sum())
    return n, k


def certified_thresholds(p, y, alpha=0.015, beta=0.20, delta=0.05, min_n=20):
    """Least-conservative certified (t_S, t_R): largest certified t_S, smallest certified t_R>t_S.

    Returns (t_S, t_R, info). t_S=None / t_R=2.0 mean that zone could not be certified (abstain).
    """
    p = np.asarray(p)
    y = np.asarray(y)
    cands = np.unique(np.round(p, 4))
    tS, bound_s, emp_s = None, 1.0, None
    for t in cands:  # increasing: keep the largest certified t_S
        n, k = zone_counts(y, p, t, "S")
        if n >= min_n and bentkus_pvalue(k, n, alpha) <= delta:
            tS = float(t)
            bound_s = certified_upper_bound(k, n, delta)
            emp_s = k / n
    tR, bound_r, emp_r = 2.0, 1.0, None
    for t in cands:  # increasing: first certified t_R above t_S
        if tS is not None and t <= tS:
            continue
        n, k = zone_counts(y, p, t, "R")
        if n >= min_n and bentkus_pvalue(k, n, beta) <= delta:
            tR = float(t)
            bound_r = certified_upper_bound(k, n, delta)
            emp_r = k / n
            break
    info = {"t_S": tS, "t_R": tR, "bound_s": bound_s, "bound_r": bound_r,
            "emp_s": emp_s, "emp_r": emp_r,
            "s_certified": tS is not None, "r_certified": tR < 2.0}
    return tS, tR, info


def apply_policy(p, tS, tR):
    p = np.asarray(p)
    tS = -1.0 if tS is None else tS
    return np.where(p <= tS, "S", np.where(p >= tR, "R", "U"))


def zone_error(y, calls, zone: str):
    y = np.asarray(y)
    calls = np.asarray(calls)
    sel = calls == zone
    n = int(sel.sum())
    if n == 0:
        return None, 0
    k = int((y[sel] == 1).sum()) if zone == "S" else int((y[sel] == 0).sum())
    return k / n, n


def weighted_zone_error(y, p, t, zone, pi_tr, pi_user):
    """Class-constant prevalence-weighted empirical error (approximate certificate at pi_user)."""
    y = np.asarray(y)
    p = np.asarray(p)
    w = np.where(y == 1, pi_user / pi_tr, (1 - pi_user) / (1 - pi_tr))
    sel = p <= t if zone == "S" else p >= t
    nw = float(w[sel].sum())
    if nw <= 0:
        return None
    kw = float((w[sel] * (y[sel] == 1)).sum()) if zone == "S" else float((w[sel] * (y[sel] == 0)).sum())
    return kw / nw