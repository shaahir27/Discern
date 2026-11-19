"""Tier B model, nested calibration, thresholds (AGENTS.md section 9). Golden implementation."""
from __future__ import annotations

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve
from sklearn.model_selection import GroupKFold


def logit(p):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))


def fit_nested(X, y, groups, C=0.3, max_iter=4000):
    n_groups = len(np.unique(groups))
    oof = np.zeros(len(y))
    for tr, te in GroupKFold(n_splits=min(5, n_groups)).split(X, y, groups):
        lr = LogisticRegression(C=C, penalty="l2", max_iter=max_iter).fit(X[tr], y[tr])
        oof[te] = lr.predict_proba(X[te])[:, 1]
    platt = LogisticRegression(C=1e6).fit(logit(oof).reshape(-1, 1), y)
    poof = platt.predict_proba(logit(oof).reshape(-1, 1))[:, 1]
    final = LogisticRegression(C=C, penalty="l2", max_iter=max_iter).fit(X, y)
    return final, platt, poof


def score(final, platt, X):
    raw = final.predict_proba(X)[:, 1]
    return raw, platt.predict_proba(logit(raw).reshape(-1, 1))[:, 1]


def thresholds(poof, y, npv=0.985, ppv=0.80, min_n=20):
    cands = np.unique(np.round(poof, 4))
    best, tb = -9, 0.5
    for t in cands:
        pr = poof >= t
        sens = (pr & (y == 1)).sum() / max((y == 1).sum(), 1)
        spec = ((~pr) & (y == 0)).sum() / max((y == 0).sum(), 1)
        if sens + spec - 1 > best:
            best, tb = sens + spec - 1, t
    lo = None
    for t in cands[::-1]:
        sel = poof <= t
        if sel.sum() >= min_n and (1 - y[sel].mean()) >= npv:
            lo = t
            break
    hi = None
    for t in cands:
        sel = poof >= t
        if sel.sum() >= min_n and y[sel].mean() >= ppv:
            hi = t
            break
    lo = -1.0 if lo is None else lo
    hi = 2.0 if hi is None else hi
    hi = max(hi, lo + 1e-6)
    return {"t_bin": float(tb), "lo": float(lo), "hi": float(hi)}


def sens_at_spec(y, p, spec_t):
    fpr, tpr, _ = roc_curve(y, p)
    ok = (1 - fpr) >= spec_t
    return float(tpr[ok].max())