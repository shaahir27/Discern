"""Shared nested out-of-fold runner for arbitrary fold lists."""
from __future__ import annotations

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold

from .features import design, vocab_from
from .model import fit_nested, logit, score, thresholds


def nested_oof(core, a, y, g, folds, fset="F2", C=0.3, max_iter=4000,
               npv=0.985, ppv=0.80, min_n=20, vocab_min=3):
    """folds: iterable of (train_idx, test_idx). Returns calibrated OOF P, ternary T, thresholds."""
    P = np.zeros(len(y))
    T = np.empty(len(y), dtype=object)
    th = []
    for tr, te in folds:
        vocab = vocab_from(core, tr, vocab_min)
        X = design(core, a, vocab, fset).values
        final, platt, poof = fit_nested(X[tr], y[tr], g[tr], C, max_iter)
        t = thresholds(poof, y[tr], npv, ppv, min_n)
        th.append(t)
        _, p = score(final, platt, X[te])
        P[te] = p
        T[te] = np.where(p <= t["lo"], "S", np.where(p >= t["hi"], "R", "U"))
    return P, T, th


def _raw_score(est, X):
    """Continuous log-odds-like score: logit(proba) for probabilistic models, margin for SVM."""
    if hasattr(est, "predict_proba"):
        return logit(est.predict_proba(X)[:, 1])
    return est.decision_function(X)


def fit_nested_est(X, y, groups, make_est):
    """Golden nested protocol generalized to any estimator factory (Platt on the raw score)."""
    n_groups = len(np.unique(groups))
    oof = np.zeros(len(y))
    for tr, te in GroupKFold(n_splits=min(5, n_groups)).split(X, y, groups):
        est = make_est()
        est.fit(X[tr], y[tr])
        oof[te] = _raw_score(est, X[te])
    platt = LogisticRegression(C=1e6).fit(oof.reshape(-1, 1), y)
    poof = platt.predict_proba(oof.reshape(-1, 1))[:, 1]
    final = make_est().fit(X, y)
    return final, platt, poof


def score_est(final, platt, X):
    raw = _raw_score(final, X)
    return raw, platt.predict_proba(raw.reshape(-1, 1))[:, 1]


def nested_oof_est(core, a, y, g, folds, make_est, fset="F2",
                   npv=0.985, ppv=0.80, min_n=20, vocab_min=3):
    """Like nested_oof but with an arbitrary estimator. folds: iterable of (train_idx, test_idx)."""
    P = np.zeros(len(y))
    T = np.empty(len(y), dtype=object)
    th = []
    for tr, te in folds:
        vocab = vocab_from(core, tr, vocab_min)
        X = design(core, a, vocab, fset).values
        final, platt, poof = fit_nested_est(X[tr], y[tr], g[tr], make_est)
        t = thresholds(poof, y[tr], npv, ppv, min_n)
        th.append(t)
        _, p = score_est(final, platt, X[te])
        P[te] = p
        T[te] = np.where(p <= t["lo"], "S", np.where(p >= t["hi"], "R", "U"))
    return P, T, th