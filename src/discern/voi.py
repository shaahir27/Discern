"""G4 / P-E15: value-of-information lab queue.

score(set) = n_U_carriers x expected beta-binomial entropy reduction after revealing one label.
Retrospective label-reveal simulation: reveal top-K sets' existing labels, unlock variants via a
binomial-test proxy for the catalogue rule, measure %U resolved. Baselines: count-only, random.
Retrospective (assumes catalogued variants transfer across lineages); not prospective.
"""
from __future__ import annotations

import numpy as np
from scipy.special import betaln, digamma
from scipy.stats import binomtest


def beta_entropy(a: float, b: float) -> float:
    return float(betaln(a, b) - (a - 1) * digamma(a) - (b - 1) * digamma(b) + (a + b - 2) * digamma(a + b))


def expected_info_gain(a: float, b: float) -> float:
    """Mutual information between one Bernoulli label and the Beta(a,b) resistance probability."""
    q = a / (a + b)
    return beta_entropy(a, b) - (q * beta_entropy(a + 1, b) + (1 - q) * beta_entropy(a, b + 1))


def variant_unlock(n: int, k: int, background: float = 0.1, alpha: float = 0.05, min_n: int = 3) -> str:
    """Binomial-test proxy for the catalogue unlock rule (documented approximation of strict_unlock)."""
    if n < min_n:
        return "U"
    p = binomtest(int(k), int(n), background).pvalue
    if p < alpha:
        return "R" if k / n > background else "S"
    return "U"


def rank_sets(set_df, strategy: str, seed: int = 0):
    """set_df: index=key, columns n (U carriers), eig (per-label info gain). Return ordered keys."""
    if strategy == "count":
        return list(set_df.sort_values("n", ascending=False).index)
    if strategy == "info_gain":
        return list(set_df.assign(score=set_df.n * set_df.eig).sort_values("score", ascending=False).index)
    if strategy == "random":
        rng = np.random.default_rng(seed)
        keys = list(set_df.index)
        rng.shuffle(keys)
        return keys
    raise ValueError(strategy)