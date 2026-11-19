"""G3 / M3: empirical-Bayes beta-binomial shrinkage of per-variant resistance proportions.

Rare variants are shrunk toward their mechanism-class mean (DBD/mid/late/full_del/other),
estimated by method of moments. A standard estimator choice (not an invention) — ablate, don't
headline. Produces raw vs shrunk proportions + Wilson CIs for evidence cards.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def wilson(k: int, n: int, z: float = 1.96):
    if n == 0:
        return 0.0, 1.0
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return float(max(0.0, center - half)), float(min(1.0, center + half))


def beta_prior_moments(k: np.ndarray, n: np.ndarray):
    """Method-of-moments Beta(a,b) prior from carrier counts k and totals n."""
    k = np.asarray(k, float)
    n = np.asarray(n, float)
    N = n.sum()
    if N == 0:
        return 1.0, 1.0
    mu = k.sum() / N
    mu = min(max(mu, 1e-6), 1 - 1e-6)
    p = np.divide(k, n, out=np.zeros_like(k), where=n > 0)
    w = n / N
    var_total = float(np.sum(w * (p - mu) ** 2))
    var_sampling = float(np.sum(w * mu * (1 - mu) / np.maximum(n, 1)))
    tau2 = max(var_total - var_sampling, 1e-6)
    conc = mu * (1 - mu) / tau2 - 1
    conc = max(conc, 1.0)
    return float(mu * conc), float((1 - mu) * conc)


def variant_shrinkage(m: pd.DataFrame, eff: pd.DataFrame | None = None) -> pd.DataFrame:
    """Per core-gene variant: n_carriers, n_resistant, raw vs class-shrunk proportion, Wilson CI."""
    mm = m[(~m.IS_SYN) & (~m.IS_NULL) & (m.FRSf >= 0.1)]
    if "y" not in mm.columns:
        raise ValueError("mutations must carry y (merge samples.y first)")
    grp = mm.groupby(["GENE", "MUT"]).agg(n=("UNIQUEID", "nunique"),
                                         k=("y", lambda s: int(s.sum()))).reset_index()
    cls_map = {}
    if eff is not None and len(eff):
        cls_map = eff.set_index("MUT").region.to_dict()
    grp["mech_class"] = grp.apply(
        lambda r: cls_map.get(r["MUT"], "other") if r["GENE"] == "Rv0678" else r["GENE"], axis=1)
    rows = []
    for cls, g in grp.groupby("mech_class"):
        a, b = beta_prior_moments(g.k.values, g.n.values)
        for _, r in g.iterrows():
            raw = r.k / r.n if r.n else 0.0
            shrunk = (r.k + a) / (r.n + a + b)
            lo, hi = wilson(int(r.k), int(r.n))
            rows.append({"GENE": r.GENE, "MUT": r.MUT, "mech_class": cls,
                         "n_carriers": int(r.n), "n_resistant": int(r.k),
                         "raw_prop": float(raw), "shrunk_prop": float(shrunk),
                         "prior_a": a, "prior_b": b, "wilson_lo": lo, "wilson_hi": hi,
                         "grade_flip_0.5": bool((raw >= 0.5) != (shrunk >= 0.5))})
    return pd.DataFrame(rows).sort_values(["mech_class", "n_carriers"], ascending=[True, False]).reset_index(drop=True)