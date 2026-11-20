"""G4 / P-E15: value-of-information lab queue + retrospective label-reveal simulation.

Ranks catalogue-U variant sets by info-gain (n_U x beta-binomial entropy reduction) vs count-only
vs random; reveals top-K sets' existing labels, unlocks variants (binomial proxy), measures %U
resolved. Priors use ONLY non-U (firm) evidence — no leakage of U labels into ranking.
"""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import json
import numpy as np
import pandas as pd
import yaml

from discern.catalogue import cat_muts
from discern.results_io import write_result
from discern.voi import expected_info_gain, rank_sets, variant_unlock

cfg = yaml.safe_load(open("config.yaml"))
proc = Path(cfg["paths"]["processed"])
s = pd.read_parquet(proc / "samples.parquet")
m = pd.read_parquet(proc / "mutations.parquet")
cat = pd.read_parquet(proc / "e03_oof_calls.parquet")
y_by = s.set_index("UNIQUEID").y
U_ids = list(cat[cat.cat_call == "U"].UNIQUEID)
firm_ids = set(cat[cat.cat_call != "U"].UNIQUEID)
total_U = len(U_ids)

# non-U (firm) evidence per variant -> ranking prior (no U-label leakage)
mf = cat_muts(m[m.UNIQUEID.isin(firm_ids) & (m.FRSf >= 0.1)]).drop_duplicates(["UNIQUEID", "MUT"])
mf["y"] = mf.UNIQUEID.map(y_by)
fe = mf.groupby("MUT").agg(nf=("UNIQUEID", "nunique"), kf=("y", "sum"))
nf, kf = fe.nf.to_dict(), fe.kf.to_dict()

# U samples -> variant set + variant list
mu = cat_muts(m[m.UNIQUEID.isin(U_ids) & (m.FRSf >= 0.1)]).drop_duplicates(["UNIQUEID", "MUT"])
sample_vars = mu.groupby("UNIQUEID").MUT.apply(list).to_dict()
key_by = mu.groupby("UNIQUEID").MUT.apply(lambda x: "|".join(sorted(x))).to_dict()
uset = pd.DataFrame({"UNIQUEID": U_ids})
uset["key"] = uset.UNIQUEID.map(lambda i: key_by.get(i, ""))
uset["y"] = uset.UNIQUEID.map(y_by)

set_rows = []
for key, g in uset.groupby("key"):
    variants = key.split("|") if key else []
    a = 1 + sum(kf.get(v, 0) for v in variants)
    b = 9 + sum(nf.get(v, 0) - kf.get(v, 0) for v in variants)
    set_rows.append({"key": key, "n": len(g), "nR": int(g.y.sum()),
                     "members": list(g.UNIQUEID), "variants": variants,
                     "eig": expected_info_gain(a, b), "prior_a": a, "prior_b": b})
set_df = pd.DataFrame(set_rows).set_index("key")
members = set_df.members.to_dict()
KS = [1, 2, 5, 10, 20, 50]


def simulate(order, K):
    revealed = set(order[:K])
    n_rev, k_rev = {}, {}
    revealed_samples = set()
    for key in revealed:
        for uid in members[key]:
            revealed_samples.add(uid)
            yy = int(y_by[uid])
            for v in sample_vars.get(uid, []):
                n_rev[v] = n_rev.get(v, 0) + 1
                k_rev[v] = k_rev.get(v, 0) + yy
    firm_vars = set()
    for v in set(n_rev) | set(nf):
        nn = nf.get(v, 0) + n_rev.get(v, 0)
        kk = kf.get(v, 0) + k_rev.get(v, 0)
        if variant_unlock(nn, kk) != "U":
            firm_vars.add(v)
    resolved = 0
    for uid in U_ids:
        if key_by.get(uid, "") in revealed or any(v in firm_vars for v in sample_vars.get(uid, [])):
            resolved += 1
    return resolved / total_U, resolved, len(revealed_samples)


rows, obs = [], {}
for strategy in ("info_gain", "count", "random"):
    order = rank_sets(set_df, strategy, seed=cfg["seeds"]["bootstrap"])
    for K in KS:
        pct, nres, nrev = simulate(order, K)
        rows.append({"strategy": strategy, "K": K, "pct_U_resolved": pct,
                     "n_resolved": nres, "n_samples_revealed": nrev,
                     "n_unrevealed_counted_resolved": nres - nrev})
        obs[f"{strategy}_pct_K{K}"] = pct
curves = pd.DataFrame(rows)
obs.update({"total_U": total_U, "n_distinct_sets": int(len(set_df)),
            "info_gain_minus_count_K10": obs["info_gain_pct_K10"] - obs["count_pct_K10"],
            "info_gain_minus_count_K50": obs["info_gain_pct_K50"] - obs["count_pct_K50"]})
top = set_df.assign(score=set_df.n * set_df.eig).sort_values("score", ascending=False).head(20)
top = top.reset_index()[["key", "n", "nR", "eig", "score", "prior_a", "prior_b"]]

write_result("voi_queue", "e15_label_reveal_sim", obs,
             tables={"resolution_curves": curves, "info_gain_ranking": top},
             method="Retrospective comparison of variant sets. K counts selected sets. Each selected "
                    "set reveals every member's saved lab label. The simulated total resolved includes "
                    "those directly revealed samples and unrevealed samples counted by the proxy.",
             interpretation="The strategies compare equal numbers of variant sets. Those sets contain "
                            "different numbers of samples.",
             caveats="The comparison uses unequal numbers of sample labels at the same K. "
                     "Unrevealed samples counted as resolved include cases the proxy can already "
                     "count using firm-case evidence alone. These totals cannot establish gains "
                     "at an equal lab testing budget or production catalogue performance.",
             reproduce="python eval/e15_voi.py",
             supports="Descriptive archive simulation. The resolution proxy does not amend the catalogue.")
print(curves.pivot(index="K", columns="strategy", values="pct_U_resolved").to_string())
print("info_gain - count: K10 %+.3f  K50 %+.3f" % (obs["info_gain_minus_count_K10"], obs["info_gain_minus_count_K50"]))
