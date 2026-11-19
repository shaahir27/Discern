"""E12: variant priority list + evidence cards among catalogue-U samples (AGENTS.md section 11)."""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import json
import math
import numpy as np
import pandas as pd
import yaml

from discern.catalogue import cat_muts
from discern.results_io import write_result

cfg = yaml.safe_load(open("config.yaml"))
proc = Path(cfg["paths"]["processed"])
s = pd.read_parquet(proc / "samples.parquet")
m = pd.read_parquet(proc / "mutations.parquet")
cat = pd.read_parquet(proc / "e03_oof_calls.parquet")
U_ids = set(cat[cat.cat_call == "U"].UNIQUEID)
y_by = s.set_index("UNIQUEID").y
src_by = s.set_index("UNIQUEID").SRC


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


# variant sets among U samples (core variants used by catalogue)
mu = cat_muts(m[m.UNIQUEID.isin(U_ids) & (m.FRSf >= 0.1)])
keys = mu.groupby("UNIQUEID").MUT.apply(lambda x: "|".join(sorted(x)))
uset = pd.DataFrame({"UNIQUEID": list(U_ids)})
uset["key"] = uset.UNIQUEID.map(keys).fillna("")
uset["y"] = uset.UNIQUEID.map(y_by)
grp = uset.groupby("key").agg(n=("UNIQUEID", "size"), nR=("y", "sum")).sort_values("n", ascending=False)
total_U = len(uset)
cov = {k: float(grp.n.head(k).sum() / total_U) for k in (5, 10, 20, 50)}
top_sets = grp.head(20).reset_index()
top_sets["pct_resistant"] = top_sets.nR / top_sets.n

# per-variant evidence cards (pooled across ALL samples) + by source
mall = cat_muts(m[m.FRSf >= 0.1]).copy()
mall["y"] = mall.UNIQUEID.map(y_by)
mall["SRC"] = mall.UNIQUEID.map(src_by)
ev_rows = []
for mut, d in mall.groupby("MUT"):
    n = d.UNIQUEID.nunique()
    nr = int(d[d.y == 1].UNIQUEID.nunique())
    lo, hi = wilson(nr, n)
    row = {"variant": mut, "n_carriers": n, "n_resistant": nr, "pct_resistant": nr / n if n else 0,
           "ci95_low": lo, "ci95_high": hi}
    for src in ("CRyPTIC", "NICD"):
        ds = d[d.SRC == src]
        row[f"{src}_n"] = int(ds.UNIQUEID.nunique())
        row[f"{src}_r"] = int(ds[ds.y == 1].UNIQUEID.nunique())
    ev_rows.append(row)
ev = pd.DataFrame(ev_rows).sort_values("n_carriers", ascending=False)

obs = {"n_U": total_U, "n_distinct_sets": int(len(grp)),
       "top5_coverage": cov[5], "top10_coverage": cov[10], "top20_coverage": cov[20], "top50_coverage": cov[50],
       "largest_set": str(grp.index[0]), "largest_set_n": int(grp.n.iloc[0]),
       "largest_set_pctR": float(grp.nR.iloc[0] / grp.n.iloc[0])}
exp = {"n_distinct_sets": 332, "top5_coverage": 0.28, "top10_coverage": 0.34,
       "top20_coverage": 0.40, "top50_coverage": 0.50}
tol = {"n_distinct_sets": 15, "top5_coverage": 0.03, "top10_coverage": 0.03,
       "top20_coverage": 0.03, "top50_coverage": 0.03}
out = write_result("priority_list", "e12_variant_priority", obs, expected=exp, tolerance=tol,
                   tables={"top_sets": top_sets, "evidence_cards": ev},
                   method="Among catalogue-U samples, group by sorted core-variant set; rank by carrier count; "
                          "top-K coverage. Per-variant Wilson 95% CI + by-source counts.",
                   interpretation="332 distinct sets cause 661 referrals; top-10 sets cover ~34%. Phenotyping a "
                                  "handful of variants would cut referrals sharply (descriptive).",
                   caveats="NICD enriches resistance; evidence cards split by source. Descriptive, not causal.",
                   reproduce="python eval/e12_priority.py",
                   supports="Referral concentration (basis for the VOI queue novelty).")
print(json.dumps({"observed": obs, "pass_all": json.loads(out.read_text())["pass_all"]}, indent=1))
print("largest sets:\n", top_sets.head(5).to_string(index=False))