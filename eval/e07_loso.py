"""E07: leave-one-site-out S2 (AGENTS.md section 11). No pooled AUC reported."""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import json
import numpy as np
import pandas as pd
import yaml

from discern.features import sample_tables
from discern.results_io import write_result
from discern.runner import nested_oof
from discern.splits import loso_folds

cfg = yaml.safe_load(open("config.yaml"))
proc = Path(cfg["paths"]["processed"])
s = pd.read_parquet(proc / "samples.parquet")
m = pd.read_parquet(proc / "mutations.parquet")
y, ids, g = s.y.values, s.UNIQUEID.values, s.SUBJID.values
core, a = sample_tables(m, ids, cfg["features"]["frs_min"], tuple(cfg["features"]["dbd_range"]))

folds = [(tr, te) for _, tr, te in loso_folds(s)]
P, T, th = nested_oof(core, a, y, g, folds, C=cfg["model"]["C"])
# binary prediction using each site's fold Youden cut
pb = np.zeros(len(s), dtype=bool)
for k, (site, tr, te) in enumerate(loso_folds(s)):
    pb[te] = P[te] >= th[k]["t_bin"]

sens = float((pb & (y == 1)).sum() / max((y == 1).sum(), 1))
spec = float(((~pb) & (y == 0)).sum() / max((y == 0).sum(), 1))
site_rows = []
for site in sorted(s.SITEID.unique()):
    sel = s.SITEID.values == site
    n, nR = int(sel.sum()), int((sel & (y == 1)).sum())
    ss = float(((~pb) & (y == 0) & sel).sum() / max(((y == 0) & sel).sum(), 1))
    se = float((pb & (y == 1) & sel).sum() / nR) if nR >= 25 else None
    site_rows.append({"site": int(site), "n": n, "n_R": nR, "spec": ss,
                      "sens_if_nR_ge25": se, "spec_reported": n >= 300})
df = pd.DataFrame(site_rows)
big = df[df.spec_reported]
obs = {"pooled_sens": sens, "pooled_spec": spec,
       "per_site_spec_min_n300": float(big.spec.min()), "per_site_spec_max_n300": float(big.spec.max())}
exp = {"pooled_sens": 0.788, "pooled_spec": 0.954}
tol = {"pooled_sens": 0.02, "pooled_spec": 0.01}
out = write_result("lr_f2_baseline", "e07_loso", obs, expected=exp, tolerance=tol,
                   tables={"per_site": df},
                   method="Leave-one-site-out (23 sites); nested calibration per training set; "
                          "pooled binary sens/spec; per-site spec for n>=300; per-site sens only if nR>=25.",
                   interpretation="Pooled sens/spec across held-out sites. Per-site specificity stable; "
                                  "per-site sensitivity unreliable (most sites <25 R).",
                   caveats="Do NOT report pooled LOSO AUC (mixes 23 calibrations). Site 5 (inferred India, "
                           "unverified) has ~5 R -> false-alarm audit only.",
                   reproduce="python eval/e07_loso.py",
                   supports="Cross-site generalization (honest).")
print(json.dumps({"observed": obs, "pass_all": json.loads(out.read_text())["pass_all"]}, indent=1))