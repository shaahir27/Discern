"""E10: prevalence-shift correction across sources (AGENTS.md section 11)."""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import json
import numpy as np
import pandas as pd
import yaml

from discern.features import sample_tables
from discern.prevalence import adjust, brier, ece
from discern.results_io import write_result
from discern.runner import nested_oof
from discern.splits import cross_source_folds

cfg = yaml.safe_load(open("config.yaml"))
proc = Path(cfg["paths"]["processed"])
s = pd.read_parquet(proc / "samples.parquet")
m = pd.read_parquet(proc / "mutations.parquet")
y, ids, g = s.y.values, s.UNIQUEID.values, s.SUBJID.values
core, a = sample_tables(m, ids, cfg["features"]["frs_min"], tuple(cfg["features"]["dbd_range"]))
pi = {src: float(s[s.SRC == src].y.mean()) for src in ("CRyPTIC", "NICD")}

EXP = {"NICD_to_CRyPTIC": (0.288, 0.013, 0.102, 0.017), "CRyPTIC_to_NICD": (0.568, 0.101, 0.501, 0.124)}
TOL = {"NICD_to_CRyPTIC": (0.01, 0.01, 0.005, 0.005), "CRyPTIC_to_NICD": (0.03, 0.03, 0.02, 0.02)}
rows, obs, exp, tol = [], {}, {}, {}
for name, tr, te in cross_source_folds(s):
    src_tr, src_te = name.split("_to_")
    P, _, _ = nested_oof(core, a, y, g, [(tr, te)], C=cfg["model"]["C"])
    p = P[te]
    yt = y[te]
    p_adj = adjust(p, pi[src_tr], pi[src_te])
    ece_raw, ece_adj = ece(yt, p), ece(yt, p_adj)
    br_raw, br_adj = brier(yt, p), brier(yt, p_adj)
    rows.append({"direction": name, "pi_train": pi[src_tr], "pi_test": pi[src_te],
                 "ece_raw": ece_raw, "ece_adj": ece_adj, "brier_raw": br_raw, "brier_adj": br_adj})
    obs[f"{name}_ece_raw"], obs[f"{name}_ece_adj"] = ece_raw, ece_adj
    obs[f"{name}_brier_raw"], obs[f"{name}_brier_adj"] = br_raw, br_adj
    exp[f"{name}_ece_raw"], exp[f"{name}_ece_adj"] = EXP[name][0], EXP[name][1]
    exp[f"{name}_brier_raw"], exp[f"{name}_brier_adj"] = EXP[name][2], EXP[name][3]
    tol[f"{name}_ece_raw"], tol[f"{name}_ece_adj"] = TOL[name][0], TOL[name][1]
    tol[f"{name}_brier_raw"], tol[f"{name}_brier_adj"] = TOL[name][2], TOL[name][3]
    print(f"{name}: ECE {ece_raw:.3f}->{ece_adj:.3f} Brier {br_raw:.3f}->{br_adj:.3f}", flush=True)

out = write_result("prevalence_correction", "e10_prevalence_shift", obs, expected=exp, tolerance=tol,
                   tables={"prevalence_shift": pd.DataFrame(rows)},
                   method="Cross-source OOF probabilities; odds-ratio prevalence correction (section 10); "
                          "ECE (10 quantile bins) + Brier raw vs adjusted.",
                   interpretation="Calibrated probabilities do NOT transfer across prevalence unless corrected. "
                                  "NICD->CRyPTIC ECE 0.288->0.013; reverse improves less (only 83 CRyPTIC R).",
                   caveats="At low prevalence the PPV-80% R threshold is rarely reached — correct behaviour.",
                   reproduce="python eval/e10_prevalence_shift.py",
                   supports="Prevalence correction validated across labs (math is Saerens 2002; application here).")
print("pass_all:", json.loads(out.read_text())["pass_all"])