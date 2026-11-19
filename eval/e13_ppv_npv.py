"""E13: PPV/NPV at realistic prevalence (AGENTS.md section 11). Bayes from CV sens/spec."""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import pandas as pd

from discern.prevalence import ppv_npv
from discern.results_io import write_result

# Operating points per spec E13 (rounded CV sens/spec)
MODEL = (0.83, 0.95)
CAT = (0.79, 0.984)
PREV = [0.01, 0.02, 0.05, 0.10]
rows = []
obs = {}
exp = {}
tol = {}
exp_ppv_model = {0.01: 0.14, 0.02: 0.25, 0.05: 0.47, 0.10: 0.65}
exp_npv_model = {0.01: 0.998, 0.02: 0.996, 0.05: 0.991, 0.10: 0.981}
exp_ppv_cat = {0.01: 0.33, 0.02: 0.50, 0.05: 0.72, 0.10: 0.85}
exp_npv_cat = {0.01: 0.998, 0.02: 0.996, 0.05: 0.989, 0.10: 0.977}
for pi in PREV:
    pm, nm = ppv_npv(*MODEL, pi)
    pc, nc = ppv_npv(*CAT, pi)
    rows.append({"prevalence": pi, "model_ppv": pm, "model_npv": nm, "cat_ppv": pc, "cat_npv": nc})
    tag = str(pi)
    obs[f"model_ppv_{tag}"], obs[f"model_npv_{tag}"] = pm, nm
    obs[f"cat_ppv_{tag}"], obs[f"cat_npv_{tag}"] = pc, nc
    exp[f"model_ppv_{tag}"], exp[f"model_npv_{tag}"] = exp_ppv_model[pi], exp_npv_model[pi]
    exp[f"cat_ppv_{tag}"], exp[f"cat_npv_{tag}"] = exp_ppv_cat[pi], exp_npv_cat[pi]
    for k in (f"model_ppv_{tag}", f"cat_ppv_{tag}"):
        tol[k] = 0.02
    for k in (f"model_npv_{tag}", f"cat_npv_{tag}"):
        tol[k] = 0.003
out = write_result("predictive_values", "e13_ppv_npv", obs, expected=exp, tolerance=tol,
                   tables={"ppv_npv": pd.DataFrame(rows)},
                   method="Bayes' rule at prevalence 1/2/5/10% using CV sens/spec: model (0.83/0.95), "
                          "catalogue (0.79/0.984).",
                   interpretation="At 2% prevalence ~half of catalogue R calls are false alarms (PPV 50%). "
                                  "R is a flag for confirmation; S is reliable (NPV ~99.6%).",
                   caveats="Illustrative operating points; real lab prevalence must replace the placeholder.",
                   reproduce="python eval/e13_ppv_npv.py",
                   supports="Prevalence-aware predictive values (clinical framing).")
import json
print("pass_all:", json.loads(out.read_text())["pass_all"])
print(pd.DataFrame(rows).to_string(index=False))