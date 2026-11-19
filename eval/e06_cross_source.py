"""E06: cross-source S3, nested inside training source (AGENTS.md section 11)."""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import json
import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import roc_auc_score

from discern.features import sample_tables
from discern.results_io import write_result
from discern.runner import nested_oof
from discern.splits import cross_source_folds

cfg = yaml.safe_load(open("config.yaml"))
proc = Path(cfg["paths"]["processed"])
s = pd.read_parquet(proc / "samples.parquet")
m = pd.read_parquet(proc / "mutations.parquet")
y, ids, g = s.y.values, s.UNIQUEID.values, s.SUBJID.values
core, a = sample_tables(m, ids, cfg["features"]["frs_min"], tuple(cfg["features"]["dbd_range"]))

EXP = {"CRyPTIC_to_NICD": (0.91, 0.844, 0.963), "NICD_to_CRyPTIC": (0.83, 0.675, 0.953)}
rows, obs, exp, tol = [], {}, {}, {}
for name, tr, te in cross_source_folds(s):
    P, T, th = nested_oof(core, a, y, g, [(tr, te)], C=cfg["model"]["C"])
    p = P[te]
    yt = y[te]
    tb = th[0]["t_bin"]
    pb = p >= tb
    auc = float(roc_auc_score(yt, p))
    sens = float((pb & (yt == 1)).sum() / max((yt == 1).sum(), 1))
    spec = float(((~pb) & (yt == 0)).sum() / max((yt == 0).sum(), 1))
    rows.append({"direction": name, "n_test": len(te), "n_R_test": int(yt.sum()),
                 "auc": auc, "sens": sens, "spec": spec})
    obs[f"{name}_auc"], obs[f"{name}_sens"], obs[f"{name}_spec"] = auc, sens, spec
    exp[f"{name}_auc"], exp[f"{name}_sens"], exp[f"{name}_spec"] = EXP[name]
    for k in ("auc", "sens", "spec"):
        tol[f"{name}_{k}"] = 0.03
    print(f"{name}: AUC {auc:.3f} sens {sens:.3f} spec {spec:.3f}", flush=True)

out = write_result("lr_f2_baseline", "e06_cross_source", obs, expected=exp, tolerance=tol,
                   tables={"cross_source": pd.DataFrame(rows)},
                   method="Train on one SRC, test on the other; nested calibration inside training source; "
                          "binary at per-fold Youden cut.",
                   interpretation="Accuracy depends on lab/source: NICD->CRyPTIC drops (few CRyPTIC resistant "
                                  "training germs, lineage/prevalence shift).",
                   caveats="NICD test has 783 R; CRyPTIC test only 83 R (wide CI).",
                   reproduce="python eval/e06_cross_source.py",
                   supports="Cross-lab generalization limit (honest).")
print("pass_all:", json.loads(out.read_text())["pass_all"])