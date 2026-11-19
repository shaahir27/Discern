"""E05: feature-set ablation on S1 (AGENTS.md section 11)."""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import json
import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

from discern.features import sample_tables
from discern.results_io import write_result
from discern.runner import nested_oof

cfg = yaml.safe_load(open("config.yaml"))
proc = Path(cfg["paths"]["processed"])
s = pd.read_parquet(proc / "samples.parquet")
m = pd.read_parquet(proc / "mutations.parquet")
y, ids, g = s.y.values, s.UNIQUEID.values, s.SUBJID.values
folds = list(GroupKFold(n_splits=5).split(s, groups=s.SUBJID))

EXP = {("F0", 0.1): (0.791, 0.503), ("F1", 0.1): (0.905, 0.672),
       ("F2", 0.1): (0.909, 0.684), ("F2", 0.75): (0.799, 0.548)}
rows = []
for (fset, frs), (ea, ep) in EXP.items():
    core, a = sample_tables(m, ids, frs, tuple(cfg["features"]["dbd_range"]))
    P, _, _ = nested_oof(core, a, y, g, folds, fset=fset, C=cfg["model"]["C"])
    auc, auprc = float(roc_auc_score(y, P)), float(average_precision_score(y, P))
    rows.append({"fset": fset, "frs_min": frs, "auc": auc, "auprc": auprc,
                 "exp_auc": ea, "exp_auprc": ep, "auc_within_0.015": abs(auc - ea) <= 0.015,
                 "auprc_within_0.015": abs(auprc - ep) <= 0.015})
    print(f"{fset} frs{frs}: AUC {auc:.3f} AUPRC {auprc:.3f}", flush=True)

df = pd.DataFrame(rows)
obs = {f"{r['fset']}_frs{r['frs_min']}_auc": r["auc"] for r in rows}
obs.update({f"{r['fset']}_frs{r['frs_min']}_auprc": r["auprc"] for r in rows})
exp = {f"{k[0]}_frs{k[1]}_auc": v[0] for k, v in EXP.items()}
exp.update({f"{k[0]}_frs{k[1]}_auprc": v[1] for k, v in EXP.items()})
tol = {k: 0.015 for k in exp}
out = write_result("lr_f2_baseline", "e05_ablation", obs, expected=exp, tolerance=tol,
                   tables={"ablation": df},
                   method="S1 nested CV; F0 vocab only, F1 +type/region, F2 +mmpL5/interaction; "
                          "F2 frs_min 0.75 ignores minor alleles.",
                   interpretation="Type/region features ~+0.11 AUC; minor alleles ~+0.11 AUC; "
                                  "mmpL5 interaction ~+0.004 AUC.",
                   caveats="Ablation, not tuning. frs 0.75 drops minor-allele signal.",
                   reproduce="python eval/e05_ablation.py",
                   supports="Quantifies feature value (minor alleles published by Fowler lab).")
print("pass_all:", json.loads(out.read_text())["pass_all"])