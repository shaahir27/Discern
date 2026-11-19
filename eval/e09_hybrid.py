"""E09: hybrid on catalogue-U samples (AGENTS.md section 11). Exploratory; never auto-call."""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import json
import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import roc_auc_score

from discern.results_io import write_result

cfg = yaml.safe_load(open("config.yaml"))
proc = Path(cfg["paths"]["processed"])
s = pd.read_parquet(proc / "samples.parquet")
cat = pd.read_parquet(proc / "e03_oof_calls.parquet")
oof = pd.read_parquet(proc / "e04_oof.parquet")
th = json.loads((Path("results/models/lr_model/e04_s1_patient_cv/results.json")).read_text())["extra"]["thresholds_per_fold"]

df = s[["UNIQUEID", "y"]].merge(cat[["UNIQUEID", "cat_call", "fold"]], on="UNIQUEID").merge(
    oof[["UNIQUEID", "p"]], on="UNIQUEID")
U = df[df.cat_call == "U"].copy()
auc_in_U = float(roc_auc_score(U.y, U.p))
lo = np.array([th[int(f)]["lo"] for f in U.fold])
hi = np.array([th[int(f)]["hi"] for f in U.fold])
resolve_S = U.p <= lo
resolve_R = U.p >= hi
remain = ~(resolve_S | resolve_R)
obs = {
    "n_U": len(U), "n_U_resistant": int(U.y.sum()), "auc_ml_within_U": auc_in_U,
    "resolved_S": int(resolve_S.sum()), "resolved_S_resistant": int(U.y[resolve_S].sum()),
    "resolved_R": int(resolve_R.sum()), "resolved_R_susceptible": int((1 - U.y)[resolve_R].sum()),
    "remain_U": int(remain.sum()), "remain_U_resistant": int(U.y[remain].sum()),
}
exp = {"auc_ml_within_U": 0.744, "resolved_S": 173, "resolved_S_resistant": 17,
       "resolved_R": 126, "resolved_R_susceptible": 51, "remain_U": 362, "remain_U_resistant": 114}
tol = {"auc_ml_within_U": 0.03, "resolved_S": 15, "resolved_S_resistant": 15,
       "resolved_R": 15, "resolved_R_susceptible": 15, "remain_U": 15, "remain_U_resistant": 15}
out = write_result("hybrid_exploratory", "e09_hybrid", obs, expected=exp, tolerance=tol,
                   tables={"U_samples": U},
                   method="Among catalogue-U samples, apply per-fold model lo/hi to resolve some. "
                          "ML AUC within U measures ranking value.",
                   interpretation="EXPLORATORY: model ranks U samples (AUC ~0.74) but resolving them is unsafe "
                                  "(resolved-S group ~9.8% resistant; resolved-R ~40% susceptible). NEVER auto-call.",
                   caveats="Used only to RANK the unresolved queue, not to change calls.",
                   reproduce="python eval/e09_hybrid.py",
                   supports="Queue ranking value; rejects unsafe auto-resolution (honest).")
print(json.dumps({"observed": obs, "pass_all": json.loads(out.read_text())["pass_all"]}, indent=1))