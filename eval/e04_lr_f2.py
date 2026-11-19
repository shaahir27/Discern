"""E04: Tier B calibrated logistic regression (F2) on S1, nested (AGENTS.md sections 8-9)."""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import json
import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import average_precision_score, roc_auc_score

from discern.features import design, sample_tables, vocab_from
from discern.model import fit_nested, score, sens_at_spec, thresholds
from discern.results_io import write_result

cfg = yaml.safe_load(open("config.yaml"))
proc = Path(cfg["paths"]["processed"])
s = pd.read_parquet(proc / "samples.parquet")
m = pd.read_parquet(proc / "mutations.parquet")
y = s.y.values
ids = s.UNIQUEID.values
g = s.SUBJID.values
folds = np.load(proc / "s1_folds.npy")

core, a = sample_tables(m, ids, cfg["features"]["frs_min"], tuple(cfg["features"]["dbd_range"]))
P = np.zeros(len(s))
R = np.zeros(len(s))
T = np.empty(len(s), dtype=object)
th = []
for k in range(5):
    te = folds == k
    tr = ~te
    vocab = vocab_from(core, tr, cfg["features"]["vocab_min_count"])
    X = design(core, a, vocab, "F2").values
    final, platt, poof = fit_nested(X[tr], y[tr], g[tr], cfg["model"]["C"], cfg["model"]["max_iter"])
    t = thresholds(poof, y[tr], cfg["thresholds"]["npv_target"], cfg["thresholds"]["ppv_target"], cfg["thresholds"]["min_rows"])
    th.append(t)
    raw, p = score(final, platt, X[te])
    P[te] = p
    R[te] = raw
    T[te] = np.where(p <= t["lo"], "S", np.where(p >= t["hi"], "R", "U"))
    print(f"fold {k}: n_vocab={len(vocab)} n_feat={X.shape[1]} lo={t['lo']:.4f} hi={t['hi']:.4f} t_bin={t['t_bin']:.4f}", flush=True)

tb = np.array([th[k]["t_bin"] for k in folds])
pb = P >= tb
firm = T != "U"
nR = max((y == 1).sum(), 1)
nS = max((y == 0).sum(), 1)
obs = {
    "auc": float(roc_auc_score(y, P)),
    "auprc": float(average_precision_score(y, P)),
    "binary_sens": float((pb & (y == 1)).sum() / nR),
    "binary_spec": float(((~pb) & (y == 0)).sum() / nS),
    "ternary_firm": float(firm.mean()),
    "ternary_sens_among_firm": float(((T == "R") & (y == 1)).sum() / max(((y == 1) & firm).sum(), 1)),
    "ternary_spec_among_firm": float(((T == "S") & (y == 0)).sum() / max(((y == 0) & firm).sum(), 1)),
    "sens_at_spec_0.95": sens_at_spec(y, P, 0.95),
    "sens_at_spec_0.97": sens_at_spec(y, P, 0.97),
    "sens_at_spec_0.98": sens_at_spec(y, P, 0.98),
    "sens_at_spec_0.985": sens_at_spec(y, P, 0.985),
    "n_features_all_data": int(design(core, a, vocab_from(core, np.ones(len(s), bool)), "F2").shape[1]),
}
np.save(proc / "e04_oof_p.npy", P)
np.save(proc / "e04_oof_raw.npy", R)
pd.DataFrame({"UNIQUEID": ids, "fold": folds, "p": P, "raw": R, "ternary": T}).to_parquet(proc / "e04_oof.parquet")

EXP = {"auc": 0.909, "auprc": 0.684, "binary_sens": 0.835, "binary_spec": 0.948, "ternary_firm": 0.947,
       "sens_at_spec_0.95": 0.830, "sens_at_spec_0.97": 0.782, "sens_at_spec_0.98": 0.652, "sens_at_spec_0.985": 0.604}
TOL = {"auc": 0.01, "auprc": 0.02, "binary_sens": 0.02, "binary_spec": 0.01, "ternary_firm": 0.01,
       "sens_at_spec_0.95": 0.03, "sens_at_spec_0.97": 0.03, "sens_at_spec_0.98": 0.03, "sens_at_spec_0.985": 0.03}
out = write_result(
    "lr_f2_baseline", "e04_s1_patient_cv", obs,
    expected=EXP, tolerance=TOL,
    tables={"metrics": pd.DataFrame([obs]), "thresholds_per_fold": pd.DataFrame(th)},
    method="L2 LogisticRegression C=0.3; nested patient-grouped Platt calibration; thresholds from inner OOF "
           "(NPV>=0.985, PPV>=0.80, Youden binary); F2 features; vocabulary>=3 training carriers.",
    interpretation="Calibrated LR matches the catalogue at matched operating points; it does not beat it. "
                   "Adds ranking + calibrated probability.",
    caveats="ternary sens/spec among firm (0.735/0.989) are descriptive (no tolerance). Model != better accuracy.",
    reproduce="python eval/e03_catalogue_cv.py && python eval/e04_lr_f2.py",
    supports="Baseline reproduction; model adds calibration/ranking, not superior accuracy.",
    extra={"thresholds_per_fold": th},
)
rec = json.loads(out.read_text())
print(json.dumps({"observed": {k: v for k, v in obs.items() if k != "n_features_all_data"},
                  "n_features_all_data": obs["n_features_all_data"], "pass_all": rec["pass_all"]}, indent=1))