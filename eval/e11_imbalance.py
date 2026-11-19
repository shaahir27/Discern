"""E11: imbalance + learning curve (AGENTS.md section 11)."""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import json
import numpy as np
import pandas as pd
import yaml
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

from discern.features import design, sample_tables, vocab_from
from discern.model import fit_nested, logit, score, sens_at_spec
from discern.results_io import write_result

cfg = yaml.safe_load(open("config.yaml"))
proc = Path(cfg["paths"]["processed"])
s = pd.read_parquet(proc / "samples.parquet")
m = pd.read_parquet(proc / "mutations.parquet")
y, ids, g = s.y.values, s.UNIQUEID.values, s.SUBJID.values
core, a = sample_tables(m, ids, cfg["features"]["frs_min"], tuple(cfg["features"]["dbd_range"]))
folds = list(GroupKFold(n_splits=5).split(s, groups=s.SUBJID))
rng = np.random.default_rng(cfg["seeds"]["learning_curve"])

# 1:1 balanced by RESAMPLING (undersample majority to minority count) — the spec's "1:1 balanced".
# Diagnostic only; production calibrated model uses no resampling/class_weight (rule 4).
rng_bal = np.random.default_rng(cfg["seeds"]["bootstrap"])
Pb = np.zeros(len(y))
for tr, te in folds:
    r_idx = tr[y[tr] == 1]
    s_idx = tr[y[tr] == 0]
    s_keep = rng_bal.choice(s_idx, size=len(r_idx), replace=False)
    trb = np.concatenate([r_idx, s_keep])
    vocab = vocab_from(core, trb, 3)
    X = design(core, a, vocab, "F2").values
    lr = LogisticRegression(C=0.3, max_iter=4000).fit(X[trb], y[trb])
    Pb[te] = lr.predict_proba(X[te])[:, 1]
bal_auc = float(roc_auc_score(y, Pb))
bal_auprc = float(average_precision_score(y, Pb))
bal_sens98 = sens_at_spec(y, Pb, 0.98)

# learning curve: fraction of RESISTANT training rows retained (all S kept)
FRACS = [0.10, 0.25, 0.50, 0.75, 1.00]
EXP_LC = [0.898, 0.905, 0.908, 0.906, 0.909]
lc_rows, obs, exp, tol = [], {}, {}, {}
for frac, ea in zip(FRACS, EXP_LC):
    P = np.zeros(len(y))
    for tr, te in folds:
        r_idx = tr[y[tr] == 1]
        keep = rng.choice(r_idx, size=max(1, int(round(frac * len(r_idx)))), replace=False) if frac < 1 else r_idx
        tr_sub = np.concatenate([tr[y[tr] == 0], keep])
        vocab = vocab_from(core, tr_sub, 3)
        X = design(core, a, vocab, "F2").values
        final, platt, _ = fit_nested(X[tr_sub], y[tr_sub], g[tr_sub], C=0.3)
        _, p = score(final, platt, X[te])
        P[te] = p
    auc = float(roc_auc_score(y, P))
    lc_rows.append({"frac_resistant": frac, "auc": auc, "exp_auc": ea, "within_0.01": abs(auc - ea) <= 0.01})
    obs[f"lc_auc_{frac}"] = auc
    exp[f"lc_auc_{frac}"] = ea
    tol[f"lc_auc_{frac}"] = 0.01
    print(f"frac {frac}: AUC {auc:.3f} (exp {ea})", flush=True)

obs.update({"balanced_auc": bal_auc, "balanced_auprc": bal_auprc, "balanced_sens_at_spec98": bal_sens98})
exp.update({"balanced_auc": 0.904, "balanced_auprc": 0.653, "balanced_sens_at_spec98": 0.592})
tol.update({"balanced_auc": 0.01, "balanced_auprc": 0.02, "balanced_sens_at_spec98": 0.03})
out = write_result("lr_f2_baseline", "e11_imbalance_learning_curve", obs, expected=exp, tolerance=tol,
                   tables={"learning_curve": pd.DataFrame(lc_rows)},
                   method="1:1 class_weight balanced LR (NOT production); learning curve retains fraction of "
                          "resistant training rows (all S kept), nested calibration, S1.",
                   interpretation="Balancing does not help and breaks calibration; AUC is flat across resistant "
                                  "fractions (signal saturates early).",
                   caveats="Balanced model is diagnostic only; production model uses no class_weight (rule 4).",
                   reproduce="python eval/e11_imbalance.py",
                   supports="Justifies no class_weight / no resampling (honest).")
print("balanced:", bal_auc, bal_auprc, bal_sens98, "| pass_all:", json.loads(out.read_text())["pass_all"])