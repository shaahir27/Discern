"""Step A: model-class ablation on identical F2 features / S1 patient-grouped folds / nested protocol.

Exploratory (no pre-registered acceptance numbers). Hyperparameters fixed a priori (D007), not tuned.
Golden LR-L2 is the reference (E04); candidates: L1 LR, RandomForest, HistGradientBoosting, linear SVM.
"""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import json
import numpy as np
import pandas as pd
import yaml
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.svm import SVC

from discern.features import sample_tables
from discern.model import sens_at_spec
from discern.results_io import write_result
from discern.runner import nested_oof_est

cfg = yaml.safe_load(open("config.yaml"))
proc = Path(cfg["paths"]["processed"])
s = pd.read_parquet(proc / "samples.parquet")
m = pd.read_parquet(proc / "mutations.parquet")
y, ids, g = s.y.values, s.UNIQUEID.values, s.SUBJID.values
folds_arr = np.load(proc / "s1_folds.npy")
folds = [(np.where(folds_arr != k)[0], np.where(folds_arr == k)[0]) for k in range(5)]
core, a = sample_tables(m, ids, cfg["features"]["frs_min"], tuple(cfg["features"]["dbd_range"]))
TH = cfg["thresholds"]

CANDIDATES = {
    "lr_l1": (lambda: LogisticRegression(C=0.3, penalty="l1", solver="liblinear", max_iter=4000),
              "L1 LogisticRegression C=0.3, liblinear"),
    "rf_model": (lambda: RandomForestClassifier(n_estimators=300, min_samples_leaf=5, max_features="sqrt",
                                                n_jobs=-1, random_state=0),
                 "RandomForest 300 trees, min_samples_leaf=5, max_features=sqrt"),
    "hgb_model": (lambda: HistGradientBoostingClassifier(max_iter=200, learning_rate=0.1, max_leaf_nodes=15,
                                                         l2_regularization=1.0, random_state=0),
                  "HistGradientBoosting 200 iters, lr=0.1, 15 leaves, l2=1.0"),
    "svm_model": (lambda: SVC(C=0.3, kernel="linear", probability=False, max_iter=4000),
                  "linear SVM C=0.3 (decision_function -> Platt)"),
}


def metrics_from(P, T, th):
    tb = np.array([th[k]["t_bin"] for k in folds_arr])
    pb = P >= tb
    firm = T != "U"
    nR, nS = max((y == 1).sum(), 1), max((y == 0).sum(), 1)
    return {
        "auc": float(roc_auc_score(y, P)),
        "auprc": float(average_precision_score(y, P)),
        "binary_sens": float((pb & (y == 1)).sum() / nR),
        "binary_spec": float(((~pb) & (y == 0)).sum() / nS),
        "ternary_firm": float(firm.mean()),
        "sens_at_spec_0.95": sens_at_spec(y, P, 0.95),
        "sens_at_spec_0.97": sens_at_spec(y, P, 0.97),
        "sens_at_spec_0.98": sens_at_spec(y, P, 0.98),
        "sens_at_spec_0.985": sens_at_spec(y, P, 0.985),
    }


# Sanity: generalized runner must reproduce golden E04 exactly for LR-L2.
lr_factory = lambda: LogisticRegression(C=0.3, penalty="l2", max_iter=4000)
Plr, Tlr, thlr = nested_oof_est(core, a, y, g, folds, lr_factory,
                                npv=TH["npv_target"], ppv=TH["ppv_target"], min_n=TH["min_rows"])
e04 = json.loads(Path("results/models/lr_model/e04_s1_patient_cv/results.json").read_text())["observed"]
lr_auc = float(roc_auc_score(y, Plr))
assert abs(lr_auc - e04["auc"]) < 1e-6, f"runner fidelity broken: {lr_auc} vs {e04['auc']}"
print(f"runner fidelity OK: LR-L2 AUC {lr_auc:.6f} == E04 {e04['auc']:.6f}", flush=True)

rows = [{"model": "lr_model", "auc": e04["auc"], "auprc": e04["auprc"],
         "binary_sens": e04["binary_sens"], "binary_spec": e04["binary_spec"],
         "ternary_firm": e04["ternary_firm"], "sens_at_spec_0.95": e04["sens_at_spec_0.95"],
         "sens_at_spec_0.97": e04["sens_at_spec_0.97"], "sens_at_spec_0.98": e04["sens_at_spec_0.98"],
         "sens_at_spec_0.985": e04["sens_at_spec_0.985"], "config": "L2 LR C=0.3 (golden E04)"}]

for name, (factory, desc) in CANDIDATES.items():
    P, T, th = nested_oof_est(core, a, y, g, folds, factory,
                              npv=TH["npv_target"], ppv=TH["ppv_target"], min_n=TH["min_rows"])
    obs = metrics_from(P, T, th)
    np.save(proc / f"e05b_{name}_p.npy", P)
    write_result(name, "e05b_model_cv", obs,
                 tables={"metrics": pd.DataFrame([obs]), "thresholds_per_fold": pd.DataFrame(th)},
                 method=f"{desc}. Identical F2 features, S1 patient-grouped 5-fold, nested patient-grouped "
                        f"Platt calibration + thresholds (NPV>={TH['npv_target']}, PPV>={TH['ppv_target']}, Youden).",
                 interpretation=f"Exploratory model-class candidate. AUC {obs['auc']:.4f}, "
                                f"sens@spec0.985 {obs['sens_at_spec_0.985']:.4f} vs golden LR "
                                f"{e04['auc']:.4f}/{e04['sens_at_spec_0.985']:.4f}.",
                 caveats="No pre-registered acceptance numbers; hyperparameters fixed a priori (D007), not tuned. "
                         "Descriptive comparison only.",
                 reproduce="python eval/e05b_model_comparison.py",
                 supports="Model-class ablation: does a non-linear model beat calibrated LR on identical features/splits?",
                 extra={"thresholds_per_fold": th, "config": desc})
    rows.append({"model": name, **obs, "config": desc})
    print(f"{name}: AUC {obs['auc']:.4f} AUPRC {obs['auprc']:.4f} sens@spec0.985 {obs['sens_at_spec_0.985']:.4f}", flush=True)

df = pd.DataFrame(rows).sort_values("auc", ascending=False).reset_index(drop=True)
best = df.iloc[0]
beats_lr = best["model"] != "lr_model" and (best["auc"] - e04["auc"]) > 0.005
decision = (f"Best by AUC: {best['model']} ({best['auc']:.4f}). "
            + ("A candidate beats golden LR by >0.005 AUC — review before shipping."
               if beats_lr else "No candidate beats golden LR by >0.005 AUC; KEEP golden LR (F2) as production."))
obs_cmp = {f"{r['model']}_{k}": r[k] for r in rows for k in
           ("auc", "auprc", "binary_sens", "binary_spec", "sens_at_spec_0.97", "sens_at_spec_0.985")}
write_result("head_to_head", "e05b_comparison", obs_cmp,
             tables={"comparison": df},
             method="Same-data model-class comparison on identical 11,868 samples, F2 features, S1 patient-grouped "
                    "nested CV. Ranked by AUC; clinical operating point sens@spec0.985 also reported.",
             interpretation=decision,
             caveats="Exploratory; hyperparameters fixed a priori (D007), not tuned. match != beat. "
                     "AUC differences <0.005 are within noise.",
             reproduce="python eval/e05b_model_comparison.py",
             supports="Justifies the model class choice (calibrated LR) with a same-data ablation.")
print("\n=== COMPARISON (ranked by AUC) ===")
print(df[["model", "auc", "auprc", "binary_sens", "binary_spec", "sens_at_spec_0.985"]].to_string(index=False))
print("\nDECISION:", decision)