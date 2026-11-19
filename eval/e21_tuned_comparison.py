"""E21 retrospective model-class comparison on fixed S1 and KM cohorts (D012–D013)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from joblib import parallel_config
from sklearn.base import BaseEstimator, clone
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.svm import SVC

sys.path.insert(0, "src")

from discern.features import CORE, design, sample_tables, vocab_from
from discern.grammar import classify
from discern.model import sens_at_spec, thresholds
from discern.results_io import write_result
from discern.runner import _raw_score, score_est

LR_GRID = {"C": [0.03, 0.1, 0.3, 1.0, 3.0]}
CANDIDATES = {
    "lr_l2_tuned": (LogisticRegression(penalty="l2", max_iter=4000), LR_GRID),
    "lr_l1_tuned": (
        LogisticRegression(penalty="l1", solver="liblinear", max_iter=4000, random_state=1),
        LR_GRID,
    ),
    "rf_tuned": (
        RandomForestClassifier(random_state=1, n_jobs=-1),
        {"n_estimators": [300, 600], "min_samples_leaf": [1, 3, 5, 10], "max_features": ["sqrt", 0.5]},
    ),
    "hgb_tuned": (
        HistGradientBoostingClassifier(random_state=1),
        {
            "learning_rate": [0.05, 0.1], "max_leaf_nodes": [15, 31],
            "l2_regularization": [0, 1, 5], "max_iter": [200, 400],
        },
    ),
    "svm_tuned": (
        SVC(probability=False, max_iter=4000, gamma="scale"),
        {"C": [0.1, 0.3, 1, 3], "kernel": ["linear", "rbf"]},
    ),
}
METHOD = (
    "F2 with a vocabulary of variants carried by at least three training samples. "
    "Select by ROC AUC on three patient-grouped inner folds, then fit Platt scaling on raw "
    "out-of-fold scores from those folds and derive the golden thresholds from the calibrated scores."
)
CAVEATS = (
    "Research use only. The selected settings and Platt calibration use the same inner "
    "folds. The outer and KM results are held-out evaluations."
)
KM_CAVEATS = (
    "KM is a sensitivity analysis with undocumented assay and critical concentration, "
    "medium-quality phenotypes, and few resistant isolates. Sensitivity at fixed specificity "
    "describes the held-out ROC curve."
)
REPRODUCE = "C:/Users/shaah/miniconda3/envs/discern/python.exe eval/e21_tuned_comparison.py"


def fit_tuned(
    estimator: BaseEstimator,
    grid: dict,
    X: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray,
) -> tuple[BaseEstimator, LogisticRegression, dict, dict]:
    inner = GroupKFold(3)
    search = GridSearchCV(
        estimator, grid, scoring="roc_auc", cv=inner, n_jobs=-1,
        refit=False, error_score="raise",
    )
    with parallel_config(backend="loky", inner_max_num_threads=1):
        search.fit(X, y, groups=groups)
    params = search.best_params_
    best = clone(estimator).set_params(**params)

    raw_oof = np.empty((len(y), 1))
    for tr, te in inner.split(X, y, groups):
        est = clone(best).fit(X[tr], y[tr])
        raw_oof[te, 0] = _raw_score(est, X[te])
    platt = LogisticRegression(C=1e6).fit(raw_oof, y)
    poof = platt.predict_proba(raw_oof)[:, 1]
    limits = thresholds(poof, y, npv=0.985, ppv=0.80, min_n=20)
    return best.fit(X, y), platt, limits, params


def km_tables(train_ids: np.ndarray, features: dict) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    ext = Path("data/external_cryptic_v34")
    dst = pd.read_parquet(
        ext / "DST_MEASUREMENTS.parquet", columns=["SOURCE", "DRUG", "PHENOTYPE"],
        filters=[("SOURCE", "==", "KMVALIDATION2025"), ("DRUG", "==", "BDQ")],
    ).reset_index()
    km = dst.drop_duplicates("UNIQUEID")
    ids, y = km.UNIQUEID.to_numpy(), km.PHENOTYPE.to_numpy() == "R"
    assert len(ids) == 792 and y.sum() == 35 and km.PHENOTYPE.isin(["S", "R"]).all()
    assert not np.intersect1d(ids, train_ids).size, "KM overlaps golden training"

    raw = pd.read_parquet(
        ext / "MUTATIONS.parquet",
        columns=["GENE", "IS_MINOR", "FRS", "MINOR_MUTATION", "MUTATION"],
        filters=[("GENE", "in", [*CORE, "mmpL5", "mmpS5"])],
    ).reset_index()
    raw = raw.loc[raw.UNIQUEID.isin(ids)].reset_index(drop=True)
    minor = raw.IS_MINOR.astype(bool)
    suffix = raw.MUTATION.where(~minor, raw.MINOR_MUTATION)
    muts = raw[["UNIQUEID", "GENE"]].assign(
        MUT=raw.GENE + "@" + suffix,
        FRSf=raw.FRS.where(minor, 1.0),
    )
    flags = pd.DataFrame(suffix.map(classify).tolist(), index=muts.index)
    core, aggregates = sample_tables(
        pd.concat([muts, flags], axis=1), ids, features["frs_min"], tuple(features["dbd_range"]),
    )
    return core, aggregates, y


def metrics(
    y: np.ndarray,
    probability: np.ndarray,
    binary: np.ndarray,
    spec_targets: tuple[float, ...],
) -> dict[str, float]:
    return {
        "auc": float(roc_auc_score(y, probability)),
        "binary_sens": float(binary[y == 1].mean()),
        "binary_spec": float((~binary[y == 0]).mean()),
        **{f"sens_at_spec_{sp}": sens_at_spec(y, probability, sp) for sp in spec_targets},
    }


def rank_results(rows: list[dict], golden_name: str = "lr_model") -> tuple[pd.DataFrame, list[str]]:
    ranked = pd.DataFrame(rows).sort_values(["S1_auc", "KM_auc", "model"], ascending=[False, False, True])
    for cohort in ("S1", "KM"):
        ranked[f"{cohort}_rank"] = ranked[f"{cohort}_auc"].rank(method="min", ascending=False).astype(int)
    golden = ranked.loc[ranked.model == golden_name].iloc[0]
    ranked["S1_auc_delta"] = ranked.S1_auc - golden.S1_auc
    ranked["KM_auc_delta"] = ranked.KM_auc - golden.KM_auc
    qualifies = (ranked.S1_auc > golden.S1_auc + 0.005) & (ranked.KM_auc >= golden.KM_auc)
    return ranked, ranked.loc[qualifies, "model"].tolist()


def main() -> None:
    cfg = yaml.safe_load(Path("config.yaml").read_text(encoding="utf-8"))
    proc, results = (Path(cfg["paths"][key]) for key in ("processed", "results"))
    features = cfg["features"]
    samples = pd.read_parquet(proc / "samples.parquet")
    muts = pd.read_parquet(proc / "mutations.parquet")
    fold_ids = np.load(proc / "s1_folds.npy", allow_pickle=False)
    ids, y, groups = (samples[key].to_numpy() for key in ("UNIQUEID", "y", "SUBJID"))
    assert len(ids) == len(np.unique(ids)) == len(fold_ids) == 11_868 and y.sum() == 866
    assert set(fold_ids) == set(range(5)), "Expected all five saved S1 folds"
    assert pd.DataFrame({"group": groups, "fold": fold_ids}).groupby("group").fold.nunique().eq(1).all()

    core, aggregates = sample_tables(muts, ids, features["frs_min"], tuple(features["dbd_range"]))
    fold_data = []
    for fold in range(5):
        te = fold_ids == fold
        vocab = vocab_from(core, ~te, 3)
        X = design(core, aggregates, vocab, "F2").to_numpy()
        fold_data.append((te, X))
    full_vocab = vocab_from(core, np.ones(len(y), dtype=bool), 3)
    X_full = design(core, aggregates, full_vocab, "F2").to_numpy()
    print("Reading the fixed KM cohort and checking training overlap", flush=True)
    km_core, km_aggregates, km_y = km_tables(ids, features)
    X_km = design(km_core, km_aggregates, full_vocab, "F2").to_numpy()
    print("KM overlap check passed. Starting the preregistered grids.", flush=True)

    cv_baseline = json.loads(
        (results / "models/lr_model/e04_s1_patient_cv/results.json").read_text(encoding="utf-8")
    )["observed"]
    km_baseline = json.loads(
        (results / "approaches/external_km/e_km_validation/results.json").read_text(encoding="utf-8")
    )["observed"]
    rows = [{
        "model": "lr_model", "S1_auc": cv_baseline["auc"],
        "S1_sens_at_spec_0.985": cv_baseline["sens_at_spec_0.985"],
        "KM_auc": km_baseline["model_auc"],
        "KM_sens_at_spec_0.95": km_baseline["model_sens_at_spec_0.95"],
        "best_params": {"C": cfg["model"]["C"], "penalty": "l2"},
    }]
    for name, (estimator, grid) in CANDIDATES.items():
        P = np.empty(len(y))
        binary, firm = np.empty(len(y), dtype=bool), np.empty(len(y), dtype=bool)
        folds = []
        for fold, (te, X) in enumerate(fold_data):
            print(f"{name} S1 fold {fold + 1}/5", flush=True)
            tr = ~te
            final, platt, limits, params = fit_tuned(estimator, grid, X[tr], y[tr], groups[tr])
            _, probability = score_est(final, platt, X[te])
            P[te] = probability
            binary[te] = probability >= limits["t_bin"]
            firm[te] = (probability <= limits["lo"]) | (probability >= limits["hi"])
            folds.append({"fold": fold, **params, **limits})
        cv_obs = {
            **metrics(y, P, binary, (0.95, 0.97, 0.985)),
            "auprc": float(average_precision_score(y, P)), "ternary_firm": float(firm.mean()),
        }
        write_result(
            name, "e21_tuned_cv", cv_obs, results_root=results,
            tables={"metrics": pd.DataFrame([cv_obs]), "best_params_per_fold": pd.DataFrame(folds)},
            extra={"best_params_per_fold": folds, "grid": grid},
            method=f"{METHOD} Evaluation uses the five saved S1 patient-grouped test folds.",
            caveats=CAVEATS, reproduce=REPRODUCE,
        )

        print(f"{name} all golden training for KM", flush=True)
        final, platt, limits, params = fit_tuned(estimator, grid, X_full, y, groups)
        _, km_probability = score_est(final, platt, X_km)
        km_obs = {
            **metrics(km_y, km_probability, km_probability >= limits["t_bin"], (0.95,)),
            "n_isolates": len(km_y), "n_resistant": int(km_y.sum()), "zero_training_overlap": True,
        }
        write_result(
            name, "e21_km_external", km_obs, results_root=results,
            tables={"metrics": pd.DataFrame([km_obs]), "best_params": pd.DataFrame([params])},
            extra={"best_params": params, "thresholds": limits},
            method=f"{METHOD} Refit using the entire golden cohort and evaluate on KMVALIDATION2025.",
            caveats=f"{CAVEATS} {KM_CAVEATS}", reproduce=REPRODUCE,
        )
        rows.append({
            "model": name, "S1_auc": cv_obs["auc"],
            "S1_sens_at_spec_0.985": cv_obs["sens_at_spec_0.985"],
            "KM_auc": km_obs["auc"], "KM_sens_at_spec_0.95": km_obs["sens_at_spec_0.95"],
            "best_params": params,
        })

    ranked, qualifying = rank_results(rows)
    if qualifying:
        decision = f"D013 candidates for E19 review: {', '.join(qualifying)}."
    else:
        decision = "Keep golden LR (F2). No tuned candidate met both D013 criteria."
    comparison = ranked.assign(best_params=ranked.best_params.map(lambda p: json.dumps(p, sort_keys=True)))
    write_result(
        "head_to_head", "e21_comparison",
        {"S1_leader": ranked.iloc[0].model, "KM_leader": ranked.loc[ranked.KM_rank == 1, "model"].tolist()},
        results_root=results, tables={"comparison": comparison},
        extra={"rankings": ranked.to_dict("records"), "qualifying": qualifying, "decision": decision},
        method=(
            "Rank S1 and KM separately, using the saved golden LR references. "
            "D013 requires an S1 AUC gain greater than 0.005 and KM AUC at least as high as golden. "
            "The comparison shows the settings selected on all golden training for the KM evaluation."
        ),
        interpretation=decision, caveats=f"{CAVEATS} {KM_CAVEATS}", reproduce=REPRODUCE,
    )
    print(comparison.to_string(index=False), flush=True)
    print(decision, flush=True)


if __name__ == "__main__":
    main()
