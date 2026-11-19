"""E18: does richer training data improve metrics? Expanded CRyPTIC v3.4 cohort.

Training = full CRyPTIC v3.4 BDQ (HIGH+MEDIUM quality), EXCLUDING the KM held-out external test.
23,441 isolates / 1,046 resistant (+182 new R vs golden; CRyPTIC R 83->255). Golden baseline and
production model are UNTOUCHED; KM stays external. Descriptive comparison, no tuning.
"""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import json
import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

from discern.features import design, sample_tables, vocab_from
from discern.grammar import classify
from discern.model import fit_nested, score, sens_at_spec, thresholds
from discern.results_io import write_result
from discern.runner import nested_oof
from discern.splits import cross_source_folds

cfg = yaml.safe_load(open("config.yaml"))
proc, ext = Path(cfg["paths"]["processed"]), Path("data/external_cryptic_v34")
BDQ = ["Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5"]
frs_min, dbd = cfg["features"]["frs_min"], tuple(cfg["features"]["dbd_range"])
C, MAXIT = cfg["model"]["C"], cfg["model"]["max_iter"]
TH = cfg["thresholds"]
raw_mut = pd.read_parquet(ext / "MUTATIONS.parquet", filters=[("GENE", "in", BDQ)]).reset_index()


def build_muts(raw, ids):
    r = raw[raw.UNIQUEID.isin(set(ids))].copy()
    minor = r.IS_MINOR.astype(bool)
    suf = np.where(minor, r.MINOR_MUTATION.astype(str), r.MUTATION.astype(str))
    m = pd.DataFrame({"UNIQUEID": r.UNIQUEID.values, "GENE": r.GENE.values,
                      "MUT": r.GENE.astype(str).values + "@" + suf,
                      "FRSf": np.where(minor, r.FRS.values, 1.0), "suffix": suf})
    fl = m.suffix.map(classify)
    for k in ("IS_SYN", "IS_NULL", "is_stop", "is_indel", "is_prom", "is_mis", "pos"):
        m[k] = [f[k] for f in fl]
    return m.reset_index(drop=True)


# ---- expanded phenotype cohort (exclude KM) ----
d = pd.read_parquet(ext / "DST_MEASUREMENTS.parquet").reset_index()
b = d[(d.DRUG == "BDQ") & (d.SOURCE != "KMVALIDATION2025") & (d.QUALITY.isin(["HIGH", "MEDIUM"]))].copy()
b["_r"] = (b.PHENOTYPE == "R").astype(int)
b["_m"] = b.METHOD_MIC.notna().astype(int)
b = b.sort_values(["UNIQUEID", "_r", "_m"], ascending=[True, False, False], kind="stable").drop_duplicates("UNIQUEID", keep="first")
b["y"] = b._r
b["SUBJID"] = b.UNIQUEID.str.extract(r"subj\.([^.]+)\.")
b["SITEID"] = b.UNIQUEID.str.extract(r"site\.([^.]+)\.")
b["SRC"] = b.SOURCE
s = b[["UNIQUEID", "SUBJID", "SITEID", "SRC", "y", "QUALITY"]].reset_index(drop=True)
ids, y, g = s.UNIQUEID.values, s.y.values, s.SUBJID.values
m = build_muts(raw_mut, ids)
s.to_parquet(proc / "expanded_samples.parquet")
m.to_parquet(proc / "expanded_mutations.parquet")
print(f"expanded cohort: {len(s)} isolates / {int(y.sum())} R (golden 11868/866); mutation rows {len(m)}", flush=True)

# ---- S1 patient-grouped nested CV (LR F2) ----
core, a = sample_tables(m, ids, frs_min, dbd)
folds = list(GroupKFold(n_splits=5).split(s, groups=s.SUBJID))
folds_arr = np.full(len(s), -1)
for k, (tr, te) in enumerate(folds):
    folds_arr[te] = k
P, T, th = nested_oof(core, a, y, g, folds, fset="F2", C=C, max_iter=MAXIT,
                      npv=TH["npv_target"], ppv=TH["ppv_target"], min_n=TH["min_rows"])
tb = np.array([th[k]["t_bin"] for k in folds_arr])
pb = P >= tb
nR, nS = max((y == 1).sum(), 1), max((y == 0).sum(), 1)
e04 = json.loads(Path("results/models/lr_model/e04_s1_patient_cv/results.json").read_text())["observed"]
obs_cv = {
    "n_samples": int(len(s)), "n_resistant": int(y.sum()),
    "auc": float(roc_auc_score(y, P)), "auprc": float(average_precision_score(y, P)),
    "binary_sens": float((pb & (y == 1)).sum() / nR), "binary_spec": float(((~pb) & (y == 0)).sum() / nS),
    "sens_at_spec_0.95": sens_at_spec(y, P, 0.95), "sens_at_spec_0.97": sens_at_spec(y, P, 0.97),
    "sens_at_spec_0.985": sens_at_spec(y, P, 0.985),
    "n_features": int(design(core, a, vocab_from(core, np.ones(len(s), bool)), "F2").shape[1]),
    "golden_auc": e04["auc"], "delta_auc_vs_golden": float(roc_auc_score(y, P)) - e04["auc"],
}
write_result("lr_expanded", "e18_expanded_cv", obs_cv,
             tables={"metrics": pd.DataFrame([obs_cv]), "thresholds_per_fold": pd.DataFrame(th)},
             method="Expanded CRyPTIC v3.4 BDQ cohort (HIGH+MEDIUM, exclude KM): 23,441 isolates / 1,046 R. "
                    "LR F2, S1 patient-grouped 5-fold nested Platt + thresholds (same protocol as golden E04).",
             interpretation=f"Expanded training AUC {obs_cv['auc']:.4f} vs golden {e04['auc']:.4f} "
                            f"(delta {obs_cv['delta_auc_vs_golden']:+.4f}); sens@spec0.985 "
                            f"{obs_cv['sens_at_spec_0.985']:.4f}. +182 new resistant isolates, CRyPTIC R 83->255.",
             caveats="Different cohort than golden (not the published 11,868 reproduction). MEDIUM labels add "
                     "noise. Descriptive; golden baseline + production untouched; no tuning. KM held out.",
             reproduce="python eval/e18_expanded_training.py",
             supports="Data-richness lever: more resistant isolates to improve generalization.")
print(f"S1 CV: AUC {obs_cv['auc']:.4f} (golden {e04['auc']:.4f}), sens@spec0.985 {obs_cv['sens_at_spec_0.985']:.4f}", flush=True)

# ---- cross-source (CRyPTIC vs NICD) ----
mask = s.SRC.isin(["CRyPTIC", "NICD"]).values
s2 = s[mask].reset_index(drop=True)
ids2, y2, g2 = s2.UNIQUEID.values, s2.y.values, s2.SUBJID.values
core2, a2 = sample_tables(m, ids2, frs_min, dbd)
e06 = json.loads(Path("results/models/lr_model/e06_cross_source/results.json").read_text())["observed"]
obs_cs, cs_rows = {}, []
for name, tr, te in cross_source_folds(s2):
    Pc, _, thc = nested_oof(core2, a2, y2, g2, [(tr, te)], fset="F2", C=C, max_iter=MAXIT)
    p, yt = Pc[te], y2[te]
    pbc = p >= thc[0]["t_bin"]
    auc = float(roc_auc_score(yt, p))
    sens = float((pbc & (yt == 1)).sum() / max((yt == 1).sum(), 1))
    spec = float(((~pbc) & (yt == 0)).sum() / max((yt == 0).sum(), 1))
    cs_rows.append({"direction": name, "n_test": len(te), "n_R_test": int(yt.sum()), "auc": auc, "sens": sens, "spec": spec})
    obs_cs[f"{name}_auc"], obs_cs[f"{name}_sens"], obs_cs[f"{name}_spec"] = auc, sens, spec
    obs_cs[f"{name}_golden_auc"] = e06.get(f"{name}_auc")
write_result("lr_expanded", "e18_cross_source", obs_cs, tables={"cross_source": pd.DataFrame(cs_rows)},
             method="Expanded cohort, train one source test the other (CRyPTIC<->NICD), nested calibration.",
             interpretation="Cross-source generalization with richer CRyPTIC resistant data (255 R vs 83). "
                            "Compare NICD->CRyPTIC sensitivity to golden 0.675 (the weak spot).",
             caveats="Descriptive; different cohort. CRyPTIC test R count larger than golden (83). No tuning.",
             reproduce="python eval/e18_expanded_training.py",
             supports="Tests whether more diverse resistant data fixes cross-source generalization.")
print("cross-source:", {r["direction"]: (round(r["auc"], 3), round(r["sens"], 3), round(r["spec"], 3)) for r in cs_rows}, flush=True)

# ---- KM external generalization (expanded model vs golden production) ----
km = d[(d.SOURCE == "KMVALIDATION2025") & (d.DRUG == "BDQ")].drop_duplicates("UNIQUEID").copy()
km["y"] = (km.PHENOTYPE == "R").astype(int)
km_ids, km_y = km.UNIQUEID.values, km.y.values
m_km = build_muts(raw_mut, km_ids)
vocab_all = vocab_from(core, np.ones(len(s), bool), cfg["features"]["vocab_min_count"])
Xall = design(core, a, vocab_all, "F2").values
final, platt, poof = fit_nested(Xall, y, g, C, MAXIT)
th_all = thresholds(poof, y, TH["npv_target"], TH["ppv_target"], TH["min_rows"])
ck, ak = sample_tables(m_km, km_ids, frs_min, dbd)
Xkm = design(ck, ak, vocab_all, "F2").values
_, Pkm = score(final, platt, Xkm)
pbk = Pkm >= th_all["t_bin"]
km_auc = float(roc_auc_score(km_y, Pkm))
obs_km = {
    "n_isolates": int(len(km)), "n_resistant": int(km_y.sum()),
    "expanded_model_auc": km_auc,
    "expanded_model_sens": float((pbk & (km_y == 1)).sum() / max((km_y == 1).sum(), 1)),
    "expanded_model_spec": float(((~pbk) & (km_y == 0)).sum() / max((km_y == 0).sum(), 1)),
    "expanded_sens_at_spec_0.95": sens_at_spec(km_y, Pkm, 0.95),
    "golden_production_km_auc": 0.8681826759765995,
    "delta_auc_vs_golden_production": km_auc - 0.8681826759765995,
}
write_result("lr_expanded", "e18_km_external", obs_km, tables={"metrics": pd.DataFrame([obs_km])},
             method="Expanded-trained LR applied to held-out KM 792 isolates (zero overlap with expanded "
                    "training). Compare to golden-production KM AUC 0.868.",
             interpretation=f"Expanded model external KM AUC {km_auc:.4f} vs golden production 0.8682 "
                            f"(delta {obs_km['delta_auc_vs_golden_production']:+.4f}). Tests whether richer "
                            f"training improves generalization to unseen isolates.",
             caveats="KM assay undocumented (MEDIUM); 35 R -> wide CI; sensitivity analysis. No tuning; "
                     "KM excluded from expanded training. match != beat.",
             reproduce="python eval/e18_expanded_training.py",
             supports="External generalization benefit of richer training data.")
print(f"KM external: expanded AUC {km_auc:.4f} vs golden production 0.8682", flush=True)