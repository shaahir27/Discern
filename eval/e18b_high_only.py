"""E18b: HIGH-quality-only expanded training (isolates data quantity from MEDIUM label noise).

Reuses expanded parquets from E18. Compares HIGH-only expansion (16,947 isolates / 206 R) to golden
(11,868/866) and to the HIGH+MEDIUM arm. Golden + production untouched; KM held out. Descriptive.
"""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import json
import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

from discern.features import design, sample_tables, vocab_from
from discern.model import fit_nested, score, sens_at_spec, thresholds
from discern.results_io import write_result
from discern.runner import nested_oof
from discern.splits import cross_source_folds

cfg = yaml.safe_load(open("config.yaml"))
proc, ext = Path(cfg["paths"]["processed"]), Path("data/external_cryptic_v34")
frs_min, dbd = cfg["features"]["frs_min"], tuple(cfg["features"]["dbd_range"])
C, MAXIT, TH = cfg["model"]["C"], cfg["model"]["max_iter"], cfg["thresholds"]

s_all = pd.read_parquet(proc / "expanded_samples.parquet")
m = pd.read_parquet(proc / "expanded_mutations.parquet")
s = s_all[s_all.QUALITY == "HIGH"].reset_index(drop=True)
ids, y, g = s.UNIQUEID.values, s.y.values, s.SUBJID.values
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
auc = float(roc_auc_score(y, P))
obs = {"n_samples": int(len(s)), "n_resistant": int(y.sum()), "auc": auc,
       "binary_sens": float((pb & (y == 1)).sum() / nR), "binary_spec": float(((~pb) & (y == 0)).sum() / nS),
       "sens_at_spec_0.95": sens_at_spec(y, P, 0.95), "sens_at_spec_0.985": sens_at_spec(y, P, 0.985),
       "golden_auc": 0.908539, "high_medium_auc": 0.8381}

mask = s.SRC.isin(["CRyPTIC", "NICD"]).values
s2 = s[mask].reset_index(drop=True)
core2, a2 = sample_tables(m, s2.UNIQUEID.values, frs_min, dbd)
cs = {}
if {"CRyPTIC", "NICD"} <= set(s2.SRC):
    for name, tr, te in cross_source_folds(s2):
        if len(te) == 0 or len(tr) == 0:
            continue
        Pc, _, thc = nested_oof(core2, a2, s2.y.values, s2.SUBJID.values, [(tr, te)], fset="F2", C=C, max_iter=MAXIT)
        p, yt = Pc[te], s2.y.values[te]
        pbc = p >= thc[0]["t_bin"]
        cs[f"{name}_auc"] = float(roc_auc_score(yt, p))
        cs[f"{name}_sens"] = float((pbc & (yt == 1)).sum() / max((yt == 1).sum(), 1))
        cs[f"{name}_spec"] = float(((~pbc) & (yt == 0)).sum() / max((yt == 0).sum(), 1))
else:
    cs["cross_source_note"] = "NICD is MEDIUM quality; absent from HIGH-only arm -> cross-source N/A"
obs.update(cs)

# KM external
d = pd.read_parquet(ext / "DST_MEASUREMENTS.parquet").reset_index()
km = d[(d.SOURCE == "KMVALIDATION2025") & (d.DRUG == "BDQ")].drop_duplicates("UNIQUEID").copy()
km["y"] = (km.PHENOTYPE == "R").astype(int)
raw_mut = pd.read_parquet(ext / "MUTATIONS.parquet", filters=[("GENE", "in", ["Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5"])]).reset_index()
from discern.grammar import classify
r = raw_mut[raw_mut.UNIQUEID.isin(set(km.UNIQUEID))].copy()
minor = r.IS_MINOR.astype(bool)
suf = np.where(minor, r.MINOR_MUTATION.astype(str), r.MUTATION.astype(str))
mk = pd.DataFrame({"UNIQUEID": r.UNIQUEID.values, "GENE": r.GENE.values, "MUT": r.GENE.astype(str).values + "@" + suf,
                   "FRSf": np.where(minor, r.FRS.values, 1.0), "suffix": suf})
fl = mk.suffix.map(classify)
for k in ("IS_SYN", "IS_NULL", "is_stop", "is_indel", "is_prom", "is_mis", "pos"):
    mk[k] = [f[k] for f in fl]
vocab_all = vocab_from(core, np.ones(len(s), bool), cfg["features"]["vocab_min_count"])
final, platt, poof = fit_nested(design(core, a, vocab_all, "F2").values, y, g, C, MAXIT)
th_all = thresholds(poof, y, TH["npv_target"], TH["ppv_target"], TH["min_rows"])
ck, ak = sample_tables(mk.reset_index(drop=True), km.UNIQUEID.values, frs_min, dbd)
_, Pkm = score(final, platt, design(ck, ak, vocab_all, "F2").values)
obs["km_external_auc"] = float(roc_auc_score(km.y.values, Pkm))
obs["golden_production_km_auc"] = 0.8681826759765995

write_result("lr_expanded", "e18b_high_only", obs, tables={"metrics": pd.DataFrame([obs])},
             method="HIGH-quality-only expanded CRyPTIC v3.4 BDQ cohort (exclude KM). LR F2, S1 patient-grouped "
                    "nested CV + cross-source + KM external. Isolates data quantity from MEDIUM label noise.",
             interpretation=f"HIGH-only expansion ({len(s)}/{int(y.sum())} R): S1 AUC {auc:.4f} "
                            f"(golden 0.9085, HIGH+MEDIUM 0.8381); KM external AUC {obs['km_external_auc']:.4f} "
                            f"(golden production 0.8682). Cross-source: "
                            f"{cs.get('NICD_to_CRyPTIC_sens', 'N/A (NICD is MEDIUM, excluded from HIGH-only)')}.",
             caveats="Descriptive; different cohort than golden reproduction. KM assay undocumented; 35 R wide CI. "
                     "Golden + production untouched; no tuning. Tests if CLEAN richer data helps generalization.",
             reproduce="python eval/e18b_high_only.py",
             supports="Data-richness ablation: quality-controlled expansion vs golden.")
print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in obs.items()}, indent=1))