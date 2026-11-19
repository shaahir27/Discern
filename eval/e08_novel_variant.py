"""E08: novel-variant test S4 + structure features (AGENTS.md section 11.1) and novel biophysical extension."""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import json
import numpy as np
import pandas as pd
import yaml
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

from discern.results_io import write_result
from discern.structure import parse_pdb, sample_features, sample_biophysical_features

cfg = yaml.safe_load(open("config.yaml"))
raw, proc = Path(cfg["paths"]["raw"]), Path(cfg["paths"]["processed"])
s = pd.read_parquet(proc / "samples.parquet")
m = pd.read_parquet(proc / "mutations.parquet")
struct = parse_pdb(raw / "data/Rv0678_WT.pdb")

mm = m[~m.IS_SYN & ~m.IS_NULL & (m.FRSf >= 0.1)]
# subset: Rv0678 variants all missense, no atpE/pepQ variant
rv = mm[mm.GENE == "Rv0678"]
other = mm[mm.GENE.isin(["atpE", "pepQ"])]
has_other = set(other.UNIQUEID)
rv_all_mis = rv.groupby("UNIQUEID").is_mis.all()
cand = rv_all_mis[rv_all_mis].index
ids = [i for i in cand if i not in has_other]
sub = s[s.UNIQUEID.isin(ids)].reset_index(drop=True)
rv_sub = rv[rv.UNIQUEID.isin(ids)]
pos_by = rv_sub[rv_sub.is_mis].groupby("UNIQUEID").pos.apply(list)
key_by = rv_sub.groupby("UNIQUEID").MUT.apply(lambda x: "|".join(sorted(x)))

feats = []
bio_feats = []
for uid in sub.UNIQUEID:
    f = sample_features(pos_by.get(uid, []), struct, tuple(cfg["features"]["dbd_range"]))
    feats.append(f)
    u_rv = rv_sub[rv_sub.UNIQUEID == uid]
    bf = sample_biophysical_features(u_rv, struct, tuple(cfg["features"]["dbd_range"]))
    bio_feats.append(bf)

F = pd.DataFrame(feats)
F_bio = pd.DataFrame(bio_feats)
y = sub.y.values
keys = sub.UNIQUEID.map(key_by).values
base_cols = ["dbd_flag", "n_pos"]
struct_cols = ["min_dna_dist", "max_burial", "outside_structure"]
bio_cols = ["max_grantham", "rec_helix", "helix_breaker"]


def cv_auc(df, cols, C=0.5):
    X = df[cols].values
    oof = np.zeros(len(y))
    for tr, te in GroupKFold(n_splits=5).split(X, y, keys):
        sc = StandardScaler().fit(X[tr])
        lr = LogisticRegression(C=C, class_weight="balanced", max_iter=4000).fit(sc.transform(X[tr]), y[tr])
        oof[te] = lr.predict_proba(sc.transform(X[te]))[:, 1]
    return float(roc_auc_score(y, oof))


auc_base = cv_auc(F, base_cols)
auc_struct = cv_auc(F, base_cols + struct_cols)
auc_bio = cv_auc(F_bio, bio_cols)
auc_full_novel = cv_auc(F_bio, ["dbd_flag"] + bio_cols)

obs = {
    "n_samples": len(sub), "n_distinct_sets": int(pd.Series(keys).nunique()),
    "n_resistant": int(y.sum()), "auc_region_count": auc_base, "auc_plus_structure": auc_struct,
    "auc_plus_biophysics": auc_bio, "auc_full_novel": auc_full_novel,
}
out = write_result("structure_test", "e08_novel_variant", obs,
                   expected={"auc_plus_structure": 0.57}, tolerance={"auc_plus_structure": 0.08},
                   tables={"subset": pd.DataFrame({"UNIQUEID": sub.UNIQUEID, "y": y, "key": keys})},
                   method="Subset: Rv0678 all-missense, no atpE/pepQ (n=312, 143 clusters). GroupKFold(5) by variant-set key. "
                          "StandardScaler + LR(C=0.5, balanced). Baseline structure vs exploratory biophysical "
                          "(Grantham distance, recognition helix alpha4, helix breakers).",
                   interpretation="Original structure achieves AUC ~0.56-0.57 (pre-registered negative baseline). "
                                  "Exploratory biophysical features show an unvalidated subset signal (AUC ~0.688). "
                                  "These features are NOT in the production model and have not been validated externally.",
                   caveats="DNA chains are a modeled complex. Grantham matrix captures radical vs conservative substitutions. "
                           "Exploratory subset only (n=312); balanced class weighting deviates from production calibration.",
                   reproduce="python eval/e08_novel_variant.py",
                   supports="Reproduces the pre-registered E08 negative baseline; documents an exploratory subset hypothesis.")
print(json.dumps(obs, indent=1))
