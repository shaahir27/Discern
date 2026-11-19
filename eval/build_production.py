"""Build production artefacts (AGENTS.md M5): train on ALL 11,868 samples; save per-model under
models/lr_model/ (D006). Catalogue + calibrated LR + evidence store + metadata + OOF replay file.
"""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import hashlib
import json
from datetime import datetime, timezone

import joblib
import numpy as np
import pandas as pd
import yaml

from discern.catalogue import build_catalogue, cat_muts
from discern.features import design, sample_tables, vocab_from
from discern.model import fit_nested, thresholds
from discern.results_io import _git_commit, _manifest_hash, _versions

cfg = yaml.safe_load(open("config.yaml"))
raw, proc = Path(cfg["paths"]["raw"]), Path(cfg["paths"]["processed"])
models_root = Path(cfg["paths"]["models"])
s = pd.read_parquet(proc / "samples.parquet")
m = pd.read_parquet(proc / "mutations.parquet")
y, ids, g = s.y.values, s.UNIQUEID.values, s.SUBJID.values
mdir = models_root / "lr_model"
mdir.mkdir(parents=True, exist_ok=True)
wild = raw / cfg["catalogue"]["wildcards"]
frs_min = cfg["features"]["frs_min"]
dbd = tuple(cfg["features"]["dbd_range"])

# ---- production catalogue (all samples) ----
build_catalogue(s, m, str(mdir / "catalogue_bdq.csv"), str(wild))
cat_version = {"name": "BDQ-production", "version": "1.0", "reference": "NC_000962.3",
               "built_utc": datetime.now(timezone.utc).isoformat(), "n_samples": int(len(s)),
               "settings": cfg["catalogue"], "data_manifest_sha256": _manifest_hash(Path("results"))}
(mdir / "catalogue_version.json").write_text(json.dumps(cat_version, indent=2), encoding="utf-8")

# ---- production model (all samples, nested calibration) ----
core, a = sample_tables(m, ids, frs_min, dbd)
vocab = vocab_from(core, np.ones(len(s), bool), cfg["features"]["vocab_min_count"])
Xdf = design(core, a, vocab, "F2")
X = Xdf.values
final, platt, poof = fit_nested(X, y, g, cfg["model"]["C"], cfg["model"]["max_iter"])
th = thresholds(poof, y, cfg["thresholds"]["npv_target"], cfg["thresholds"]["ppv_target"], cfg["thresholds"]["min_rows"])
pi_tr = float(y.mean())
artefact = {"final": final, "platt": platt, "vocab": vocab, "feature_names": list(Xdf.columns),
            "thresholds": th, "pi_tr": pi_tr, "fset": "F2", "frs_min": frs_min, "dbd": dbd,
            "C": cfg["model"]["C"], "n_features": int(X.shape[1])}
joblib.dump(artefact, mdir / "bdq_lr.joblib")

# ---- evidence store (per core-gene variant, pooled + by source) ----
y_by, src_by = s.set_index("UNIQUEID").y, s.set_index("UNIQUEID").SRC
mall = cat_muts(m[m.FRSf >= frs_min]).drop_duplicates(["UNIQUEID", "MUT"]).copy()
mall["y"] = mall.UNIQUEID.map(y_by)
mall["SRC"] = mall.UNIQUEID.map(src_by)


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 1.0
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return float(max(0, c - h)), float(min(1, c + h))


ev = []
for mut, d in mall.groupby("MUT"):
    n = int(d.UNIQUEID.nunique())
    nr = int(d[d.y == 1].UNIQUEID.nunique())
    lo, hi = wilson(nr, n)
    row = {"variant": mut, "gene": mut.split("@")[0], "mutation": mut.split("@", 1)[1],
           "n_carriers": n, "n_resistant": nr, "pct_resistant": nr / n if n else 0.0,
           "ci95_low": lo, "ci95_high": hi}
    for src in ("CRyPTIC", "NICD"):
        ds = d[d.SRC == src]
        row[f"{src}_n"] = int(ds.UNIQUEID.nunique())
        row[f"{src}_r"] = int(ds[ds.y == 1].UNIQUEID.nunique())
    ev.append(row)
ev_df = pd.DataFrame(ev).sort_values("n_carriers", ascending=False)
ev_df.to_csv(mdir / "evidence_store.csv", index=False)

# ---- metadata ----
def fhash(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16] if Path(p).exists() else None


metadata = {"model": "lr_model", "package_versions": _versions(), "git_commit": _git_commit(),
            "data_manifest_sha256": _manifest_hash(Path("results")), "pi_tr": pi_tr,
            "n_samples": int(len(s)), "n_resistant": int(y.sum()), "thresholds": th,
            "n_features": int(X.shape[1]), "feature_names": list(Xdf.columns),
            "catalogue_file": "catalogue_bdq.csv", "model_file": "bdq_lr.joblib",
            "evidence_file": "evidence_store.csv",
            "results_hashes": {"e03": fhash("results/approaches/catalogue/e03_s1_catalogue_cv/results.json"),
                               "e04": fhash("results/models/lr_model/e04_s1_patient_cv/results.json")},
            "built_utc": datetime.now(timezone.utc).isoformat()}
(mdir / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

# ---- OOF replay file (held-out calls/probs, NOT production) ----
cat = pd.read_parquet(proc / "e03_oof_calls.parquet")
oof = pd.read_parquet(proc / "e04_oof.parquet")
replay = s[["UNIQUEID", "SRC", "SITEID", "y"]].merge(cat[["UNIQUEID", "fold", "cat_call"]], on="UNIQUEID").merge(
    oof[["UNIQUEID", "p"]], on="UNIQUEID")
replay.to_parquet("results/oof_predictions.parquet")

# ---- demo files ----
demo = Path("demo")
demo.mkdir(exist_ok=True)
mv = m[m.FRSf >= frs_min][["UNIQUEID", "GENE", "MUT", "FRSf"]].copy()
mv["mutation"] = mv.MUT.str.split("@").str[-1]
pick = replay.groupby("cat_call", group_keys=False).apply(lambda d: d.head(4), include_groups=False)
pick_ids = list(pick.UNIQUEID)[:12]
held = mv[mv.UNIQUEID.isin(pick_ids)].merge(replay[["UNIQUEID", "cat_call", "p", "y"]], on="UNIQUEID")
held = held.rename(columns={"UNIQUEID": "sample_id", "GENE": "gene", "FRSf": "frs"})
held[["sample_id", "gene", "mutation", "frs", "cat_call", "p", "y"]].to_csv(demo / "held_out_samples.csv", index=False)
ex = mv[mv.UNIQUEID.isin(pick_ids[:4])].rename(columns={"UNIQUEID": "sample_id", "GENE": "gene", "FRSf": "frs"})
ex[["sample_id", "gene", "mutation", "frs"]].to_csv(demo / "example_batch.csv", index=False)

print(f"production built: n_features={X.shape[1]} vocab={len(vocab)} pi_tr={pi_tr:.4f} thresholds={th}")
print(f"evidence store: {len(ev_df)} variants; replay: {len(replay)} samples; demo held-out: {len(pick_ids)}")
print(json.dumps({k: metadata[k] for k in ("n_samples", "n_resistant", "n_features", "thresholds")}, indent=1))