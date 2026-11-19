"""E03: catalogue cross-validation on S1, build per fold (AGENTS.md sections 6-7)."""
import sys
import time
from pathlib import Path

sys.path.insert(0, "src")
import numpy as np
import pandas as pd
import yaml
from sklearn.model_selection import GroupKFold

from discern.catalogue import apply_catalogue, build_catalogue
from discern.metrics import ternary_metrics
from discern.results_io import write_result

cfg = yaml.safe_load(open("config.yaml"))
raw, proc = Path(cfg["paths"]["raw"]), Path(cfg["paths"]["processed"])
s = pd.read_parquet(proc / "samples.parquet")
m = pd.read_parquet(proc / "mutations.parquet")
wild = raw / cfg["catalogue"]["wildcards"]
cache_dir = proc / "cache"
cache_dir.mkdir(parents=True, exist_ok=True)

folds = np.full(len(s), -1)
calls = np.empty(len(s), dtype=object)
for k, (tr, te) in enumerate(GroupKFold(n_splits=5).split(s, groups=s.SUBJID)):
    t0 = time.time()
    assert set(s.SUBJID.iloc[tr]).isdisjoint(set(s.SUBJID.iloc[te]))
    rc = build_catalogue(s.iloc[tr], m, str(cache_dir / f"cat_fold{k}.csv"), str(wild))
    calls[te] = apply_catalogue(rc, s.iloc[te], m)
    folds[te] = k
    print(f"fold {k}: train {len(tr)} test {len(te)} {time.time()-t0:.0f}s", flush=True)

np.save(proc / "s1_folds.npy", folds)
pd.DataFrame({"UNIQUEID": s.UNIQUEID, "fold": folds, "cat_call": calls}).to_parquet(proc / "e03_oof_calls.parquet")
obs = ternary_metrics(s.y, calls)
EXP = {"firm": 0.944, "sens_firm": 0.792, "spec_firm": 0.984, "sens_U_as_S": 0.604,
       "spec_U_as_S": 0.984, "sens_U_as_R": 0.842, "spec_U_as_R": 0.943, "n_U": 661, "n_U_resistant": 206}
TOL = {"firm": 0.01, "sens_firm": 0.015, "spec_firm": 0.004, "sens_U_as_S": 0.015, "spec_U_as_S": 0.004,
       "sens_U_as_R": 0.015, "spec_U_as_R": 0.01, "n_U": 15, "n_U_resistant": 15}
out = write_result(
    "catalogue_tierA", "e03_s1_catalogue_cv", obs,
    expected=EXP, tolerance=TOL, tables={"metrics": pd.DataFrame([obs])},
    method="catomatic BinaryBuilder per training fold (Binomial, background 0.1, two tails, strict_unlock, "
           "frs 0.1, no_mmpL5 wildcards); piezo predict on test; patient-grouped 5-fold.",
    interpretation="Rebuilt catalogue under patient-grouped CV (stricter than the paper's shuffled CV).",
    caveats="In-sample (built on all data) would be optimistic; this is out-of-fold.",
    reproduce="python eval/e03_catalogue_cv.py",
    supports="Baseline reproduction of Adlard 2025 catalogue performance.",
)
import json
rec = json.loads(out.read_text())
print(json.dumps({"observed": obs, "pass_all": rec["pass_all"]}, indent=1))