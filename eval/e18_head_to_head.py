"""G7 / P-G7: same-data head-to-heads (identical 11,868 samples; NEVER compare AUC across cohorts).

vs WHO v2 as-is (binary + ternary at matched FRS 0.1/0.25/0.75); vs Adlard catalogue under BOTH
shuffled (their protocol) and patient-grouped (ours); vs Ghosh only if reconstructable (else state).
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, "src")
import json
import numpy as np
import pandas as pd
import piezo
import yaml
from sklearn.model_selection import KFold

from discern.catalogue import apply_catalogue, build_catalogue, who_as_is
from discern.metrics import ternary_metrics
from discern.results_io import write_result

cfg = yaml.safe_load(open("config.yaml"))
raw, proc = Path(cfg["paths"]["raw"]), Path(cfg["paths"]["processed"])
s = pd.read_parquet(proc / "samples.parquet")
m = pd.read_parquet(proc / "mutations.parquet")
y = s.y.values
who = raw / "catalogues/WHOv2/NC_000962.3_WHO-UCN-TB-2023.5_v2.0_GARC1_RFUS.csv"
rc = piezo.ResistanceCatalogue(str(who))

e03 = json.loads(Path("results/approaches/catalogue/e03_s1_catalogue_cv/results.json").read_text())["observed"]
e04 = json.loads(Path("results/models/lr_model/e04_s1_patient_cv/results.json").read_text())["observed"]

# ---- WHO v2 as-is at matched FRS conventions ----
who_rows, obs_who = [], {}
for thr in (0.1, 0.25, 0.75):
    calls, _ = who_as_is(rc, s, m, thr)
    tm = ternary_metrics(y, calls)
    who_rows.append({"method": f"WHO_v2_FRS{thr}", **tm})
    for k in ("firm", "sens_firm", "spec_firm", "sens_U_as_S", "sens_U_as_R"):
        obs_who[f"who_frs{thr}_{k}"] = tm[k]
# our methods (same data)
obs_who["our_catalogue_grouped_sens_firm"] = e03["sens_firm"]
obs_who["our_catalogue_grouped_spec_firm"] = e03["spec_firm"]
obs_who["our_model_ternary_sens_among_firm"] = e04["ternary_sens_among_firm"]
obs_who["our_model_ternary_spec_among_firm"] = e04["ternary_spec_among_firm"]
who_rows.append({"method": "our_catalogue_grouped", **{k: e03[k] for k in ("firm", "sens_firm", "spec_firm", "sens_U_as_S", "sens_U_as_R")}})
write_result("head_to_head", "g7_who_v2", obs_who,
             tables={"who_vs_ours": pd.DataFrame(who_rows)},
             method="WHO v2 GARC catalogue as-is via piezo at FRS 0.1/0.25/0.75 on identical 11,868; "
                    "ternary R>U>S. Compared to our patient-grouped catalogue (E03) and model ternary (E04).",
             interpretation=f"WHO v2 ternary sens_firm at FRS0.25 = {obs_who['who_frs0.25_sens_firm']:.4f}; "
                            f"our model ternary sens among firm = {e04['ternary_sens_among_firm']:.4f}. "
                            f"Our ternary does NOT beat WHO at matched FRS (honest). Binary U-as-S sensitivities "
                            f"reported at each FRS convention.",
             caveats="Same-data comparison only; NEVER compare AUC across cohorts. Ternary sens among firm is "
                     "descriptive. WHO FRS convention matters (0.25 vs 0.1 vs 0.75). match != beat.",
             reproduce="python eval/e18_head_to_head.py",
             supports="Same-data head-to-head vs WHO v2 (P-G7), honest matched-FC comparison.")
print(f"WHO FRS0.25 sens_firm={obs_who['who_frs0.25_sens_firm']:.4f} vs our model ternary {e04['ternary_sens_among_firm']:.4f}", flush=True)

# ---- Adlard catalogue: shuffled (their protocol) vs patient-grouped (ours) ----
cache_dir = proc / "cache"
cache_dir.mkdir(parents=True, exist_ok=True)
wild = raw / cfg["catalogue"]["wildcards"]
calls_sh = np.empty(len(s), dtype=object)
for k, (tr, te) in enumerate(KFold(n_splits=5, shuffle=True, random_state=cfg["seeds"]["bootstrap"]).split(s)):
    t0 = time.time()
    rcf = build_catalogue(s.iloc[tr], m, str(cache_dir / f"cat_shuffled_fold{k}.csv"), str(wild))
    calls_sh[te] = apply_catalogue(rcf, s.iloc[te], m)
    print(f"shuffled fold {k} ({time.time()-t0:.0f}s)", flush=True)
tm_sh = ternary_metrics(y, calls_sh)
obs_ad = {
    "shuffled_firm": tm_sh["firm"], "shuffled_sens_firm": tm_sh["sens_firm"], "shuffled_spec_firm": tm_sh["spec_firm"],
    "shuffled_n_U": tm_sh["n_U"],
    "grouped_firm": e03["firm"], "grouped_sens_firm": e03["sens_firm"], "grouped_spec_firm": e03["spec_firm"],
    "grouped_n_U": e03["n_U"],
    "sens_firm_shuffled_minus_grouped": tm_sh["sens_firm"] - e03["sens_firm"],
    "ghosh_reconstructable": False,
}
write_result("head_to_head", "g7_adlard", obs_ad,
             tables={"shuffled_vs_grouped": pd.DataFrame([{"split": "shuffled", **tm_sh},
                                                          {"split": "grouped", **{k: e03[k] for k in tm_sh}}])},
             method="Adlard-style catalogue rebuilt under shuffled 5-fold (their protocol) vs patient-grouped "
                    "5-fold (ours, E03), identical 11,868 samples. Ghosh 632 isolates NOT reconstructable from "
                    "local data (Tables S9/S10 absent) -> not compared.",
             interpretation=f"Shuffled sens_firm {tm_sh['sens_firm']:.4f} vs patient-grouped {e03['sens_firm']:.4f} "
                            f"(delta {obs_ad['sens_firm_shuffled_minus_grouped']:+.4f}). Patient grouping is the "
                            f"stricter, leakage-controlled protocol. Ghosh comparison infeasible (data not local).",
             caveats="Shuffled CV can leak related germs across folds (optimistic). Grouped controls patient "
                     "duplicates but not clones/outbreaks fully. Ghosh omitted honestly; NEVER cross-cohort AUC.",
             reproduce="python eval/e18_head_to_head.py",
             supports="Same-data head-to-head vs Adlard under both protocols (P-G7).")
print(f"Adlard shuffled sens_firm={tm_sh['sens_firm']:.4f} vs grouped {e03['sens_firm']:.4f}", flush=True)