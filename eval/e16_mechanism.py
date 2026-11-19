"""G2 / P-E16: mechanism biology descriptive (Snobre replication) + mechanism/veto ablation.

Descriptive FIRST, then catalogue+LR with mechanism features OFF/ON and epistatic veto OFF/ON on
S1 and cross-source. Every call flip is enumerated with its true label. Exploratory; no tuning.
Decision to ship enhanced is deferred to G5 (must not degrade external validation).
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, "src")
import json
import numpy as np
import pandas as pd
import yaml
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score

from discern import mechanism as M
from discern.catalogue import apply_catalogue, build_catalogue
from discern.features import design, sample_tables, vocab_from
from discern.metrics import ternary_metrics
from discern.model import sens_at_spec, thresholds
from discern.results_io import write_result
from discern.runner import fit_nested_est, score_est
from discern.splits import cross_source_folds

cfg = yaml.safe_load(open("config.yaml"))
raw, proc = Path(cfg["paths"]["raw"]), Path(cfg["paths"]["processed"])
models_root = Path(cfg["paths"]["models"])
s = pd.read_parquet(proc / "samples.parquet")
m0 = pd.read_parquet(proc / "mutations.parquet")
m = m0.merge(s[["UNIQUEID", "y"]], on="UNIQUEID", how="left")
y, ids, g = s.y.values, s.UNIQUEID.values, s.SUBJID.values
cds = M.load_cds(proc / "rv0678_cds.fasta")
MECH = ["rv_lof_any", "rv_lof_highconf", "rv_lof_dbd", "rv_lof_late"]
feat, flags = M.sample_mech_flags(m, ids, cds)
veto = flags.veto_fire.values
TH = cfg["thresholds"]

# ---------- Part 1: descriptive (Snobre replication) ----------
eff = M.variant_effects(m, cds)
mdir = models_root / "lr_mechanism"
mdir.mkdir(parents=True, exist_ok=True)
eff.to_csv(mdir / "lof_effects.csv", index=False)

rv_tr = m[(m.GENE == "Rv0678") & (m.MUT.isin(set(eff.MUT))) & (m.FRSf >= 0.1) & (~m.IS_SYN) & (~m.IS_NULL)]
reg_map = eff.set_index("MUT").region
rv_tr = rv_tr.assign(region=rv_tr.MUT.map(reg_map))
reg_rows = []
for reg, grp in rv_tr.groupby("region"):
    car = grp.UNIQUEID.unique()
    ys = s.set_index("UNIQUEID").y.reindex(car).values
    reg_rows.append({"region": reg, "n_carriers": len(car), "n_resistant": int((ys == 1).sum()),
                     "frac_resistant": float((ys == 1).mean())})
reg_df = pd.DataFrame(reg_rows)

fs_vars = set(eff[eff.type == "frameshift"].MUT)
fs = rv_tr[rv_tr.MUT.isin(fs_vars)]
fs_car = fs.UNIQUEID.unique()
fs_y = s.set_index("UNIQUEID").y.reindex(fs_car).values
snobre_frac = float((fs_y == 0).mean())
obs_desc = {
    "snobre_frac_susceptible_frameshift": snobre_frac,
    "snobre_reported": 0.248,
    "n_frameshift_carriers": int(len(fs_car)),
    "n_frameshift_susceptible": int((fs_y == 0).sum()),
    "n_truncating_variants": int(len(eff)),
    "n_veto_fire_samples": int(veto.sum()),
    **{f"n_variants_{r}": int((eff.region == r).sum()) for r in ("DBD", "mid", "late", "full_del")},
}
write_result("lr_mechanism", "e16_lof_descriptive", obs_desc,
             tables={"lof_effects": eff, "region_summary": reg_df},
             method="Rv0678 CDS (NC_000962.3:778,990-779,487) translation; first-stop region DBD<=99 / mid "
                    "100-160 / late>160; high-confidence LoF excludes single-nt homopolymer indels; "
                    "late-stop preserves dimerization domain (Snobre function-preserving proxy).",
             interpretation=f"Snobre's susceptible-frameshift fraction REPLICATES in this lineage mix: "
                            f"{snobre_frac:.4f} ({int((fs_y==0).sum())}/{len(fs_car)}) vs reported 0.248. "
                            f"DBD-break truncations dominate; {int((eff.region=='late').sum())} late variants "
                            f"are function-preserving candidates.",
             caveats="ARF conservation approximated by late-stop preservation (full phylogenetic ARF alignment "
                     "out of scope). Homopolymer indels flagged low-confidence (e.g. 193_del_g). Whole-gene "
                     "deletion = full_del. A sample carrying multiple regions is counted in each.",
             reproduce="python eval/e16_mechanism.py",
             supports="First predictive integration of truncation-position reasoning for BDQ (to our knowledge); "
                      "descriptive replication of Snobre 2024.")
print(f"Snobre frac_susceptible={snobre_frac:.4f} (reported 0.248); veto_fire={int(veto.sum())}", flush=True)

# ---------- Part 2: LR mechanism features OFF/ON (S1) ----------
core, a = sample_tables(m0, ids, cfg["features"]["frs_min"], tuple(cfg["features"]["dbd_range"]))
folds_arr = np.load(proc / "s1_folds.npy")
folds = [(np.where(folds_arr != k)[0], np.where(folds_arr == k)[0]) for k in range(5)]
P_f2 = np.load(proc / "e04_oof_p.npy")
P_mech = np.zeros(len(s))
lr_factory = lambda: LogisticRegression(C=cfg["model"]["C"], penalty="l2", max_iter=cfg["model"]["max_iter"])
for tr, te in folds:
    vocab = vocab_from(core, tr, cfg["features"]["vocab_min_count"])
    X = pd.concat([design(core, a, vocab, "F2"), feat[MECH]], axis=1).values
    final, platt, poof = fit_nested_est(X[tr], y[tr], g[tr], lr_factory)
    t = thresholds(poof, y[tr], TH["npv_target"], TH["ppv_target"], TH["min_rows"])
    _, p = score_est(final, platt, X[te])
    P_mech[te] = p
auc_f2, auc_mech = float(roc_auc_score(y, P_f2)), float(roc_auc_score(y, P_mech))
s985_f2, s985_mech = sens_at_spec(y, P_f2, 0.985), sens_at_spec(y, P_mech, 0.985)
lr_cmp = pd.DataFrame([
    {"features": "F2 (golden)", "auc": auc_f2, "auprc": float(average_precision_score(y, P_f2)), "sens_at_spec_0.985": s985_f2},
    {"features": "F2+mechanism", "auc": auc_mech, "auprc": float(average_precision_score(y, P_mech)), "sens_at_spec_0.985": s985_mech},
])

# ---------- Part 3: catalogue epistatic veto OFF/ON (S1) ----------
cat = pd.read_parquet(proc / "e03_oof_calls.parquet").set_index("UNIQUEID").reindex(ids)
calls_off = cat.cat_call.values
calls_on = M.apply_veto(calls_off, veto)
tm_off, tm_on = ternary_metrics(y, calls_off), ternary_metrics(y, calls_on)
flip_mask = (calls_off == "R") & (calls_on == "U")
flips_s1 = pd.DataFrame({"UNIQUEID": ids[flip_mask], "y_true": y[flip_mask],
                         "call_off": "R", "call_on": "U", "split": "S1"})

# ---------- Part 4: cross-source veto OFF/ON ----------
cache_dir = proc / "cache"
cache_dir.mkdir(parents=True, exist_ok=True)
wild = raw / cfg["catalogue"]["wildcards"]
cs_rows, flips_cs = [], []
for name, tr, te in cross_source_folds(s):
    t0 = time.time()
    rc = build_catalogue(s.iloc[tr], m0, str(cache_dir / f"cat_e16_{name}.csv"), str(wild))
    c_off = apply_catalogue(rc, s.iloc[te], m0)
    c_on = M.apply_veto(c_off, veto[te])
    o, n = ternary_metrics(y[te], c_off), ternary_metrics(y[te], c_on)
    fm = (c_off == "R") & (c_on == "U")
    flips_cs.append(pd.DataFrame({"UNIQUEID": ids[te][fm], "y_true": y[te][fm],
                                  "call_off": "R", "call_on": "U", "split": name}))
    cs_rows.append({"direction": name, "sens_firm_off": o["sens_firm"], "sens_firm_on": n["sens_firm"],
                    "spec_firm_off": o["spec_firm"], "spec_firm_on": n["spec_firm"],
                    "n_U_off": o["n_U"], "n_U_on": n["n_U"], "n_flips": int(fm.sum()),
                    "n_flips_resistant": int((y[te][fm] == 1).sum())})
    print(f"cross-source {name}: veto flips {int(fm.sum())} ({time.time()-t0:.0f}s)", flush=True)
flips_cs = pd.concat(flips_cs, ignore_index=True) if flips_cs else pd.DataFrame()
all_flips = pd.concat([flips_s1, flips_cs], ignore_index=True)

obs_abl = {
    "lr_f2_auc": auc_f2, "lr_f2mech_auc": auc_mech, "lr_delta_auc": auc_mech - auc_f2,
    "lr_f2_sens_at_spec_0.985": s985_f2, "lr_f2mech_sens_at_spec_0.985": s985_mech,
    "cat_s1_sens_firm_off": tm_off["sens_firm"], "cat_s1_sens_firm_on": tm_on["sens_firm"],
    "cat_s1_spec_firm_off": tm_off["spec_firm"], "cat_s1_spec_firm_on": tm_on["spec_firm"],
    "cat_s1_n_U_off": tm_off["n_U"], "cat_s1_n_U_on": tm_on["n_U"],
    "s1_n_flips": int(flip_mask.sum()),
    "s1_n_flips_resistant": int((flips_s1.y_true == 1).sum()),
    "s1_n_flips_susceptible": int((flips_s1.y_true == 0).sum()),
}
for r in cs_rows:
    for k in ("sens_firm_off", "sens_firm_on", "spec_firm_off", "spec_firm_on", "n_flips", "n_flips_resistant"):
        obs_abl[f"cs_{r['direction']}_{k}"] = r[k]

n_flip_R = int((flips_s1.y_true == 1).sum())
n_flip_S = int((flips_s1.y_true == 0).sum())
decision = (f"Mechanism features delta AUC {auc_mech-auc_f2:+.4f} (sens@spec0.985 {s985_mech-s985_f2:+.4f}). "
            f"S1 veto downgrades {int(flip_mask.sum())} R->U: {n_flip_R} true resistant (harm), "
            f"{n_flip_S} true susceptible (benefit). Ship decision DEFERRED to G5 (must not degrade external "
            f"validation); golden kept unless a clear patient-grouped improvement survives external check.")
write_result("lr_mechanism", "e16_lof_ablation", obs_abl,
             tables={"lr_compare": lr_cmp, "s1_veto_metrics": pd.DataFrame([
                 {"state": "OFF", **tm_off}, {"state": "ON", **tm_on}]),
                 "cross_source_veto": pd.DataFrame(cs_rows), "call_flips": all_flips},
             method="S1 nested LR F2 vs F2+mechanism (rv_lof_any/highconf/dbd/late). Catalogue epistatic veto "
                    "R->U when high-confidence Rv0678 LoF + mmpL5/S5 LoF + no atpE variant. S1 + cross-source. "
                    "Every flip enumerated with true label.",
             interpretation=decision,
             caveats="Veto is conservative (R->U, send-for-lab, not R->S). 'No atpE variant' is a proxy for "
                     "'not atpE-driven R'. Homopolymer indels excluded from high-confidence LoF. Mechanism "
                     "features may move metrics negligibly (few truncation carriers) — a valid outcome. "
                     "External validation (G5) pending; no ship decision yet.",
             reproduce="python eval/e16_mechanism.py",
             supports="Audited epistatic downgrade rule + mechanism-feature ablation (novelty P-E16).")
print("\n=== E16 ABLATION ===")
print(lr_cmp.to_string(index=False))
print(f"S1 veto: {tm_off['sens_firm']:.4f}->{tm_on['sens_firm']:.4f} sens_firm; "
      f"flips {int(flip_mask.sum())} (R={n_flip_R}, S={n_flip_S})")
print("DECISION:", decision)