"""G3 / P-E14: conformal certificate audit (baseline vs certified policy) + M3 shrinkage.

Audits S1 pooled OOF, per-fold fluctuation, out-of-site, and cross-source certificate transfer;
prevalence-weighted empirical errors. No acceptance tolerance (new science; report violations).
Also produces the beta-binomial shrinkage estimator output (descriptive; ablate, don't headline).
"""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import json
import numpy as np
import pandas as pd
import yaml

from discern import conformal as C
from discern.shrinkage import variant_shrinkage
from discern.results_io import write_result

cfg = yaml.safe_load(open("config.yaml"))
proc = Path(cfg["paths"]["processed"])
models_root = Path(cfg["paths"]["models"])
s = pd.read_parquet(proc / "samples.parquet")
m = pd.read_parquet(proc / "mutations.parquet").merge(s[["UNIQUEID", "y"]], on="UNIQUEID", how="left")
y, ids = s.y.values, s.UNIQUEID.values
P = np.load(proc / "e04_oof_p.npy")
folds = np.load(proc / "s1_folds.npy")
SRC, SITE = s.SRC.values, s.SITEID.values
pi_tr = float(y.mean())
ALPHA, BETA, DELTA = 0.015, 0.20, 0.05

e04 = json.loads(Path("results/models/lr_model/e04_s1_patient_cv/results.json").read_text())
th = e04["extra"]["thresholds_per_fold"]
lo = np.array([th[int(f)]["lo"] for f in folds])
hi = np.array([th[int(f)]["hi"] for f in folds])
T_base = np.where(P <= lo, "S", np.where(P >= hi, "R", "U"))
tS, tR, info = C.certified_thresholds(P, y, ALPHA, BETA, DELTA)
T_cert = C.apply_policy(P, tS, tR)


def policy_summary(T):
    es, ns = C.zone_error(y, T, "S")
    er, nr = C.zone_error(y, T, "R")
    return {"emp_s_error": es, "emp_r_error": er, "n_S": ns, "n_R": nr,
            "n_U": int((T == "U").sum()), "n_U_resistant": int(((T == "U") & (y == 1)).sum())}


base, cert = policy_summary(T_base), policy_summary(T_cert)
obs = {
    "alpha": ALPHA, "beta": BETA, "delta": DELTA, "pi_train": pi_tr,
    **{f"baseline_{k}": v for k, v in base.items()},
    "certified_t_S": tS, "certified_t_R": tR,
    "certified_bound_s": info["bound_s"], "certified_bound_r": info["bound_r"],
    "certified_s_certified": info["s_certified"], "certified_r_certified": info["r_certified"],
    **{f"certified_{k}": v for k, v in cert.items()},
}

# per-fold fluctuation: certify on other folds, apply to held-out fold
fT, fS, fR = [], [], []
for k in range(5):
    cal, te = folds != k, folds == k
    ts, tr, _ = C.certified_thresholds(P[cal], y[cal], ALPHA, BETA, DELTA)
    Tk = C.apply_policy(P[te], ts, tr)
    es, _ = C.zone_error(y[te], Tk, "S")
    er, _ = C.zone_error(y[te], Tk, "R")
    fT.append((ts, tr)); fS.append(es); fR.append(er)
fS = [v for v in fS if v is not None]; fR = [v for v in fR if v is not None]
tS_vals = [t[0] for t in fT if t[0] is not None]
obs.update({"fold_tS_min": min(tS_vals) if tS_vals else None,
            "fold_tS_max": max(tS_vals) if tS_vals else None,
            "fold_n_tS_certified": len(tS_vals),
            "fold_emp_s_min": min(fS) if fS else None, "fold_emp_s_max": max(fS) if fS else None,
            "fold_emp_r_min": min(fR) if fR else None, "fold_emp_r_max": max(fR) if fR else None})

# best achievable S-zone bound (quantifies how far the strict alpha certificate is from reachable)
best_s_bound, best_s_t = 1.0, None
for t in np.unique(np.round(P, 4)):
    n, k = C.zone_counts(y, P, t, "S")
    if n >= 20:
        b = C.certified_upper_bound(k, n, DELTA)
        if b < best_s_bound:
            best_s_bound, best_s_t = b, float(t)
obs["s_best_achievable_bound"] = best_s_bound
obs["s_best_achievable_t"] = best_s_t

# out-of-site audit: pooled certificate applied per site
site_rows, viol_s, viol_r = [], 0, 0
for site in sorted(np.unique(SITE)):
    sel = SITE == site
    es, ns = C.zone_error(y[sel], T_cert[sel], "S")
    er, nr = C.zone_error(y[sel], T_cert[sel], "R")
    vs = es is not None and es > info["bound_s"]
    vr = er is not None and er > info["bound_r"]
    viol_s += int(vs); viol_r += int(vr)
    site_rows.append({"site": int(site), "n": int(sel.sum()), "emp_s_error": es, "emp_r_error": er,
                      "violates_s": vs, "violates_r": vr})
obs.update({"out_of_site_n_sites": len(site_rows), "out_of_site_violations_s": viol_s,
            "out_of_site_violations_r": viol_r})

# cross-source certificate transfer: certify on A OOF, evaluate on B OOF
cs_rows = []
for a, b in (("CRyPTIC", "NICD"), ("NICD", "CRyPTIC")):
    ca, cb = SRC == a, SRC == b
    ts, tr, inf = C.certified_thresholds(P[ca], y[ca], ALPHA, BETA, DELTA)
    Tb = C.apply_policy(P[cb], ts, tr)
    es, ns = C.zone_error(y[cb], Tb, "S")
    er, nr = C.zone_error(y[cb], Tb, "R")
    vs = es is not None and es > inf["bound_s"]
    vr = er is not None and er > inf["bound_r"]
    cs_rows.append({"direction": f"{a}_to_{b}", "bound_s": inf["bound_s"], "bound_r": inf["bound_r"],
                    "emp_s_error": es, "emp_r_error": er, "violates_s": vs, "violates_r": vr,
                    "n_U": int((Tb == "U").sum())})
    obs.update({f"cs_{a}_to_{b}_bound_s": inf["bound_s"], f"cs_{a}_to_{b}_bound_r": inf["bound_r"],
                f"cs_{a}_to_{b}_emp_s_error": es, f"cs_{a}_to_{b}_emp_r_error": er,
                f"cs_{a}_to_{b}_violates_s": vs, f"cs_{a}_to_{b}_violates_r": vr})

# prevalence-weighted empirical errors (approximate certificate at user prevalence)
for pi in (0.02, 0.05):
    tag = str(pi).replace(".", "")
    ws = C.weighted_zone_error(y, P, tS, "S", pi_tr, pi) if tS is not None else None
    wr = C.weighted_zone_error(y, P, tR, "R", pi_tr, pi) if tR < 2.0 else None
    obs[f"weighted_emp_s_error_pi{tag}"] = ws
    obs[f"weighted_emp_r_error_pi{tag}"] = wr

write_result("conformal", "e14_certificate_audit", obs,
             tables={"policy_compare": pd.DataFrame([{"policy": "baseline", **base}, {"policy": "certified", **cert}]),
                     "per_fold": pd.DataFrame([{"fold": k, "tS": fT[k][0], "tR": fT[k][1]} for k in range(5)]),
                     "out_of_site": pd.DataFrame(site_rows), "cross_source": pd.DataFrame(cs_rows)},
             method="LTT/Bentkus risk control on inner OOF calibrated scores. S-zone error<=alpha=0.015, "
                    "R-zone error<=beta=0.20 at delta=0.05; certified upper bound = inverted Bentkus. "
                    "Baseline = empirical NPV/PPV thresholds. Audited S1 pooled, per-fold, out-of-site, "
                    "cross-source; class-constant prevalence weights (Tibshirani 2019).",
             interpretation=f"R-zone certified: bound {info['bound_r']:.4f} <= beta {BETA} "
                            f"(achieved {cert['emp_r_error']}). S-zone NOT certifiable at alpha={ALPHA}/"
                            f"delta={DELTA}: best achievable bound {best_s_bound:.4f} (at tS={best_s_t}); the "
                            f"strict NPV>=98.5% certificate fails given the score distribution. Certified policy "
                            f"abstains from S (tS={tS}), U-zone {base['n_U']}->{cert['n_U']} (widens, as "
                            f"pre-registered). Out-of-site S violations {viol_s}/{len(site_rows)}; cross-source "
                            f"transfer in tables. The bound, not the point target, is the guarantee.",
             caveats="Patient grouping => effective n ~ patients; marginal split-conformal validity only, not "
                     "i.i.d. R-zone bound is wide (~170 positives/fold). Prevalence-weighted certificate is "
                     "approximate at an estimated user prevalence (Barber 2022 coverage-gap). Exchangeability "
                     "breaks across labs — violations are a RESULT, motivating the per-lab prevalence slider.",
             reproduce="python eval/e14_conformal.py",
             supports="First risk-controlled triage application to BDQ resistance (to our knowledge); "
                      "distribution-free abstention certificates re-weighted to lab prevalence.")
print(f"certified tS={tS} tR={tR} bound_s={info['bound_s']:.4f} bound_r={info['bound_r']:.4f} "
      f"R_certified={info['r_certified']}; U {base['n_U']}->{cert['n_U']}; "
      f"out-of-site viol S={viol_s}/{len(site_rows)}", flush=True)

# ---------- M3 shrinkage (descriptive) ----------
eff = pd.read_csv(models_root / "lr_mechanism" / "lof_effects.csv")
sh = variant_shrinkage(m, eff)
sdir = models_root / "lr_shrinkage"
sdir.mkdir(parents=True, exist_ok=True)
sh.to_csv(sdir / "variant_shrinkage.csv", index=False)
obs_sh = {
    "n_variants": int(len(sh)),
    "mean_abs_shrinkage": float((sh.shrunk_prop - sh.raw_prop).abs().mean()),
    "max_abs_shrinkage": float((sh.shrunk_prop - sh.raw_prop).abs().max()),
    "n_grade_flips_0.5": int(sh["grade_flip_0.5"].sum()),
    "n_mech_classes": int(sh.mech_class.nunique()),
}
write_result("lr_shrinkage", "e14_shrinkage", obs_sh,
             tables={"variant_shrinkage": sh.head(200),
                     "class_priors": sh.groupby("mech_class")[["prior_a", "prior_b"]].first().reset_index()},
             method="Empirical-Bayes beta-binomial shrinkage of per-variant R proportion toward mechanism-class "
                    "mean (method-of-moments Beta prior). Wilson 95% CI on raw proportion.",
             interpretation=f"{int(len(sh))} variants shrunk; mean |shrunk-raw| {obs_sh['mean_abs_shrinkage']:.4f}; "
                            f"{obs_sh['n_grade_flips_0.5']} variants flip a 0.5 grade. Rare variants pulled toward "
                            f"class mean; standard estimator choice, not an invention.",
             caveats="Descriptive estimator output for evidence cards; ablate, don't headline. No accuracy claim. "
                     "Class priors estimated from the same data (empirical Bayes).",
             reproduce="python eval/e14_conformal.py",
             supports="M3 shrinkage rung (beta-binomial mechanism-class evidence grading).")
print(f"shrinkage: {int(len(sh))} variants, mean|delta|={obs_sh['mean_abs_shrinkage']:.4f}, "
      f"grade_flips={obs_sh['n_grade_flips_0.5']}", flush=True)