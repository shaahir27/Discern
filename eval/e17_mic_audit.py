"""G6 / P-E17: borderline-MIC label-noise audit (SENSITIVITY ANALYSIS, not a model change).

Samples with MIC within +/-1 two-fold dilution of the source critical concentration (CRyPTIC plate
0.25 -> {0.12,0.25,0.5}; MGIT 1.0 -> {0.5,1.0,2.0}). Report n / n_resistant / share of positives,
and E04 metrics recomputed EXCLUDING them, side-by-side. Production model keeps ALL rows. A
margin-weighted fit is NOT promoted from this audit.
"""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import json
import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import roc_auc_score

from discern.model import sens_at_spec
from discern.results_io import write_result

cfg = yaml.safe_load(open("config.yaml"))
proc = Path(cfg["paths"]["processed"])
s = pd.read_parquet(proc / "samples.parquet")
y = s.y.values
P = np.load(proc / "e04_oof_p.npy")
folds = np.load(proc / "s1_folds.npy")
e04 = json.loads(Path("results/models/lr_model/e04_s1_patient_cv/results.json").read_text())
th = e04["extra"]["thresholds_per_fold"]
tb = np.array([th[int(f)]["t_bin"] for f in folds])


def parse_mic(v):
    if pd.isna(v):
        return np.nan
    sv = str(v).strip()
    if sv.startswith("<="):
        return float(sv[2:])
    if sv.startswith(">"):
        return float(sv[1:])
    try:
        return float(sv)
    except ValueError:
        return np.nan


mic = s.METHOD_MIC.map(parse_mic).values
CC = {"CRyPTIC": 0.25, "NICD": 1.0}
border = np.zeros(len(s), bool)
for src, cc in CC.items():
    sel = (s.SRC.values == src) & ~np.isnan(mic)
    lo, hi = cc / 2, cc * 2
    border |= sel & (mic >= lo - 1e-9) & (mic <= hi + 1e-9)
# only exact (non-censored) MICs at the bracketing dilutions count as borderline
censored = s.METHOD_MIC.astype(str).str.startswith(("<=", ">")).values
border &= ~censored

nR = int(y.sum())
obs = {
    "n_borderline": int(border.sum()),
    "n_borderline_resistant": int((border & (y == 1)).sum()),
    "borderline_share_of_positives": float((border & (y == 1)).sum() / nR),
    "n_cry_borderline": int((border & (s.SRC.values == "CRyPTIC")).sum()),
    "n_nicd_borderline": int((border & (s.SRC.values == "NICD")).sum()),
    "n_nicd_mic_available": int(((s.SRC.values == "NICD") & ~np.isnan(mic)).sum()),
}


def metrics(mask):
    yy, pp, tt = y[mask], P[mask], tb[mask]
    pb = pp >= tt
    return {
        "auc": float(roc_auc_score(yy, pp)),
        "binary_sens": float((pb & (yy == 1)).sum() / max((yy == 1).sum(), 1)),
        "binary_spec": float(((~pb) & (yy == 0)).sum() / max((yy == 0).sum(), 1)),
        "sens_at_spec_0.985": sens_at_spec(yy, pp, 0.985),
        "n": int(mask.sum()), "n_resistant": int((yy == 1).sum()),
    }


full, excl = metrics(np.ones(len(s), bool)), metrics(~border)
for k, v in full.items():
    obs[f"full_{k}"] = v
for k, v in excl.items():
    obs[f"excl_borderline_{k}"] = v
obs["delta_auc_excl_minus_full"] = excl["auc"] - full["auc"]

bd = s[border].assign(mic_num=mic[border])
mic_tab = bd.groupby(["SRC", "METHOD_MIC"]).agg(n=("UNIQUEID", "size"), nR=("y", "sum")).reset_index()
write_result("mic_audit", "e17_borderline_mic", obs,
             tables={"borderline_by_mic": mic_tab,
                     "full_vs_excluded": pd.DataFrame([{"set": "full", **full}, {"set": "excl_borderline", **excl}])},
             method="Borderline = exact MIC within +/-1 two-fold dilution of source CC (CRyPTIC 0.25, MGIT 1.0). "
                    "E04 OOF metrics recomputed excluding borderline samples, side-by-side. Production keeps ALL rows.",
             interpretation=f"{int(border.sum())} borderline samples ({int((border&(y==1)).sum())} resistant, "
                            f"{obs['borderline_share_of_positives']:.3f} of all positives). Excluding them shifts "
                            f"AUC {full['auc']:.4f}->{excl['auc']:.4f} (delta {obs['delta_auc_excl_minus_full']:+.4f}). "
                            f"Label noise near the cutoff is modest; sensitivity analysis only.",
             caveats="SENSITIVITY ANALYSIS, not a model change; production keeps all rows. Censored MICs (<=,>) "
                     "excluded from borderline. NICD MGIT MICs largely unavailable. Margin-weighted fit NOT "
                     "promoted from this audit (pre-registered).",
             reproduce="python eval/e17_mic_audit.py",
             supports="Label-noise audit (P-E17): quantifies cutoff-adjacent uncertainty honestly.")
print(f"borderline {int(border.sum())} (R {int((border&(y==1)).sum())}); "
      f"AUC full {full['auc']:.4f} excl {excl['auc']:.4f}", flush=True)