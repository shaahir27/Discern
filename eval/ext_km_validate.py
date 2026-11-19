"""G5 / P-EXT: external validation on KMVALIDATION2025 genotypes (CRyPTIC v3.4.0 MUTATIONS).

Applies the PRODUCTION catalogue + model (trained on all 11,868) to the 792 non-overlapping KM
isolates' DNA variants; compares to their lab phenotype (35 R / 757 S). Sensitivity analysis
(assay/CC undocumented, MEDIUM quality); stratified by lineage. No tuning; report what comes out.
Requires data/external_cryptic_v34/MUTATIONS.parquet (Zenodo 15680920, CC-BY 4.0).
"""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import json

import joblib
import numpy as np
import pandas as pd
import piezo
import yaml
from km_cohort import load_km_cohort
from sklearn.metrics import roc_auc_score

from discern.catalogue import apply_catalogue
from discern.features import design, sample_tables
from discern.metrics import ternary_metrics
from discern.model import score, sens_at_spec
from discern.results_io import write_result

cfg = yaml.safe_load(Path("config.yaml").read_text(encoding="utf-8"))
proc = Path(cfg["paths"]["processed"])
ext = Path("data/external_cryptic_v34")
models_dir = Path(cfg["paths"]["models"]) / "lr_model"
km, m = load_km_cohort(ext, proc)

ids = km.UNIQUEID.values
y = km.y.values
s_ext = km[["UNIQUEID", "y"]].copy()

# ---- production catalogue + model ----
rc = piezo.ResistanceCatalogue(str(models_dir / "catalogue_bdq.csv"))
calls = apply_catalogue(rc, s_ext, m, frs_min=cfg["features"]["frs_min"])
tm = ternary_metrics(y, calls)

art = joblib.load(models_dir / "bdq_lr.joblib")
core, a = sample_tables(m, ids, art["frs_min"], tuple(art["dbd"]))
X = design(core, a, art["vocab"], "F2").values
raw_p, P = score(art["final"], art["platt"], X)
tb = art["thresholds"]["t_bin"]
pb = P >= tb
nR, nS = max((y == 1).sum(), 1), max((y == 0).sum(), 1)
auc = float(roc_auc_score(y, P)) if len(np.unique(y)) > 1 else None

obs = {
    "n_isolates": len(km), "n_resistant": int(y.sum()), "zero_training_overlap": True,
    "n_with_variants": int(m.UNIQUEID.nunique()),
    "cat_firm": tm["firm"], "cat_sens_firm": tm["sens_firm"], "cat_spec_firm": tm["spec_firm"],
    "cat_sens_U_as_S": tm["sens_U_as_S"], "cat_spec_U_as_S": tm["spec_U_as_S"],
    "cat_sens_U_as_R": tm["sens_U_as_R"], "cat_spec_U_as_R": tm["spec_U_as_R"],
    "cat_n_U": tm["n_U"], "cat_n_U_resistant": tm["n_U_resistant"],
    "model_auc": auc,
    "model_binary_sens": float((pb & (y == 1)).sum() / nR),
    "model_binary_spec": float(((~pb) & (y == 0)).sum() / nS),
    "model_sens_at_spec_0.95": sens_at_spec(y, P, 0.95) if auc else None,
    "model_sens_at_spec_0.97": sens_at_spec(y, P, 0.97) if auc else None,
}
# lineage stratification (catalogue ternary)
lin_rows = []
for lin, grp in km.groupby("lineage"):
    sel = km.lineage.values == lin
    if (y[sel] == 1).sum() == 0:
        continue
    tml = ternary_metrics(y[sel], calls[sel])
    lin_rows.append({"lineage": lin, "n": int(sel.sum()), "n_resistant": int((y[sel] == 1).sum()),
                     "sens_firm": tml["sens_firm"], "spec_firm": tml["spec_firm"], "n_U": tml["n_U"]})
    obs[f"lin_{lin}_sens_firm"] = tml["sens_firm"]
    obs[f"lin_{lin}_spec_firm"] = tml["spec_firm"]

write_result("external_km", "e_km_validation", obs,
             tables={"metrics": pd.DataFrame([obs]), "lineage_stratification": pd.DataFrame(lin_rows)},
             method="Production catalogue (built on all 11,868) + production calibrated LR applied to "
                    "KMVALIDATION2025 DNA variants from CRyPTIC v3.4.0 MUTATIONS.parquet (Zenodo 15680920, "
                    "CC-BY 4.0). 792 isolates, ZERO training overlap. Compared to lab phenotype (35 R/757 S).",
             interpretation=f"EXTERNAL test on {len(km)} non-overlapping isolates: catalogue sens_firm "
                            f"{tm['sens_firm']:.4f} / spec_firm {tm['spec_firm']:.4f}; model AUC {auc:.4f}, "
                            f"binary sens {obs['model_binary_sens']:.4f} / spec {obs['model_binary_spec']:.4f}. "
                            f"Stratified by lineage in tables.",
             caveats="SENSITIVITY ANALYSIS, not headline proof: KM assay/CC UNDOCUMENTED (MEDIUM quality). "
                     "Only 35 resistant -> wide sensitivity CI. No Indian external resistant isolates exist. "
                     "Production model has NOT seen these isolates (zero overlap verified). match != beat.",
             reproduce="python eval/ext_km_validate.py",
             supports="First external evaluation of a BDQ triage tool on non-overlapping isolates (KM cohort).")
print(json.dumps({"n": len(km), "nR": int(y.sum()), "cat_sens_firm": tm["sens_firm"],
                  "cat_spec_firm": tm["spec_firm"], "model_auc": auc,
                  "model_sens": obs["model_binary_sens"], "model_spec": obs["model_binary_spec"]}, indent=1))
