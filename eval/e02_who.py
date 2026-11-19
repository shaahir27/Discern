"""E02: WHO v2 catalogue applied as-is (AGENTS.md section 7)."""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import pandas as pd
import piezo
import yaml

from discern.catalogue import who_as_is
from discern.metrics import ternary_metrics
from discern.results_io import write_result

cfg = yaml.safe_load(open("config.yaml"))
raw, proc = Path(cfg["paths"]["raw"]), Path(cfg["paths"]["processed"])
s = pd.read_parquet(proc / "samples.parquet")
m = pd.read_parquet(proc / "mutations.parquet")
who = raw / "catalogues/WHOv2/NC_000962.3_WHO-UCN-TB-2023.5_v2.0_GARC1_RFUS.csv"
rc = piezo.ResistanceCatalogue(str(who))

EXP = {
    0.75: {"firm": 0.968, "sens_firm": 0.561, "spec_firm": 0.987, "sens_U_as_S": 0.497},
    0.1: {"firm": 0.961, "sens_firm": 0.801, "spec_firm": 0.981, "sens_U_as_S": 0.679},
}
for thr in (0.75, 0.1):
    calls, cache = who_as_is(rc, s, m, thr)
    obs = ternary_metrics(s.y, calls)
    obs["variant_call_counts"] = pd.Series(cache).value_counts().to_dict()
    tol = {k: 0.005 for k in EXP[thr]}
    out = write_result(
        "who_v2", f"e02_who_as_is_frs{thr}", obs,
        expected=EXP[thr], tolerance=tol,
        tables={"metrics": pd.DataFrame([obs])},
        method=f"WHO v2 GARC catalogue via piezo; non-syn/non-null variants FRSf>={thr}; "
               "exception->S; aggregate R>U>S.",
        interpretation=f"WHO v2 as-is at FRS>={thr}, pooled. Anchor sens_U_as_S must match published.",
        caveats="U counted as S for the anchor sensitivity. NICD enriches resistance; read per-source.",
        reproduce="python eval/e02_who.py",
        supports="Baseline reproduction (published WHO v2 performance).",
    )
    import json
    rec = json.loads(out.read_text())
    print(f"thr {thr}: sens_U_as_S={obs['sens_U_as_S']:.4f} pass={rec['pass_all']}")