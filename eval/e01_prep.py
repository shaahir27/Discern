"""E01: data prep exact counts (AGENTS.md section 5)."""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import pandas as pd
import yaml

from discern.prep import E01_EXPECTED, e01_counts, prepare
from discern.results_io import write_result

cfg = yaml.safe_load(open("config.yaml"))
raw, proc = Path(cfg["paths"]["raw"]), Path(cfg["paths"]["processed"])
s, m = prepare(raw, proc)
obs = e01_counts(s, m)
tol = {k: 0 for k in E01_EXPECTED}  # exact counts
out = write_result(
    "_data_prep", "e01_prep", obs,
    expected=E01_EXPECTED, tolerance=tol,
    tables={"counts": pd.DataFrame([obs])},
    method="AGENTS.md §5: CRyPTIC QUALITY=HIGH + all NICD; dedup prefer R then MIC; inner-join GENOMES; "
           "mutations GENE@MUT with MINOR_MUTATION when minor, FRSf; drop synonymous/null.",
    interpretation="Exact reproduction of the published cohort counts (11,868 samples / 866 resistant).",
    caveats="median_FRSf_minor (~0.375) is approximate, not an exact acceptance row. "
            "Repo DST files hold 22,347 R/S IDs; the QUALITY=HIGH + genome join yields 11,868.",
    reproduce="python eval/make_manifest.py && python eval/e01_prep.py",
    supports="Baseline reproduction (not a novelty claim).",
)
import json
rec = json.loads(out.read_text())
print(json.dumps({"observed": obs, "pass_all": rec["pass_all"], "diff": rec["diff"]}, indent=1))