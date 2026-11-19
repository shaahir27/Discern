"""Summarize uncertainty from saved evaluation predictions and counts."""

import sys
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import average_precision_score, roc_auc_score

sys.path.insert(0, "src")

from discern.features import CORE
from discern.metrics import ternary_metrics
from discern.research import BANNER, Reports
from discern.results_io import write_result

LEVEL = 0.95
REPLICATES = 500


def wilson(value: float, n: int) -> tuple[float, float]:
    p = round(value * n) / n
    z = NormalDist().inv_cdf((1 + LEVEL) / 2)
    denominator = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denominator
    radius = z * (p * (1 - p) / n + z**2 / (4 * n**2)) ** 0.5 / denominator
    return center - radius, center + radius


def main() -> None:
    cfg = yaml.safe_load(Path("config.yaml").read_text(encoding="utf-8"))
    proc = Path(cfg["paths"]["processed"])
    reports = Reports(Path("."))
    oof = pd.read_parquet(reports.root / "oof_predictions.parquet")
    patients = pd.read_parquet(proc / "samples.parquet", columns=["UNIQUEID", "SUBJID"]).set_index("UNIQUEID")
    groups = oof.join(patients, on="UNIQUEID").groupby("SUBJID", sort=False).indices
    assert sum(map(len, groups.values())) == len(oof)
    blocks = list(groups.values())
    y, calls, p = (oof[key].to_numpy() for key in ("y", "cat_call", "p"))
    def scores(indices: np.ndarray) -> tuple[float, ...]:
        metrics = ternary_metrics(y[indices], calls[indices])
        return (
            metrics["sens_firm"], metrics["spec_firm"],
            roc_auc_score(y[indices], p[indices]), average_precision_score(y[indices], p[indices]),
        )

    points = scores(np.arange(len(oof)))
    golden = reports.observed("catalogue"), reports.observed("model")
    assert np.allclose(points, [golden[0]["sens_firm"], golden[0]["spec_firm"], golden[1]["auc"], golden[1]["auprc"]])
    rng = np.random.default_rng(cfg["seeds"]["bootstrap"])
    draws = []
    for _ in range(REPLICATES):
        selected = rng.integers(0, len(blocks), size=len(blocks))
        draws.append(scores(np.concatenate([blocks[i] for i in selected])))
    limits = np.quantile(draws, [(1 - LEVEL) / 2, (1 + LEVEL) / 2], axis=0)
    rows = [
        {"metric": key, "estimate": v, "lower": float(limits[0, i]), "upper": float(limits[1, i]),
         "method": "S1 patient-group bootstrap of saved predictions"}
        for i, (key, v) in enumerate(zip(("S1 catalogue sensitivity", "S1 catalogue specificity", "S1 model AUC", "S1 model AUPRC"), points, strict=True))
    ]
    km = reports.observed("km")
    n_r = km["n_resistant"] - km["cat_n_U_resistant"]
    n_s = km["n_isolates"] - km["cat_n_U"] - n_r
    for label, key, n in (("sensitivity", "cat_sens_firm", n_r), ("specificity", "cat_spec_firm", n_s)):
        value = km[key]
        lo, hi = wilson(value, n)
        rows.append({"metric": f"KM catalogue {label}", "estimate": value, "lower": lo, "upper": hi,
                     "method": "Wilson interval among firm calls"})
    muts = pd.read_parquet(proc / "mutations.parquet", columns=["UNIQUEID", "GENE", "IS_SYN", "IS_NULL", "FRSf"])
    retained = muts[muts.GENE.isin(CORE) & ~muts.IS_SYN & ~muts.IS_NULL & (muts.FRSf >= cfg["features"]["frs_min"])]
    total_r = int((y == 1).sum())
    without = int(((y == 1) & ~oof.UNIQUEID.isin(retained.UNIQUEID)).sum())
    write_result(
        "research_archive", "g8_evaluation_context",
        {"interval_level": LEVEL, "bootstrap_replicates": REPLICATES, "n_resistant": total_r,
         "n_resistant_without_retained_core_variant": without, "fraction_without_retained_core_variant": without / total_r},
        tables={"intervals": pd.DataFrame(rows)}, extra={"intervals": rows},
        method="Patient resampling of the archived S1 predictions and Wilson intervals from saved KM firm-call totals.",
        caveats=f"{BANNER} Resampling fixed predictions measures sample uncertainty. Refitting variability is outside these intervals.",
        reproduce="python eval/g8_context.py",
    )


if __name__ == "__main__":
    main()
