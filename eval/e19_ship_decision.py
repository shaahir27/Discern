"""Apply the preregistered D011 ship rule using saved evidence and the KM audit."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import piezo
import yaml

sys.path.insert(0, "src")

from km_cohort import load_km_cohort

from discern import mechanism as M
from discern.catalogue import apply_catalogue
from discern.metrics import ternary_metrics
from discern.results_io import write_result

REPRODUCE = "C:/Users/shaah/miniconda3/envs/discern/python.exe eval/e19_ship_decision.py"
KM_CAVEAT = "Research use only. KM sensitivity analysis with 35 resistant isolates, undocumented assay and critical concentration, and medium-quality phenotypes."


def saved(results: Path, test: str) -> dict:
    return json.loads((results / test / "results.json").read_text(encoding="utf-8"))


def call_metrics(y: np.ndarray, off: np.ndarray, on: np.ndarray) -> dict:
    results = {}
    for state, calls in (("off", off), ("on", on)):
        metrics = ternary_metrics(y, calls)
        for metric in ("sens_firm", "spec_firm"):
            label = 1 if metric == "sens_firm" else 0
            results[f"cat_{metric}_{state}"] = metrics[metric] if ((y == label) & (calls != "U")).any() else None
        results[f"cat_n_U_{state}"] = metrics["n_U"]
    flips = off != on
    return {
        **results, "n_flips": int(flips.sum()),
        "n_flips_resistant": int((flips & (y == 1)).sum()),
        "n_flips_susceptible": int((flips & (y == 0)).sum()),
        "n_true_R_to_S": int(((y == 1) & (off == "R") & (on == "S")).sum()),
        "all_flips_R_to_U": bool(((off[flips] == "R") & (on[flips] == "U")).all()),
    }


def main() -> None:
    cfg = yaml.safe_load(Path("config.yaml").read_text(encoding="utf-8"))
    proc, results, models = (Path(cfg["paths"][key]) for key in ("processed", "results", "models"))
    golden = saved(results, "models/lr_model/e04_s1_patient_cv")["observed"]
    golden_km = saved(results, "approaches/external_km/e_km_validation")["observed"]
    e16 = saved(results, "models/lr_mechanism/e16_lof_ablation")["observed"]
    shrinkage = saved(results, "models/lr_shrinkage/e14_shrinkage")["observed"]
    conformal = saved(results, "approaches/conformal/e14_certificate_audit")["observed"]
    comparison = saved(results, "approaches/head_to_head/e21_comparison")
    assert not comparison["extra"]["qualifying"] and e16["lr_delta_auc"] <= 0.005

    km, muts = load_km_cohort(Path("data/external_cryptic_v34"), proc)
    ids, y = km.UNIQUEID.to_numpy(), km.y.to_numpy()
    assert (len(ids), int(y.sum())) == (golden_km["n_isolates"], golden_km["n_resistant"]) == (792, 35)
    assert km.PHENOTYPE.isin(["R", "S"]).all() and km.lineage.notna().all()
    frs_min = cfg["features"]["frs_min"]
    rc = piezo.ResistanceCatalogue(str(models / "lr_model/catalogue_bdq.csv"))
    cds = M.load_cds(proc / "rv0678_cds.fasta")
    _, flags = M.sample_mech_flags(muts, ids, cds, frs_min=frs_min)
    off = np.asarray(apply_catalogue(rc, km[["UNIQUEID", "y"]], muts, frs_min=frs_min))
    on = M.apply_veto(off, flags.veto_fire.to_numpy())
    obs = call_metrics(y, off, on)
    for key in ("sens_firm", "spec_firm", "n_U"):
        assert obs[f"cat_{key}_off"] == golden_km[f"cat_{key}"], f"Golden KM {key} drifted"
    obs.update(
        n_isolates=len(y), n_resistant=int(y.sum()), zero_training_overlap=True,
        model_auc=golden_km["model_auc"],
    )
    calls = km[["UNIQUEID", "lineage", "y"]].assign(call_off=off, call_on=on)
    flips = calls.loc[off != on]
    lineages = []
    for lineage, mask in ((lin, km.lineage.to_numpy() == lin) for lin in sorted(km.lineage.unique())):
        lineages.append({
            "lineage": lineage, "n": int(mask.sum()), "n_resistant": int(y[mask].sum()),
            **call_metrics(y[mask], off[mask], on[mask]),
        })
    assert len(flips) == sum(row["n_flips"] for row in lineages) == obs["n_flips"]
    assert obs["n_flips"] == obs["n_flips_resistant"] + obs["n_flips_susceptible"]
    obs["per_lineage"] = {row["lineage"]: {k: v for k, v in row.items() if k != "lineage"} for row in lineages}

    s1_reduction = e16["cat_s1_sens_firm_off"] - e16["cat_s1_sens_firm_on"]
    gates = {
        "km_specificity_preserved": obs["cat_spec_firm_on"] >= obs["cat_spec_firm_off"],
        "km_safety": obs["all_flips_R_to_U"] and obs["n_true_R_to_S"] == 0,
        "s1_sensitivity_reduction_within_margin": s1_reduction <= 0.005,
    }
    adopt = all(gates.values())
    pick = {"model_folder": "models/lr_model", "epistatic_veto": adopt}
    decision = "Adopt the audited Tier A veto modifier." if adopt else "Keep the golden catalogue."
    write_result(
        "ship_decision", "e19_km_veto", obs, results_root=results,
        tables={"call_flips": flips, "lineage_stratification": pd.DataFrame(lineages)},
        extra={"d011_checks": gates, "s1_sens_firm_reduction": s1_reduction, "final_pick": pick},
        method="Apply the frozen golden catalogue to the fixed KM cohort. Audit the existing veto on each call and use the saved production KM AUC as the LR anchor.",
        interpretation=decision, caveats=KM_CAVEAT, reproduce=REPRODUCE,
    )
    rows = [
        {
            "rung": "M0", "name": "Golden L2 LR (F2)",
            "S1": {
                "auc": golden["auc"], "sens_at_spec_0.985": golden["sens_at_spec_0.985"],
                "cat_sens_firm": e16["cat_s1_sens_firm_off"], "cat_spec_firm": e16["cat_s1_spec_firm_off"],
            },
            "KM": {key: golden_km[key] for key in ("model_auc", "cat_sens_firm", "cat_spec_firm")},
            "qualifies": True,
            "reason": "E21 retained the golden model class under D013.",
        },
        {
            "rung": "M1", "name": "Mechanism features",
            "S1": {key: value for key, value in e16.items() if key.startswith("lr_")}, "KM": None,
            "qualifies": False, "reason": "The S1 AUC gain was below the predeclared margin.",
        },
        {
            "rung": "M2", "name": "Epistatic veto",
            "S1": {key: value for key, value in e16.items() if not key.startswith("lr_")},
            "KM": {key: value for key, value in obs.items() if key != "per_lineage"},
            "qualifies": adopt, "reason": decision,
        },
        {
            "rung": "M3", "name": "Shrinkage", "S1": None, "KM": None, "qualifies": False,
            "reason": "Descriptive evidence-card estimator.",
            "evidence": shrinkage,
        },
        {
            "rung": "M4", "name": "Conformal abstention", "S1": conformal, "KM": None, "qualifies": False,
            "reason": "The R-zone certified. The S-zone missed its target and certificates failed under lab shift.",
        },
    ]
    for row in rows:
        row["final_pick"] = row["rung"] == "M0" or (row["rung"] == "M2" and adopt)
    table = pd.DataFrame(rows)
    for col in ("S1", "KM", "evidence"):
        table[col] = table[col].map(lambda value: json.dumps(value, sort_keys=True) if isinstance(value, dict) else "")
    write_result(
        "ship_decision", "e19_ladder_summary", {"model": pick["model_folder"], "epistatic_veto": adopt},
        results_root=results, tables={"ladder": table},
        extra={"ladder": rows, "final_pick": pick, "d011_checks": gates},
        method="Apply D011 to the saved E16, shrinkage, conformal, and E21 records and the KM veto audit.",
        interpretation=f"Keep golden L2 LR (F2). {decision}",
        caveats=f"{KM_CAVEAT} Report the certificate bounds with their confidence level and lab-shift results.",
        reproduce=REPRODUCE,
    )
    print(json.dumps({key: value for key, value in obs.items() if key != "per_lineage"}, indent=2))
    print(f"D011 checks {gates}. {decision}")


if __name__ == "__main__":
    main()
