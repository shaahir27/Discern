"""Structured results writer.

Layout (user-aligned): two sibling containers inside results/
    results/
      INDEX.md / index.json
      models/                     # trained models only (type=model)
        <model_name>/
          MODEL_REPORT.md         # consolidated report
          <test>/
            REPORT.md  results.json  *.csv
      approaches/                 # everything else (baseline/analysis/data/method/...)
        <approach>/
          SUMMARY.md              # consolidated report
          <test>/
            REPORT.md  results.json  *.csv

Routing is by approach type: type=model -> models/, otherwise -> approaches/.
results.json is the source of truth (incl. prose). Markdown is a rounded, UTF-8 rendered view.
Summaries and INDEX regenerate on every run; regenerate_all() re-renders without re-running
experiments. migrate() re-keys legacy names and places each folder in the correct container.
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

import pandas as pd

# Canonical approaches: name -> (display name, type). type drives grouping + summary filename.
APPROACHES = {
    "data_prep": ("Data preparation", "data"),
    "who_baseline": ("WHO v2 baseline", "baseline"),
    "catalogue": ("Tier A — catalogue", "baseline"),
    "lr_model": ("Tier B — model (LR F2)", "model"),
    "lr_expanded": ("Tier B — expanded-training (CRyPTIC v3.4)", "model"),
    "lr_l1": ("Tier B — L1 logistic regression", "model"),
    "rf_model": ("Tier B — random forest", "model"),
    "hgb_model": ("Tier B — histogram gradient boosting", "model"),
    "svm_model": ("Tier B — linear SVM", "model"),
    "lr_l2_tuned": ("Tier B — tuned L2 logistic regression", "model"),
    "lr_l1_tuned": ("Tier B — tuned L1 logistic regression", "model"),
    "rf_tuned": ("Tier B — tuned random forest", "model"),
    "hgb_tuned": ("Tier B — tuned histogram gradient boosting", "model"),
    "svm_tuned": ("Tier B — tuned SVM", "model"),
    "structure_test": ("Structure novel-variant test", "model"),
    "hybrid": ("Hybrid (exploratory)", "analysis"),
    "prevalence_correction": ("Prevalence correction", "analysis"),
    "priority_list": ("Variant priority list", "analysis"),
    "predictive_values": ("Predictive values (PPV/NPV)", "analysis"),
    "lr_mechanism": ("Tier B — mechanism-enhanced", "model"),
    "lr_shrinkage": ("Tier B — shrinkage estimator", "model"),
    "conformal": ("Conformal certificates", "method"),
    "voi_queue": ("Value-of-information queue", "analysis"),
    "external_km": ("External — KM validation", "validation"),
    "external_sobol": ("External — Sobol validation", "validation"),
    "external_choi": ("External — Choi validation", "validation"),
    "mic_audit": ("MIC label-noise audit", "analysis"),
    "head_to_head": ("Head-to-head comparisons", "comparison"),
    "ship_decision": ("Production ship decision", "analysis"),
}

# Legacy folder names -> canonical approach.
ALIAS = {
    "_data_prep": "data_prep", "who_v2": "who_baseline", "catalogue_tierA": "catalogue",
    "lr_f2_baseline": "lr_model", "structure_novel_variant": "structure_test",
    "hybrid_exploratory": "hybrid", "lr_mechanism_enhanced": "lr_mechanism",
    "conformal_certified": "conformal", "external_validation": "external_km",
}

TYPE_ORDER = ["model", "baseline", "method", "analysis", "validation", "comparison", "data"]
TYPE_LABELS = {"model": "Models", "baseline": "Baselines", "method": "Methods",
               "analysis": "Analyses", "validation": "External validation",
               "comparison": "Comparisons", "data": "Data"}

TEST_TITLES = {
    "e01_prep": "E01 · Data preparation",
    "e02_who_as_is_frs0.1": "E02 · WHO v2 as-is (FRS ≥ 0.1)",
    "e02_who_as_is_frs0.75": "E02 · WHO v2 as-is (FRS ≥ 0.75)",
    "e03_s1_catalogue_cv": "E03 · Catalogue CV (S1, patient-grouped)",
    "e04_s1_patient_cv": "E04 · Model CV (S1, patient-grouped)",
    "e05_ablation": "E05 · Feature ablation",
    "e06_cross_source": "E06 · Cross-source (S3)",
    "e07_loso": "E07 · Leave-one-site-out (S2)",
    "e08_novel_variant": "E08 · Novel-variant holdout (S4)",
    "e09_hybrid": "E09 · Hybrid on catalogue-U",
    "e10_prevalence_shift": "E10 · Calibration with known target prevalence",
    "e15_label_reveal_sim": "E15 · Retrospective comparison of variant sets",
    "e11_imbalance_learning_curve": "E11 · Imbalance + learning curve",
    "e12_variant_priority": "E12 · Variant priority list",
    "e13_ppv_npv": "E13 · PPV/NPV at prevalence",
}
_ACRONYMS = {"cv": "CV", "loso": "LOSO", "ppv": "PPV", "npv": "NPV", "frs": "FRS", "who": "WHO",
             "mic": "MIC", "voi": "VOI", "s1": "S1", "s2": "S2", "s3": "S3", "s4": "S4",
             "lr": "LR", "f2": "F2", "auc": "AUC", "km": "KM"}


def _approach(model: str) -> str:
    return ALIAS.get(model, model)


APPROACHES_DIR = "approaches"
MODELS_DIR = "models"


def _container(approach: str) -> str:
    return MODELS_DIR if _meta(approach)[1] == "model" else APPROACHES_DIR


def _approaches_root(results_root: Path) -> Path:
    return results_root / APPROACHES_DIR


def _models_root(results_root: Path) -> Path:
    return results_root / MODELS_DIR


def _approach_dir(results_root: Path, approach: str) -> Path:
    return results_root / _container(approach) / approach


def _container_roots(results_root: Path) -> list[Path]:
    return [_models_root(results_root), _approaches_root(results_root)]


def _meta(approach: str) -> tuple[str, str]:
    return APPROACHES.get(approach, (approach.replace("_", " ").title(), "analysis"))


def _summary_filename(approach: str) -> str:
    return "MODEL_REPORT.md" if _meta(approach)[1] == "model" else "SUMMARY.md"


def _versions() -> dict:
    out = {}
    for pkg in ("scikit-learn", "pandas", "numpy", "scipy", "catomatic", "piezo", "pyarrow"):
        try:
            out[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            out[pkg] = None
    return out


def _git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return "uncommitted"


def _manifest_hash(results_root: Path) -> str | None:
    mf = results_root.parent / "data/raw/MANIFEST.json"
    if mf.exists():
        try:
            return json.loads(mf.read_text(encoding="utf-8")).get("manifest_sha256")
        except Exception:
            return None
    return None


def _is_num(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _fmt(v) -> str:
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return f"{v:.4f}"
    return str(v)


def _badge(p) -> str:
    return "PASS ✅" if p is True else ("FAIL ❌" if p is False else "DESCRIPTIVE •")


def _mark(p) -> str:
    return "✅" if p is True else ("❌" if p is False else "—")


def _test_title(test: str) -> str:
    if test in TEST_TITLES:
        return TEST_TITLES[test]
    m = re.match(r"(e\d+|ext|g\d+)_(.*)", test)
    tag, rest = (m.group(1).upper(), m.group(2)) if m else ("", test)
    words = [_ACRONYMS.get(w.lower(), w[:1].upper() + w[1:]) for w in rest.replace("_", " ").split()]
    title = " ".join(words)
    return f"{tag} · {title}" if tag else title


def _exp_order(test: str) -> tuple:
    m = re.match(r"e(\d+)", test)
    return (0, int(m.group(1))) if m else (1, test)


def _headline(observed: dict) -> str:
    for k in ("auc", "sens_firm", "pooled_sens", "sens_U_as_S", "n_samples", "auc_plus_structure",
              "auc_ml_within_U", "n_distinct_sets"):
        if k in observed and _is_num(observed[k]):
            return f"{k}={_fmt(observed[k])}"
    nums = [(k, v) for k, v in observed.items() if _is_num(v)]
    return f"{nums[0][0]}={_fmt(nums[0][1])}" if nums else ""


def _metrics_tables(o, e, t, w) -> list[str]:
    acc = [k for k in o if k in e and _is_num(o[k])]
    desc = [k for k in o if k not in e and _is_num(o[k])]
    details = [k for k in o if not _is_num(o[k])]
    lines = []
    if acc:
        lines += ["| Metric | Value | Expected | ±Tol | Status |", "|---|---|---|---|---|"]
        lines += [f"| `{k}` | {_fmt(o[k])} | {_fmt(e[k])} | {_fmt(t.get(k, ''))} | {_mark(w.get(k))} |" for k in acc]
        lines.append("")
    if desc:
        lines += ["| Metric | Value |", "|---|---|"]
        lines += [f"| `{k}` | {_fmt(o[k])} |" for k in desc]
        lines.append("")
    if details:
        lines += [f"- **{k}:** `{o[k]}`" for k in details]
        lines.append("")
    return lines


def render_report(record: dict) -> str:
    o, e, t, w = record["observed"], record.get("expected", {}), record.get("tolerance", {}), record.get("within_tol", {})
    lines = [
        f"# {_test_title(record['test'])}",
        "",
        f"> **Approach:** `{record['model']}` · **Test:** `{record['test']}` · **{_badge(record.get('pass_all'))}**",
        "",
    ]
    if record.get("interpretation"):
        lines += [record["interpretation"], ""]
    lines += ["## Metrics", ""] + _metrics_tables(o, e, t, w)
    if record.get("method"):
        lines += ["## Method", record["method"], ""]
    if record.get("supports"):
        lines += ["## Novelty claim", record["supports"], ""]
    if record.get("caveats"):
        lines += ["## Caveats / negative results", record["caveats"], ""]
    if record.get("reproduce"):
        lines += ["## Reproduce", "```bash", record["reproduce"], "```", ""]
    mf = record.get("data_manifest_sha256") or "n/a"
    ver = " · ".join(f"{k} {v}" for k, v in record.get("versions", {}).items() if v)
    lines += [
        "## Provenance",
        f"- **Date (UTC):** {record.get('timestamp', 'n/a')}",
        f"- **Git commit:** {record.get('git_commit', 'n/a')}",
        f"- **Data manifest:** `{mf[:12]}…`",
        f"- **Versions:** {ver}",
        "",
        "## Files",
        "- `results.json` (source of truth)",
    ]
    lines += [f"- `{f}`" for f in record.get("csv_files", [])]
    return "\n".join(lines) + "\n"


def render_approach_summary(approach: str, records: list[dict]) -> str:
    name, typ = _meta(approach)
    records = sorted(records, key=lambda r: _exp_order(r["test"]))
    npass = sum(1 for r in records if r.get("pass_all") is True)
    lines = [
        f"# {name}",
        "",
        f"> **Type:** {typ} · **Approach:** `{approach}` · **{npass}/{len(records)} tests pass** {_mark(npass == len(records))}",
        "",
        "## Tests in this approach",
        "",
        "| Test | Status | Headline | Detail |",
        "|---|---|---|---|",
    ]
    for r in records:
        lines.append(f"| {_test_title(r['test'])} | {_mark(r.get('pass_all'))} | {_headline(r['observed'])} | [report]({r['test']}/REPORT.md) |")
    lines.append("")
    for r in records:
        lines += [f"## {_test_title(r['test'])} — {_badge(r.get('pass_all'))}", ""]
        if r.get("interpretation"):
            lines += [r["interpretation"], ""]
        lines += _metrics_tables(r["observed"], r.get("expected", {}), r.get("tolerance", {}), r.get("within_tol", {}))
        if r.get("caveats"):
            lines += [f"_Caveats: {r['caveats']}_", ""]
    return "\n".join(lines) + "\n"


def render_index(idx: dict) -> str:
    approaches = idx["approaches"]
    all_tests = [t for a in approaches.values() for t in a["tests"].values()]
    npass = sum(1 for t in all_tests if t.get("pass_all") is True)
    lines = [
        "# Discern — Results Index",
        "",
        f"**{npass}/{len(all_tests)} tests pass** ✅ · _generated {idx.get('generated', '')}_",
        "",
        "Navigate: trained models live in `models/`, other approaches in `approaches/`. "
        "Pick one → its summary page → individual test reports.",
        "",
    ]
    by_type: dict[str, list] = {}
    for key, a in approaches.items():
        by_type.setdefault(a["type"], []).append((key, a))
    for typ in TYPE_ORDER:
        if typ not in by_type:
            continue
        lines += [f"## {TYPE_LABELS.get(typ, typ.title())}", "", "| Approach | Tests | Pass | Summary |", "|---|---|---|---|"]
        for key, a in sorted(by_type[typ]):
            tests = list(a["tests"].values())
            np_ = sum(1 for t in tests if t.get("pass_all") is True)
            lines.append(f"| {a['name']} | {len(tests)} | {np_}/{len(tests)} | [{a['summary_file']}]({a['container']}/{key}/{a['summary_file']}) |")
        lines.append("")
    return "\n".join(lines) + "\n"


def _load_index(results_root: Path) -> dict:
    p = results_root / "index.json"
    if p.exists():
        idx = json.loads(p.read_text(encoding="utf-8"))
        if "approaches" in idx:
            return idx
    return {"approaches": {}, "generated": ""}


def _rebuild_approach(results_root: Path, approach: str):
    adir = _approach_dir(results_root, approach)
    records = []
    for jf in sorted(adir.glob("*/results.json")):
        records.append(json.loads(jf.read_text(encoding="utf-8")))
    if not records:
        return
    fn = _summary_filename(approach)
    (adir / fn).write_text(render_approach_summary(approach, records), encoding="utf-8")


def _rebuild_index(results_root: Path):
    idx = {"approaches": {}, "generated": datetime.now(timezone.utc).isoformat()}
    for croot in _container_roots(results_root):
        if not croot.exists():
            continue
        for adir in sorted(p for p in croot.iterdir() if p.is_dir()):
            approach = adir.name
            container = croot.name
            name, typ = _meta(approach)
            tests = {}
            for jf in sorted(adir.glob("*/results.json")):
                r = json.loads(jf.read_text(encoding="utf-8"))
                tests[r["test"]] = {"path": f"{container}/{approach}/{r['test']}/results.json",
                                    "pass_all": r.get("pass_all"), "headline": _headline(r["observed"]),
                                    "observed": r["observed"]}
            if tests:
                idx["approaches"][approach] = {"name": name, "type": typ, "container": container,
                                               "summary_file": _summary_filename(approach), "tests": tests}
    (results_root / "index.json").write_text(json.dumps(idx, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    (results_root / "INDEX.md").write_text(render_index(idx), encoding="utf-8")


def write_result(model, test, observed, *, expected=None, tolerance=None, tables=None,
                 method="", interpretation="", caveats="", reproduce="", supports="",
                 extra=None, results_root: Path = Path("results")) -> Path:
    approach = _approach(model)
    run_dir = _approach_dir(results_root, approach) / test
    run_dir.mkdir(parents=True, exist_ok=True)
    expected, tolerance = expected or {}, tolerance or {}
    diff = {k: observed[k] - expected[k] for k in expected if k in observed and _is_num(observed[k])}
    within = {k: abs(observed[k] - expected[k]) <= tolerance[k] for k in expected if k in observed and k in tolerance}
    pass_all = all(within.values()) if within else None

    csv_files = []
    if tables:
        for nm, df in tables.items():
            fn = f"{nm}.csv"
            df.to_csv(run_dir / fn, index=False)
            csv_files.append(fn)
    if not csv_files:
        pd.DataFrame([observed]).to_csv(run_dir / "results.csv", index=False)
        csv_files = ["results.csv"]

    record = {
        "model": approach, "test": test,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "git_commit": _git_commit(), "data_manifest_sha256": _manifest_hash(results_root),
        "versions": _versions(),
        "observed": observed, "expected": expected, "tolerance": tolerance,
        "diff": diff, "within_tol": within, "pass_all": pass_all,
        "method": method, "interpretation": interpretation, "caveats": caveats,
        "reproduce": reproduce, "supports": supports, "csv_files": csv_files, "extra": extra or {},
    }
    (run_dir / "results.json").write_text(json.dumps(record, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    (run_dir / "REPORT.md").write_text(render_report(record), encoding="utf-8")
    _rebuild_approach(results_root, approach)
    _rebuild_index(results_root)
    return run_dir / "results.json"


def regenerate_all(results_root: Path = Path("results")) -> int:
    n = 0
    for croot in _container_roots(results_root):
        if not croot.exists():
            continue
        for jf in sorted(croot.glob("*/*/results.json")):
            rec = json.loads(jf.read_text(encoding="utf-8"))
            (jf.parent / "REPORT.md").write_text(render_report(rec), encoding="utf-8")
            n += 1
        for adir in sorted(p for p in croot.iterdir() if p.is_dir()):
            _rebuild_approach(results_root, adir.name)
    _rebuild_index(results_root)
    return n


def migrate(results_root: Path = Path("results")) -> list[str]:
    """Place every approach folder in its correct container and re-key results.json.

    type=model -> results/models/<name>/ ; everything else -> results/approaches/<name>/.
    Handles legacy flat names (ALIAS), legacy folders under approaches/, and models that
    were previously parked under approaches/. Idempotent; never overwrites an existing
    canonical test folder.
    """
    moved = []
    for croot in _container_roots(results_root):
        croot.mkdir(parents=True, exist_ok=True)
    sources = []
    for base in (results_root, _models_root(results_root), _approaches_root(results_root)):
        if base.exists():
            for p in sorted(base.iterdir()):
                if p.is_dir() and p.name not in (MODELS_DIR, APPROACHES_DIR):
                    sources.append(p)
    for src in sources:
        canonical = ALIAS.get(src.name, src.name)
        target = _approach_dir(results_root, canonical)
        if src.resolve() == target.resolve():
            continue
        target.mkdir(parents=True, exist_ok=True)
        for test_dir in sorted(p for p in src.iterdir() if p.is_dir()):
            dst = target / test_dir.name
            if dst.exists():
                continue
            shutil.move(str(test_dir), str(dst))
            jf = dst / "results.json"
            if jf.exists():
                rec = json.loads(jf.read_text(encoding="utf-8"))
                rec["model"] = canonical
                jf.write_text(json.dumps(rec, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
            moved.append(f"{src.relative_to(results_root)}/{test_dir.name} -> {target.relative_to(results_root)}/{test_dir.name}")
        shutil.rmtree(src)
    return moved
