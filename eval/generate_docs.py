"""Regenerate public evaluation documents from the report registry."""

import sys
from pathlib import Path

sys.path.insert(0, "src")

from discern.documentation import model_card
from discern.research import BANNER, Reports


def main() -> None:
    root = Path(".")
    reports = Reports(root)
    card = ["# Discern model card", BANNER]
    for heading, paragraphs in model_card(reports):
        card.extend([f"## {heading}", *paragraphs])
    (root / "docs/MODEL_CARD.md").write_text("\n\n".join(card) + "\n", encoding="utf-8", newline="\n")

    results = ["# Evaluation results", "", BANNER, "", "Metrics are reproduced from the indexed result JSON.", ""]
    for key, path in reports.paths.items():
        record = reports.get(key)
        status = {True: "Passed", False: "Failed", None: "Descriptive"}[record["pass_all"]]
        results.extend([f"## {key}", "", f"{status} · [Source JSON](../{path.as_posix()})", ""])
        if record["test"] in ("e08_novel_variant", "e10_prevalence_shift", "e14_certificate_audit",
                              "e15_label_reveal_sim"):
            results.extend([record["caveats"], ""])
        results.extend(["| Metric | Observed | Expected | Tolerance |", "|---|---|---|---|"])
        for metric, value in record["observed"].items():
            if isinstance(value, (str, int, float)) and "|" not in str(value):
                results.append(f"| `{metric}` | {value} | {record['expected'].get(metric, '—')} | "
                               f"{record['tolerance'].get(metric, '—')} |")
        tables = reports.tables(key, record)
        if tables:
            links = " · ".join(f"[{name}](../{path.with_name(name).as_posix()})" for name in tables)
            results.extend(["", links])
        results.append("")
    (root / "docs/RESULTS.md").write_text("\n".join(results), encoding="utf-8")


if __name__ == "__main__":
    main()
