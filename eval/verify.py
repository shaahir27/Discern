"""Verify saved acceptance results without rerunning experiments."""

import json
from math import isfinite
from pathlib import Path
from typing import Any

BANNER = "Research use only"


def read_object(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise TypeError(f"{path.name} must contain a JSON object")
    return value


def check_result(
    root: Path, model: str, test: str, entry: dict[str, Any], misses: list[str]
) -> tuple[str, str]:
    path = (root / entry["path"]).resolve()
    if not path.is_relative_to(root):
        raise ValueError("result path leads outside the results directory")
    record = read_object(path)
    if (record["model"], record["test"]) != (model, test):
        raise ValueError("report identity differs from its index entry")

    observed, expected, tolerance, stored_checks = (
        record[key] for key in ("observed", "expected", "tolerance", "within_tol")
    )
    if not all(
        isinstance(field, dict)
        for field in (observed, expected, tolerance, stored_checks)
    ):
        raise TypeError("metric fields must be JSON objects")
    if expected.keys() != tolerance.keys():
        raise ValueError("expected metrics and tolerances must have matching keys")

    checks = {}
    for metric, target in expected.items():
        value, limit = observed[metric], tolerance[metric]
        if not all(
            type(number) in (int, float) and isfinite(number)
            for number in (value, target, limit)
        ):
            raise ValueError(f"{metric} requires finite numeric values")
        if limit < 0:
            raise ValueError(f"{metric} requires a nonnegative tolerance")
        diff = value - target
        checks[metric] = abs(diff) <= limit
        if not checks[metric]:
            misses.append(
                f"| {model}/{test} | {metric} | {value!r} | {target!r} | {diff!r} | {limit!r} |"
            )

    passed = all(checks.values()) if checks else None
    if (
        record["pass_all"] is not passed
        or entry["pass_all"] is not passed
        or stored_checks.keys() != checks.keys()
        or any(stored_checks[metric] is not ok for metric, ok in checks.items())
    ):
        return "FAIL", "saved pass flags differ from recalculated checks"
    if passed is None:
        return "DESCRIPTIVE", "no acceptance criteria"
    if not passed:
        return "FAIL", f"{sum(not ok for ok in checks.values())} outside tolerance"
    return "PASS", f"{len(checks)} within tolerance"


def main() -> int:
    root = Path("results").resolve()
    print(BANNER)
    try:
        index = read_object(root / "index.json")
        entries = [
            (model, test, entry)
            for model, approach in index["approaches"].items()
            for test, entry in approach["tests"].items()
        ]
    except (OSError, AttributeError, KeyError, TypeError, ValueError) as exc:
        print(f"FAIL | results/index.json | {exc}")
        return 1

    print("Status | Test | Checks")
    misses: list[str] = []
    passed = failed = 0
    for model, test, entry in entries:
        try:
            status, detail = check_result(root, model, test, entry, misses)
        except (OSError, KeyError, TypeError, ValueError) as exc:
            status, detail = "FAIL", str(exc)
        passed += status == "PASS"
        failed += status == "FAIL"
        print(f"{status} | {model}/{test} | {detail}")

    if misses:
        report = root / "ACCEPTANCE_FAILURES.md"
        report.write_text(
            "\n".join(
                [
                    "# Acceptance failures",
                    "",
                    f"> {BANNER}",
                    "",
                    "| Test | Metric | Observed | Expected | Difference | Tolerance |",
                    "|---|---|---|---|---|---|",
                    *misses,
                    "",
                ]
            ),
            encoding="utf-8",
        )
        print("See results/ACCEPTANCE_FAILURES.md for differences.")
    if passed + failed == 0:
        print("No acceptance tests registered.")
        return 1
    print(f"{passed}/{passed + failed} acceptance tests passed.")
    return int(failed > 0)


if __name__ == "__main__":
    raise SystemExit(main())
