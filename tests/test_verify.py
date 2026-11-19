import json
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "eval" / "verify.py"


def write_fixture(
    root: Path, *, observed: int = 11, tolerance: dict | None = None
) -> Path:
    results = root / "results"
    path = Path("approaches/example/accuracy/results.json")
    report = results / path
    report.parent.mkdir(parents=True)
    report.write_text(
        json.dumps(
            {
                "model": "example",
                "test": "accuracy",
                "observed": {"count": observed},
                "expected": {"count": 10},
                "tolerance": {"count": 1} if tolerance is None else tolerance,
                "within_tol": {"count": True},
                "pass_all": True,
            }
        ),
        encoding="utf-8",
    )
    index = {
        "approaches": {
            "example": {
                "tests": {"accuracy": {"path": path.as_posix(), "pass_all": True}}
            }
        }
    }
    (results / "index.json").write_text(json.dumps(index), encoding="utf-8")
    return report


def run_verifier(root: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT)],
        cwd=root,
        capture_output=True,
        text=True,
        check=False,
    )


def test_accepts_value_at_tolerance_boundary(tmp_path: Path) -> None:
    write_fixture(tmp_path)

    run = run_verifier(tmp_path)

    assert run.returncode == 0, run.stderr
    assert "example/accuracy" in run.stdout
    assert "1/1 acceptance tests passed" in run.stdout
    assert not (tmp_path / "results/ACCEPTANCE_FAILURES.md").exists()


def test_recomputes_tolerance_and_reports_exact_diff(tmp_path: Path) -> None:
    write_fixture(tmp_path, observed=13)

    run = run_verifier(tmp_path)

    assert run.returncode == 1, run.stderr
    assert "FAIL" in run.stdout
    report = (tmp_path / "results/ACCEPTANCE_FAILURES.md").read_text(encoding="utf-8")
    assert "| example/accuracy | count | 13 | 10 | 3 | 1 |" in report
    assert "Research use only" in report


def test_marks_descriptive_records_without_counting_them(tmp_path: Path) -> None:
    write_fixture(tmp_path)
    results = tmp_path / "results"
    path = "approaches/example/audit/results.json"
    report = results / path
    report.parent.mkdir()
    report.write_text(
        json.dumps(
            {
                "model": "example",
                "test": "audit",
                "observed": {"summary": "retrospective"},
                "expected": {},
                "tolerance": {},
                "within_tol": {},
                "pass_all": None,
            }
        ),
        encoding="utf-8",
    )
    index = json.loads((results / "index.json").read_text(encoding="utf-8"))
    index["approaches"]["example"]["tests"]["audit"] = {"path": path, "pass_all": None}
    (results / "index.json").write_text(json.dumps(index), encoding="utf-8")

    run = run_verifier(tmp_path)

    assert run.returncode == 0, run.stderr
    assert "DESCRIPTIVE" in run.stdout
    assert "1/1 acceptance tests passed" in run.stdout


@pytest.mark.parametrize(
    "broken", ["missing_report", "missing_tolerance", "stale_index"]
)
def test_rejects_incomplete_or_inconsistent_acceptance_data(
    tmp_path: Path, broken: str
) -> None:
    report = write_fixture(
        tmp_path, tolerance={} if broken == "missing_tolerance" else None
    )
    if broken == "missing_report":
        report.unlink()
    if broken == "stale_index":
        index_path = tmp_path / "results/index.json"
        index = json.loads(index_path.read_text(encoding="utf-8"))
        index["approaches"]["example"]["tests"]["accuracy"]["pass_all"] = None
        index_path.write_text(json.dumps(index), encoding="utf-8")

    run = run_verifier(tmp_path)

    assert run.returncode == 1, run.stderr
    assert "example/accuracy" in run.stdout
    assert "FAIL" in run.stdout


def test_fails_if_no_acceptance_tests_are_registered(tmp_path: Path) -> None:
    results = tmp_path / "results"
    results.mkdir()
    (results / "index.json").write_text('{"approaches": {}}', encoding="utf-8")

    run = run_verifier(tmp_path)

    assert run.returncode == 1, run.stderr
    assert "No acceptance tests" in run.stdout
