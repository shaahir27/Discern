"""Read immutable research reports and a fixed public replay cohort."""

import csv
import io
import json
import os
from collections import Counter
from decimal import Decimal, InvalidOperation
from itertools import islice
from pathlib import Path
from typing import Any

import pandas as pd

from .evidence import EvidenceStore

BANNER = "Research use only. For evaluation and education. Clinical decisions require laboratory confirmation."
MODE = "held-out in cross-validation"
CALL_LABELS = {"R": "Resistant", "S": "Susceptible", "U": "Send for lab test"}
REPORTS = {
    "prep": "data_prep/e01_prep",
    "catalogue": "catalogue/e03_s1_catalogue_cv",
    "model": "lr_model/e04_s1_patient_cv",
    "mechanism": "lr_mechanism/e16_lof_ablation",
    "km": "external_km/e_km_validation",
    "conformal": "conformal/e14_certificate_audit",
    "priority": "priority_list/e12_variant_priority",
    "ship": "ship_decision/e19_ladder_summary",
    "intervals": "research_archive/g8_evaluation_context",
}
MAX_CSV_BYTES = 16_384
PRIVATE_HEADERS = {"uniqueid", "sample_id", "subjid", "accession", "run_accession"}
PROFILE_COLUMNS = ["sample_id", "gene", "mutation", "frs"]
OPTIONAL_PROFILE_COLUMNS = {"patient_id", "collection_date", "assay"}


class UploadError(ValueError):
    """A user-correctable upload failure that carries no submitted values."""

    def __init__(
        self, code: str, message: str, *, line: int | None = None,
        column: str | None = None, status: int = 422,
    ):
        super().__init__(message)
        self.code = code
        self.line = line
        self.column = column
        self.status = status

    def __str__(self) -> str:
        if self.line is None:
            return super().__str__()
        location = f"CSV line {self.line}"
        if self.column is not None:
            location += f", column {self.column}"
        return f"{location}. {super().__str__()}"


def workspace() -> Path:
    return Path(os.environ.get("DISCERN_ROOT", ".")).resolve()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def csv_text(columns: list[str], rows: list[dict]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, [*columns, "banner"], extrasaction="ignore", lineterminator="\n")
    writer.writeheader()
    writer.writerows({**row, "banner": BANNER} for row in rows)
    return output.getvalue()


def read_batch(payload: bytes) -> tuple[list[str], list[tuple[int, dict[str, str]]]]:
    if len(payload) > MAX_CSV_BYTES:
        raise UploadError("file_too_large", "Choose a smaller CSV of recorded sample IDs.", status=413)
    try:
        text = payload.decode("utf-8-sig")
    except UnicodeError:
        raise UploadError("invalid_encoding", "Save the file as UTF-8 CSV and try again.") from None
    reader = csv.reader(io.StringIO(text, newline=""), strict=True)
    rows = []
    line = 1
    try:
        fields = next(reader, [])
        if not fields and not text.strip():
            raise UploadError("empty_csv", "Add at least one recorded sample ID.")
        headers = fields[:-1] if fields[-1:] == ["banner"] else fields
        is_valid_header = (
            headers == ["sample_id"]
            or headers == PROFILE_COLUMNS
            or (len(headers) >= 4 and headers[:4] == PROFILE_COLUMNS and set(headers[4:]).issubset(OPTIONAL_PROFILE_COLUMNS))
        )
        if not is_valid_header:
            raise UploadError(
                "invalid_header", "Use the columns in the example CSV or a list with a sample_id column.", line=1,
            )
        req_headers = headers[:4] if len(headers) >= 4 and headers[:4] == PROFILE_COLUMNS else headers
        while (values := next(reader, None)) is not None:
            if values:
                if len(values) > len(fields):
                    raise UploadError("extra_cells", "Remove the extra cells from this row.", line=line + 1)
                row = dict(zip(headers, (value.strip() for value in values)))
                for column in req_headers:
                    if not row.get(column):
                        raise UploadError("missing_cell", "Fill in this required cell.", line=line + 1, column=column)
                rows.append((line + 1, row))
            line = reader.line_num
    except csv.Error:
        raise UploadError("malformed_csv", "Check the CSV quoting and separators.", line=reader.line_num or 1) from None
    if not rows:
        raise UploadError("empty_csv", "Add at least one recorded sample ID.")
    return headers, rows


class Reports:
    def __init__(self, root: Path):
        self.root = root / "results"
        index = read_json(self.root / "index.json")
        self.approaches = index["approaches"]
        self.paths = {
            f"{approach}/{test}": self.root / record["path"]
            for approach, section in self.approaches.items()
            for test, record in section["tests"].items()
        }

    def get(self, key: str) -> dict:
        return read_json(self.paths[key])

    def observed(self, name: str) -> dict:
        return self.get(REPORTS[name])["observed"]

    def tables(self, key: str, report: dict) -> list[str]:
        names = []
        for name in report["csv_files"]:
            if name != Path(name).name or name == "accessions.csv":
                continue
            with self.paths[key].with_name(name).open(encoding="utf-8-sig", newline="") as stream:
                header = next(csv.reader(stream), [])
            if not (PRIVATE_HEADERS & {column.casefold() for column in header}):
                names.append(name)
        return names

    def table(self, key: str, name: str, page: int = 1, size: int = 25) -> tuple[list[str], list[list[str]], bool]:
        with self.paths[key].with_name(name).open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.reader(stream)
            columns = next(reader, [])
            start = (page - 1) * size
            rows = list(islice(reader, start, start + size + 1))
        return columns, rows[:size], len(rows) > size


class Archive:
    def __init__(self, root: Path):
        self.reports = Reports(root)
        model = root / "models/lr_model"
        self.evidence = EvidenceStore(model / "evidence_store.csv")
        self.versions = read_json(model / "metadata.json")["package_versions"]
        demo = pd.read_csv(root / "demo/held_out_samples.csv", dtype={"frs": str}, keep_default_na=False)
        oof = pd.read_parquet(self.reports.root / "oof_predictions.parquet").set_index("UNIQUEID")
        prep = self.reports.observed("prep")
        if (
            not oof.index.is_unique or oof.index.hasnans or oof.index.empty
            or not oof.y.isin([0, 1]).all() or not oof.p.between(0, 1).all()
            or not oof.cat_call.isin(CALL_LABELS).all()
            or not pd.api.types.is_integer_dtype(oof.fold) or oof.fold.isna().any()
            or oof.SRC.isna().any() or oof.fold.nunique() < 2
            or len(oof) != prep["n_samples"] or oof.y.sum() != prep["n_resistant"]
        ):
            raise ValueError("Replay archive differs from the evaluated cohort.")
        self.score_contexts = {}
        for fold, held_out in oof.groupby("fold").y:
            reference_n = len(oof) - len(held_out)
            reference_r = int(oof.y.sum() - held_out.sum())
            self.score_contexts[fold] = {
                "source_report": REPORTS["model"],
                "reference_prevalence": reference_r / reference_n,
                "reference_n": reference_n,
                "reference_n_resistant": reference_r,
                "uncertainty": None,
                "uncertainty_status": "not_archived",
            }
        self.groups = demo.groupby("sample_id", sort=False)
        self.roster = {}
        for sid, group in self.groups:
            row = oof.loc[sid]
            if (
                not group.cat_call.eq(row.cat_call).all() or not group.y.eq(row.y).all()
                or not group.p.sub(row.p).abs().le(1e-12).all()
            ):
                raise ValueError("Replay fixture differs from the archived predictions.")
            self.roster[sid] = {
                "sample_id": sid, "call": str(row.cat_call), "call_label": CALL_LABELS[row.cat_call],
                "probability": float(row.p), "fold": int(row.fold),
                "source": str(row.SRC), "phenotype": "Resistant" if row.y else "Susceptible",
            }
        self.example_id = next(
            sid for sid, row in self.roster.items()
            if row["call"] == "R" and row["phenotype"] == "Resistant"
        )

    def profile(self, sid: str) -> list[dict]:
        return self.groups.get_group(sid)[PROFILE_COLUMNS].to_dict("records")

    def imported_ids(self, payload: bytes) -> list[str]:
        headers, rows = read_batch(payload)
        message = "We could not match the CSV to a recorded sample. Try the example CSV or check the sample IDs."
        first_lines = {}
        for line, row in rows:
            sid = row["sample_id"]
            if sid not in self.roster:
                raise UploadError("unknown_sample", message, line=line, column="sample_id")
            first_lines.setdefault(sid, line)
        if headers == PROFILE_COLUMNS or (len(headers) >= 4 and headers[:4] == PROFILE_COLUMNS):
            remaining = {
                sid: Counter((row["gene"], row["mutation"], Decimal(row["frs"])) for row in self.profile(sid))
                for sid in first_lines
            }
            for line, row in rows:
                try:
                    fraction = Decimal(row["frs"])
                except InvalidOperation:
                    fraction = Decimal("NaN")
                if not fraction.is_finite() or not 0 <= fraction <= 1:
                    raise UploadError(
                        "invalid_fraction", "Enter a number between zero and one.", line=line, column="frs",
                    )
                expected = remaining[row["sample_id"]]
                key = row["gene"], row["mutation"], fraction
                if expected[key] <= 0:
                    raise UploadError("profile_mismatch", message, line=line)
                expected[key] -= 1
            for sid, counts in remaining.items():
                if counts.total():
                    raise UploadError("incomplete_profile", message, line=first_lines[sid], column="sample_id")
        return list(first_lines)

    def card(self, variant: str) -> dict:
        return self.evidence.card(variant) or {
            "variant": variant, "n_carriers": 0, "n_resistant": 0, "pct_resistant": None,
            "ci95_low": None, "ci95_high": None, "by_source": {},
            "note": "No carrier evidence in the archived training cohort.",
        }

    def replay(self, sid: str) -> dict:
        row = self.roster[sid]
        evidence = []
        group = self.groups.get_group(sid)
        for item in group.drop_duplicates(["gene", "mutation"]).itertuples():
            evidence.append({**self.card(f"{item.gene}@{item.mutation}"), "frs": float(item.frs)})
        return {
            **row, "evidence": evidence, "mode": MODE, "banner": BANNER,
            "score_context": self.score_contexts[row["fold"]].copy(),
            "call_explanation": {
                "status": "not_archived",
                "source_report": REPORTS["catalogue"],
                "contributing_variants": None,
            },
        }

    def queue(self, ids: list[str]) -> list[dict[str, Any]]:
        rows = [self.roster[sid] for sid in ids]
        rows.sort(key=lambda row: (row["call"] != "U", -row["probability"], row["sample_id"]))
        return [{**row, "rank": rank if row["call"] == "U" else None} for rank, row in enumerate(rows, 1)]
