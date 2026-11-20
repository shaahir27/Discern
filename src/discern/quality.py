"""Live CSV quality check for uploaded variant profiles and sample rosters.

Computes completeness scores from required fields and identifies rows
with missing values, duplicates, or invalid formatting.
"""

from __future__ import annotations

import csv
import io
import re
from decimal import Decimal, InvalidOperation
from typing import Any

from .grammar import valid_mutation

REQUIRED_PROFILE_COLUMNS = ["sample_id", "gene", "mutation", "frs"]
ALLOWED_GENES = {"Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5"}
SAMPLE_ID_REGEX = re.compile(r"^[A-Za-z0-9_\-\.]{1,64}$")


def check_csv_quality(payload: bytes | str) -> dict[str, Any]:
    """Analyze a CSV payload and return a comprehensive quality report."""
    if isinstance(payload, bytes):
        try:
            text = payload.decode("utf-8-sig")
        except UnicodeError:
            return {
                "completeness_score": 0.0,
                "total_rows": 0,
                "valid_rows": 0,
                "total_cells": 0,
                "filled_cells": 0,
                "missing_cells": 0,
                "missing_count": 0,
                "duplicate_count": 0,
                "formatting_count": 1,
                "total_issues": 1,
                "issues": [{
                    "line": 1,
                    "type": "format",
                    "field": "encoding",
                    "message": "File is not valid UTF-8. Save the file with UTF-8 encoding.",
                }],
                "is_valid": False,
                "summary": "Invalid file encoding.",
                "format": "unknown",
                "headers": [],
            }
    else:
        text = payload

    stripped_text = text.strip()
    if not stripped_text:
        return {
            "completeness_score": 0.0,
            "total_rows": 0,
            "valid_rows": 0,
            "total_cells": 0,
            "filled_cells": 0,
            "missing_cells": 0,
            "missing_count": 0,
            "duplicate_count": 0,
            "formatting_count": 1,
            "total_issues": 1,
            "issues": [{
                "line": 1,
                "type": "format",
                "field": "file",
                "message": "The uploaded CSV file is empty.",
            }],
            "is_valid": False,
            "summary": "Empty CSV file.",
            "format": "unknown",
            "headers": [],
        }

    try:
        reader = csv.reader(io.StringIO(text, newline=""), strict=True)
        raw_header = next(reader, None)
    except csv.Error as err:
        return {
            "completeness_score": 0.0,
            "total_rows": 0,
            "valid_rows": 0,
            "total_cells": 0,
            "filled_cells": 0,
            "missing_cells": 0,
            "missing_count": 0,
            "duplicate_count": 0,
            "formatting_count": 1,
            "total_issues": 1,
            "issues": [{
                "line": 1,
                "type": "format",
                "field": "csv",
                "message": f"Malformed CSV structure: {err}",
            }],
            "is_valid": False,
            "summary": "Malformed CSV file.",
            "format": "unknown",
            "headers": [],
        }

    if not raw_header:
        return {
            "completeness_score": 0.0,
            "total_rows": 0,
            "valid_rows": 0,
            "total_cells": 0,
            "filled_cells": 0,
            "missing_cells": 0,
            "missing_count": 0,
            "duplicate_count": 0,
            "formatting_count": 1,
            "total_issues": 1,
            "issues": [{
                "line": 1,
                "type": "format",
                "field": "header",
                "message": "No header row found.",
            }],
            "is_valid": False,
            "summary": "Missing header row.",
            "format": "unknown",
            "headers": [],
        }

    headers = [col.strip() for col in raw_header]
    has_banner = headers[-1:] == ["banner"]
    if has_banner:
        headers = headers[:-1]

    # Determine required columns
    is_roster = headers == ["sample_id"]
    is_profile = all(col in headers for col in REQUIRED_PROFILE_COLUMNS) or set(headers).issubset(set(REQUIRED_PROFILE_COLUMNS))
    required_cols = ["sample_id"] if is_roster else REQUIRED_PROFILE_COLUMNS
    csv_format = "roster" if is_roster else ("profile" if is_profile else "custom")

    issues: list[dict[str, Any]] = []

    # Check for missing required columns in header
    for req_col in required_cols:
        if req_col not in headers:
            issues.append({
                "line": 1,
                "type": "format",
                "field": req_col,
                "message": f"Header line 1 is missing required column '{req_col}'",
            })

    rows_data: list[tuple[int, list[str]]] = []
    line_num = 1
    while True:
        try:
            row_values = next(reader, None)
            if row_values is None:
                break
            line_num = reader.line_num
            # skip completely empty lines
            if not row_values or all(val.strip() == "" for val in row_values):
                continue
            rows_data.append((line_num, [v.strip() for v in row_values]))
        except csv.Error as err:
            line_num += 1
            issues.append({
                "line": line_num,
                "type": "format",
                "field": "csv",
                "message": f"Line {line_num}: malformed CSV quoting or separator: {err}",
            })

    total_rows = len(rows_data)
    if total_rows == 0:
        issues.append({
            "line": 1,
            "type": "format",
            "field": "data",
            "message": "No data rows found in the CSV file.",
        })
        return {
            "completeness_score": 0.0,
            "total_rows": 0,
            "valid_rows": 0,
            "total_cells": 0,
            "filled_cells": 0,
            "missing_cells": 0,
            "missing_count": 0,
            "duplicate_count": 0,
            "formatting_count": len(issues),
            "total_issues": len(issues),
            "issues": issues,
            "is_valid": False,
            "summary": "No data rows found.",
            "format": csv_format,
            "headers": headers,
        }

    # Tracking for completeness and duplicates
    total_required_cells = total_rows * len(required_cols)
    filled_required_cells = 0
    seen_exact_rows: dict[tuple[str, ...], int] = {}
    seen_variants: dict[tuple[str, str, str], int] = {}
    rows_with_issues: set[int] = set()

    for line_idx, vals in rows_data:
        if has_banner and len(vals) == len(headers) + 1:
            vals = vals[:-1]
        # Check column count
        if len(vals) < len(headers):
            issues.append({
                "line": line_idx,
                "type": "format",
                "field": "columns",
                "message": f"Line {line_idx}: has {len(vals)} columns, expected {len(headers)}",
            })
            rows_with_issues.add(line_idx)
        elif len(vals) > len(headers):
            issues.append({
                "line": line_idx,
                "type": "format",
                "field": "columns",
                "message": f"Line {line_idx}: has extra {len(vals) - len(headers)} columns",
            })
            rows_with_issues.add(line_idx)

        row_dict = dict(zip(headers, vals))

        # Check required fields completeness and missing values
        for req_col in required_cols:
            val = row_dict.get(req_col, "")
            if val != "":
                filled_required_cells += 1
            else:
                issues.append({
                    "line": line_idx,
                    "type": "missing",
                    "field": req_col,
                    "message": f"Line {line_idx}: missing required value in '{req_col}'",
                })
                rows_with_issues.add(line_idx)

        # Check duplicates
        exact_key = tuple(row_dict.get(h, "") for h in headers)
        exact_duplicate = exact_key in seen_exact_rows
        if exact_duplicate:
            prior_line = seen_exact_rows[exact_key]
            issues.append({
                "line": line_idx,
                "type": "duplicate",
                "field": "row",
                "message": f"Line {line_idx}: exact duplicate of line {prior_line}",
            })
            rows_with_issues.add(line_idx)
        else:
            seen_exact_rows[exact_key] = line_idx

        # Check variant duplicates (same sample, gene, mutation)
        if csv_format == "profile" and "sample_id" in row_dict and "gene" in row_dict and "mutation" in row_dict:
            var_key = (row_dict["sample_id"], row_dict["gene"], row_dict["mutation"])
            if all(k != "" for k in var_key):
                if var_key in seen_variants and not exact_duplicate:
                    prior_line = seen_variants[var_key]
                    issues.append({
                        "line": line_idx,
                        "type": "duplicate",
                        "field": "variant",
                        "message": f"Line {line_idx}: repeated variant for this sample (earlier at line {prior_line})",
                    })
                    rows_with_issues.add(line_idx)
                elif var_key not in seen_variants:
                    seen_variants[var_key] = line_idx

        # Check formatting of specific fields
        # 1. sample_id
        sid = row_dict.get("sample_id", "")
        if sid and not SAMPLE_ID_REGEX.fullmatch(sid):
            issues.append({
                "line": line_idx,
                "type": "format",
                "field": "sample_id",
                "message": f"Line {line_idx}: use up to 64 letters, digits, dots, hyphens, or underscores",
            })
            rows_with_issues.add(line_idx)

        # 2. gene
        if "gene" in headers:
            gene = row_dict.get("gene", "")
            if gene != "" and gene not in ALLOWED_GENES:
                issues.append({
                    "line": line_idx,
                    "type": "format",
                    "field": "gene",
                    "message": f"Line {line_idx}: expected one of {', '.join(sorted(ALLOWED_GENES))}",
                })
                rows_with_issues.add(line_idx)

        # 3. mutation
        if "mutation" in headers:
            mut = row_dict.get("mutation", "")
            if mut != "":
                prefix, sep, suffix = mut.partition("@")
                if (sep and prefix != row_dict.get("gene")) or not valid_mutation(suffix if sep else prefix):
                    issues.append({
                        "line": line_idx,
                        "type": "format",
                        "field": "mutation",
                        "message": f"Line {line_idx}: invalid mutation syntax",
                    })
                    rows_with_issues.add(line_idx)

        # 4. frs
        if "frs" in headers:
            frs_str = row_dict.get("frs", "")
            if frs_str != "":
                try:
                    dec = Decimal(frs_str)
                    if not dec.is_finite() or dec < 0 or dec > 1:
                        issues.append({
                            "line": line_idx,
                            "type": "format",
                            "field": "frs",
                            "message": f"Line {line_idx}: frs must be between 0.0 and 1.0",
                        })
                        rows_with_issues.add(line_idx)
                except InvalidOperation:
                    issues.append({
                        "line": line_idx,
                        "type": "format",
                        "field": "frs",
                        "message": f"Line {line_idx}: frs must be a valid numeric fraction",
                    })
                    rows_with_issues.add(line_idx)

        # 5. collection_date
        if "collection_date" in headers:
            cdate = row_dict.get("collection_date", "").strip()
            if cdate != "":
                from .longitudinal import parse_timepoint_date
                parsed_val, _ = parse_timepoint_date(cdate)
                if parsed_val is None:
                    issues.append({
                        "line": line_idx,
                        "type": "format",
                        "field": "collection_date",
                        "message": f"Line {line_idx}: enter a YYYY-MM-DD date or a non-negative day offset",
                    })
                    rows_with_issues.add(line_idx)

    completeness_score = round(
        (filled_required_cells / total_required_cells * 100.0) if total_required_cells > 0 else 0.0,
        1,
    )

    missing_count = sum(1 for iss in issues if iss["type"] == "missing")
    duplicate_count = sum(1 for iss in issues if iss["type"] == "duplicate")
    formatting_count = sum(1 for iss in issues if iss["type"] == "format")

    valid_rows = total_rows - len(rows_with_issues)

    if not issues:
        summary = f"All {total_rows} rows passed validation with 100% completeness."
    else:
        summary = f"{len(issues)} issue(s) detected across {len(rows_with_issues)} of {total_rows} row(s)."

    return {
        "completeness_score": completeness_score,
        "total_rows": total_rows,
        "valid_rows": max(0, valid_rows),
        "total_cells": total_required_cells,
        "filled_cells": filled_required_cells,
        "missing_cells": total_required_cells - filled_required_cells,
        "missing_count": missing_count,
        "duplicate_count": duplicate_count,
        "formatting_count": formatting_count,
        "total_issues": len(issues),
        "issues": issues,
        "is_valid": len(issues) == 0,
        "summary": summary,
        "format": csv_format,
        "headers": headers,
    }
