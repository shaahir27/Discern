"""Longitudinal serial-isolate tracking and subclonal trajectory analysis."""
from __future__ import annotations

import datetime
import math
import re
from typing import Any

from .grammar import classify, suffix_of

CORE_GENES = {"Rv0678", "atpE", "pepQ"}


def parse_timepoint_date(val: Any) -> tuple[float | None, str]:
    """Parse collection date into a numeric sort key (days/ordinal) and human display string."""
    if val is None:
        return None, ""
    if isinstance(val, (int, float)):
        if isinstance(val, bool) or not math.isfinite(val) or val < 0 or val != int(val):
            return None, ""
        return float(val), f"Day {int(val)}"
    
    s = str(val).strip()
    if not s:
        return None, ""

    # Integer or Day X
    m_day = re.match(r"^(?:day\s*)?(\d+)$", s, re.IGNORECASE)
    if m_day and len(m_day.group(1)) <= 9:
        d = int(m_day.group(1))
        return float(d), f"Day {d}"

    # ISO date: YYYY-MM-DD
    m_iso = re.match(r"^(\d{4})-(\d{2})-(\d{2})$", s)
    if m_iso:
        try:
            dt = datetime.date(int(m_iso.group(1)), int(m_iso.group(2)), int(m_iso.group(3)))
            return float(dt.toordinal()), s
        except ValueError:
            return None, s

    return None, s


def analyze_patient_timelines(sample_records: list[dict]) -> dict[str, dict]:
    """Analyze multi-visit isolates per patient and identify subclonal emergence/competition."""
    by_patient: dict[str, list[dict]] = {}
    for r in sample_records:
        pid = str(r.get("patient_id") or "").strip()
        if not pid:
            continue
        by_patient.setdefault(pid, []).append(r)

    timelines: dict[str, dict] = {}

    for pid, records in by_patient.items():
        if len(records) < 2:
            continue

        # Parse date and sort visits
        parsed_records = []
        for r in records:
            d_key, d_disp = parse_timepoint_date(r.get("collection_date"))
            if d_key is None:
                continue
            variants = {}
            for variant in r.get("variants", []):
                mutation = suffix_of(variant.get("mutation", ""))
                flags = classify(mutation)
                if mutation and not (flags["IS_SYN"] or flags["IS_NULL"]):
                    key = (variant.get("gene", ""), mutation)
                    frs = float(variant.get("frs", 0.0))
                    if key not in variants or frs > variants[key]["frs"]:
                        variants[key] = {**variant, "mutation": mutation, "frs": frs}
            parsed_records.append({
                "sort_key": d_key,
                "date_display": d_disp,
                "sample_id": r.get("sample_id", ""),
                "call": r.get("call", "U"),
                "call_label": r.get("call_label", ""),
                "probability": r.get("probability"),
                "variants": list(variants.values()),
            })

        if len(parsed_records) < 2 or len({r["date_display"].startswith("Day ") for r in parsed_records}) > 1:
            continue
        parsed_records.sort(key=lambda x: (x["sort_key"], x["sample_id"]))

        # Identify all distinct variants seen across timepoints
        all_variants: dict[str, str] = {}  # mut_key -> gene
        for visit in parsed_records:
            for v in visit["variants"]:
                gene = v.get("gene", "")
                mut = v.get("mutation", "")
                if gene and mut:
                    mut_key = f"{gene}@{mut}"
                    all_variants[mut_key] = gene

        var_list = sorted(all_variants.keys())

        # Construct trajectory matrix: variant -> list of frs across visits
        trajectories: dict[str, list[float]] = {m: [] for m in var_list}
        for visit in parsed_records:
            v_map = {f"{v.get('gene')}@{v.get('mutation')}": float(v.get("frs", 0.0)) for v in visit["variants"]}
            for m in var_list:
                trajectories[m].append(v_map.get(m, 0.0))

        # Compute wild-type remainder per visit
        core_variants = [m for m in var_list if all_variants[m] in CORE_GENES]
        core_sums = [sum(trajectories[m][j] for m in core_variants) for j in range(len(parsed_records))]
        wt_series = [max(0.0, round(1.0 - value, 4)) for value in core_sums]

        # Evaluate longitudinal alerts
        alerts = []
        # 1. DE_NOVO_VARIANT_EMERGENCE
        for m in var_list:
            if all_variants[m] in CORE_GENES:
                series = trajectories[m]
                if series[0] == 0.0:
                    for j in range(1, len(series)):
                        if series[j] >= 0.10:
                            visit_date = parsed_records[j]["date_display"]
                            alerts.append({
                                "code": "DE_NOVO_VARIANT_EMERGENCE",
                                "variant": m,
                                "date": visit_date,
                                "visit_idx": j,
                                "message": f"De novo emergence of {m} at {visit_date} (FRS: {series[j]:.2f}).",
                            })
                            break

        # 2. RISING_SUBCLONE_TRAJECTORY
        for m in var_list:
            if all_variants[m] in CORE_GENES:
                series = trajectories[m]
                for j in range(1, len(series)):
                    delta = series[j] - series[j - 1]
                    crossed_threshold = (series[j - 1] < 0.75 and series[j] >= 0.75 and series[j] >= 0.10)
                    if (delta >= 0.15 and series[j] >= 0.10) or crossed_threshold:
                        visit_date = parsed_records[j]["date_display"]
                        alerts.append({
                            "code": "RISING_SUBCLONE_TRAJECTORY",
                            "variant": m,
                            "date": visit_date,
                            "visit_idx": j,
                            "message": f"Rising subclone trajectory for {m} at {visit_date} ({series[j-1]:.2f} → {series[j]:.2f}).",
                        })

        # 3. POLYCLONAL_RV0678_COMPETITION
        rv_variants = [m for m in var_list if all_variants[m] == "Rv0678" and max(trajectories[m]) >= 0.10]
        if len(rv_variants) >= 2:
            alerts.append({
                "code": "POLYCLONAL_RV0678_COMPETITION",
                "variant": ", ".join(rv_variants),
                "date": "Cohort",
                "message": f"Polyclonal Rv0678 competition detected across {len(rv_variants)} alleles ({', '.join(rv_variants)}).",
            })

        # 4. CLONAL_REPLACEMENT
        if len(rv_variants) >= 2:
            for j in range(1, len(parsed_records)):
                dropped = [m for m in rv_variants if trajectories[m][j - 1] - trajectories[m][j] >= 0.15]
                risen = [m for m in rv_variants if trajectories[m][j] - trajectories[m][j - 1] >= 0.15]
                if dropped and risen:
                    visit_date = parsed_records[j]["date_display"]
                    alerts.append({
                        "code": "CLONAL_REPLACEMENT",
                        "variant": f"{', '.join(dropped)} replaced by {', '.join(risen)}",
                        "date": visit_date,
                        "visit_idx": j,
                        "message": f"Clonal replacement at {visit_date}: {', '.join(dropped)} dropped while {', '.join(risen)} expanded.",
                    })

        timelines[pid] = {
            "patient_id": pid,
            "n_visits": len(parsed_records),
            "timepoints": parsed_records,
            "variants": var_list,
            "core_variants": core_variants,
            "trajectories": trajectories,
            "wt_series": wt_series,
            "has_overlapping_fractions": any(value > 1.000000001 for value in core_sums),
            "alerts": alerts,
            "has_alerts": len(alerts) > 0,
        }

    return timelines
