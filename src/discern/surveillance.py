"""Lightweight batch recurrence and contamination surveillance auditing."""
from __future__ import annotations

import re
from typing import Any

from .evidence import EvidenceStore
from .grammar import classify

CORE_AND_MODIFIERS = {"Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5"}


def _is_homopolymer_hotspot(gene: str, mut_suffix: str) -> bool:
    if gene != "Rv0678":
        return False
    m = re.match(r"^(\d+)_(?:ins|del)_", mut_suffix)
    if not m:
        m = re.match(r"^(\d+)(?:dup|del)", mut_suffix)
    if m:
        pos = int(m.group(1))
        if (192 <= pos <= 198) or (138 <= pos <= 144):
            return True
    return False


def batch_recurrence_audit(
    rows_by_sample: dict[str, list[dict]],
    evidence_store: EvidenceStore,
    *,
    patient_map: dict[str, str] | None = None,
) -> list[dict[str, Any]]:
    """Audit batch for recurring variants across distinct patients to detect artifact/cluster escalation."""
    patient_map = patient_map or {}

    variant_carriers: dict[tuple[str, str], dict[str, Any]] = {}

    for sid, variants in rows_by_sample.items():
        entity_id = patient_map.get(sid) or sid
        seen_in_sample = set()

        for v in variants:
            gene = v.get("gene", "")
            mut = v.get("mutation", "")
            suf = mut.split("@")[-1]
            frs = float(v.get("frs", 1.0))

            if gene not in CORE_AND_MODIFIERS:
                continue
            if frs < 0.10:
                continue

            flags = classify(suf)
            if flags.get("IS_SYN") or flags.get("IS_NULL"):
                continue

            v_key = (gene, suf)
            if v_key in seen_in_sample:
                continue
            seen_in_sample.add(v_key)

            if v_key not in variant_carriers:
                variant_carriers[v_key] = {"entities": set(), "sample_ids": []}
            variant_carriers[v_key]["entities"].add(entity_id)
            variant_carriers[v_key]["sample_ids"].append(sid)

    alerts = []

    for (gene, suf), data in variant_carriers.items():
        n_patients = len(data["entities"])
        if n_patients < 2:
            continue

        mut_full = f"{gene}@{suf}"
        card = evidence_store.card(mut_full)
        n_global_carriers = card.get("n_carriers", 0) if card else 0

        is_hp = _is_homopolymer_hotspot(gene, suf)

        if is_hp or n_global_carriers >= 50:
            category = "QC_CONTAMINATION_OR_ARTIFACT_RISK"
            msg = (
                f"{mut_full} appeared in {n_patients} patients in this batch at a known homopolymer "
                "or high-homoplasy site. Rule out Nanopore indel artifact or wet-lab cross-contamination "
                "before suspecting transmission."
            )
        elif n_global_carriers < 20:
            category = "RARE_VARIANT_WGS_ESCALATION"
            msg = (
                f"Rare variant {mut_full} (global cohort carriers: {n_global_carriers}) appeared in "
                f"{n_patients} distinct patients in this batch. Five-gene targeted profiles cannot confirm "
                "transmission alone—escalate to WGS (<=12-SNP clustering) and review lab/clinic links."
            )
        else:
            category = "RECURRENT_VARIANT_NOTICE"
            msg = f"Variant {mut_full} appeared in {n_patients} distinct patients in this batch."

        alerts.append({
            "category": category,
            "gene": gene,
            "mutation": suf,
            "variant": mut_full,
            "n_patients": n_patients,
            "sample_ids": sorted(set(data["sample_ids"])),
            "message": msg,
        })

    cat_order = {
        "QC_CONTAMINATION_OR_ARTIFACT_RISK": 0,
        "RARE_VARIANT_WGS_ESCALATION": 1,
        "RECURRENT_VARIANT_NOTICE": 2,
    }
    alerts.sort(key=lambda a: (cat_order.get(a["category"], 9), -a["n_patients"], a["variant"]))
    return alerts
