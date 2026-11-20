"""Sequencing-assay and tNGS panel profile definitions and blind-spot auditing."""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Literal

from .grammar import classify

ALL_BDQ_GENES = ["Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5"]
CORE_BDQ_GENES = ["Rv0678", "atpE", "pepQ"]


@dataclass(frozen=True)
class PanelProfile:
    key: str
    label: str
    platform: Literal["illumina", "nanopore", "mixed"]
    covered_genes: frozenset[str]
    vaf_floor: float
    notes: str


PANELS: dict[str, PanelProfile] = {
    "wgs": PanelProfile(
        key="wgs",
        label="Whole Genome Sequencing (WGS) — All 5 loci",
        platform="mixed",
        covered_genes=frozenset({"Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5"}),
        vaf_floor=0.10,
        notes="Comprehensive 5-locus coverage (default).",
    ),
    "deeplex_myc_tb": PanelProfile(
        key="deeplex_myc_tb",
        label="GenoScreen Deeplex Myc-TB (Commercial tNGS)",
        platform="illumina",
        covered_genes=frozenset({"Rv0678"}),
        vaf_floor=0.03,
        notes="Targeted Illumina panel; sequences Rv0678 only (missing atpE, pepQ, mmpL5, mmpS5).",
    ),
    "deeplex_xl": PanelProfile(
        key="deeplex_xl",
        label="Deeplex Myc-TB XL (Expanded beta tNGS)",
        platform="illumina",
        covered_genes=frozenset({"Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5"}),
        vaf_floor=0.03,
        notes="Expanded research panel covering all 5 bedaquiline resistance loci.",
    ),
    "ampore_tb": PanelProfile(
        key="ampore_tb",
        label="ONT AmPORE-TB / wf-tb-amr (Nanopore)",
        platform="nanopore",
        covered_genes=frozenset({"Rv0678", "atpE"}),
        vaf_floor=0.15,
        notes="Targeted Oxford Nanopore amplicon panel; missing pepQ and mmpL5/mmpS5.",
    ),
    "icmr_nirt_minion": PanelProfile(
        key="icmr_nirt_minion",
        label="ICMR-NIRT India MinION Pilot",
        platform="nanopore",
        covered_genes=frozenset({"Rv0678"}),
        vaf_floor=0.15,
        notes="Targeted MinION sequencing for Indian surveillance; sequences Rv0678 only.",
    ),
}


def _parse_indel_pos(suffix: str) -> int | None:
    """Extract 1-based nucleotide position from an indel mutation suffix."""
    m = re.match(r"^(-?\d+)", suffix)
    if m:
        try:
            return int(m.group(1))
        except ValueError:
            return None
    return None


def audit_sample_panel(variants: list[dict], assay_key: str = "wgs") -> dict:
    """Audit a sample's variant list against a specific sequencing assay/panel."""
    panel = PANELS.get(assay_key, PANELS["wgs"])
    warnings: list[str] = []

    if assay_key not in PANELS:
        warnings.append("UNKNOWN_ASSAY_FALLBACK_WGS")

    # Determine per-gene mutation presence
    gene_mutated: dict[str, bool] = {g: False for g in ALL_BDQ_GENES}
    has_homopolymer_risk = False
    has_below_lod = False

    for v in variants:
        gene = v.get("gene", "")
        mut = v.get("mutation", "")
        suffix = mut.split("@")[-1]
        frs = float(v.get("frs", 1.0))

        flags = classify(suffix)
        is_retained = (not flags["IS_SYN"]) and (not flags["IS_NULL"]) and (frs >= 0.10 or frs >= panel.vaf_floor)

        if gene in ALL_BDQ_GENES and is_retained and gene in panel.covered_genes:
            gene_mutated[gene] = True

        # Check LOD floor
        if 0 < frs < panel.vaf_floor:
            has_below_lod = True

        # Check homopolymer artifact risk on Rv0678
        if gene == "Rv0678" and flags.get("is_indel", False):
            pos = _parse_indel_pos(suffix)
            if pos is not None and ((192 <= pos <= 198) or (138 <= pos <= 144)) and (frs < 0.25 or panel.platform == "nanopore"):
                has_homopolymer_risk = True

    # Gene status
    gene_status: dict[str, str] = {}
    for g in ALL_BDQ_GENES:
        if g not in panel.covered_genes:
            gene_status[g] = "not_interrogated"
        elif gene_mutated[g]:
            gene_status[g] = "mutated"
        else:
            gene_status[g] = "wild_type"

    # Epistatic veto status
    has_pump_genes = frozenset({"mmpL5", "mmpS5"}).issubset(panel.covered_genes)
    epistatic_veto_status = "active" if has_pump_genes else "blind_missing_pump_genes"

    if epistatic_veto_status != "active":
        warnings.append("EPISTATIC_VETO_BLIND")

    if not frozenset({"atpE", "pepQ"}).issubset(panel.covered_genes):
        warnings.append("UNSEQUENCED_CORE_LOCI")

    if has_homopolymer_risk:
        warnings.append("HOMOPOLYMER_ARTIFACT_RISK")

    if has_below_lod:
        warnings.append("BELOW_PANEL_LOD")

    return {
        "assay_key": panel.key,
        "assay_label": panel.label,
        "platform": panel.platform,
        "vaf_floor": panel.vaf_floor,
        "covered_genes": sorted(panel.covered_genes),
        "unsequenced_genes": sorted(set(ALL_BDQ_GENES) - panel.covered_genes),
        "gene_status": gene_status,
        "epistatic_veto_status": epistatic_veto_status,
        "warnings": warnings,
        "notes": panel.notes,
    }


def compute_cohort_panel_coverage(covered_genes: frozenset[str]) -> dict:
    """Deterministic arithmetic of known-locus resistant cohort coverage.
    
    In the 11,868 reference cohort, 866 isolates are phenotypically resistant:
    - 729 isolates carry identified resistance variants in core genes:
        Rv0678 carriers: 645 isolates (88.48% of known-locus R)
        atpE carriers:     51 isolates ( 7.00% of known-locus R)
        pepQ carriers:     42 isolates ( 5.76% of known-locus R)
        (9 isolates carry co-occurring core variants)
    - 137 isolates (15.8%) carry no core-gene mutations (invisible resistance).
    """
    if len(covered_genes) == 0:
        return {
            "known_locus_sensitivity": 0.0,
            "overall_sensitivity_upper_bound": 0.0,
            "missed_known_pct": 100.0,
            "missed_loci": ["Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5"],
        }

    # Known-locus counts (non-overlapping decomposition)
    n_known_total = 729
    n_rv_only = 636
    n_atpe_only = 44
    n_pepq_only = 40
    n_multi = 9

    covered_count = 0
    if "Rv0678" in covered_genes:
        covered_count += n_rv_only
    if "atpE" in covered_genes:
        covered_count += n_atpe_only
    if "pepQ" in covered_genes:
        covered_count += n_pepq_only
    if any(g in covered_genes for g in ["Rv0678", "atpE", "pepQ"]):
        covered_count += n_multi

    sens_known = round(covered_count / n_known_total, 4)
    # Overall sensitivity upper bound among all 866 resistant strains
    sens_overall = round(covered_count / 866, 4)
    missed_pct = round((1.0 - sens_known) * 100, 2)
    missed_loci = sorted(set(ALL_BDQ_GENES) - covered_genes)

    return {
        "known_locus_sensitivity": sens_known,
        "overall_sensitivity_upper_bound": sens_overall,
        "missed_known_pct": missed_pct,
        "missed_loci": missed_loci,
    }
