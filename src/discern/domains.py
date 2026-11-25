"""Verified H37Rv protein domain architecture for all five BDQ-relevant loci.

Static offline reference metadata sourced from UniProtKB & InterPro
(Mycobacterium tuberculosis strain ATCC 25618 / H37Rv, taxon 83332):
  - Rv0678 (mmpR5): UniProt I6Y8F7 (165 aa), mechanism.py, structure.py
  - atpE (Rv1305): UniProt P9WPS1 (81 aa), Pfam PF00137
  - pepQ (Rv2535c): UniProt I6YDN6 (372 aa), Pfam PF01321 / PF00557
  - mmpL5 (Rv0676c): UniProt P9WJV1 (964 aa), InterPro IPR004869 / Pfam PF03176
  - mmpS5 (Rv0677c): UniProt P9WJS7 (142 aa), InterPro IPR008693 / Pfam PF05423

License references for bundled static annotations:
  - UniProt: https://www.uniprot.org/help/license
  - InterPro: https://www.ebi.ac.uk/about/terms-of-use/
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from .grammar import classify

_TRACK_X0 = 24.0
_TRACK_W = 592.0


@dataclass(frozen=True)
class DomainSegment:
    start: int
    end: int
    name: str
    short_label: str
    fill: str
    text_fill: str
    role: str
    is_chip: bool = True


@dataclass(frozen=True)
class GeneArchitecture:
    gene: str
    locus: str
    uniprot: str
    length: int
    summary: str
    has_model_rules: bool
    source_label: str
    segments: tuple[DomainSegment, ...]
    ticks: tuple[int, ...]
    accent: tuple[int, int, str, str] | None = None  # (start, end, color, label)


ARCHITECTURES: dict[str, GeneArchitecture] = {
    "Rv0678": GeneArchitecture(
        gene="Rv0678",
        locus="mmpR5",
        uniprot="I6Y8F7",
        length=165,
        summary="MarR-family transcriptional repressor of mmpS5–mmpL5",
        has_model_rules=True,
        source_label="Domain rules active (DBD ≤99, late >160)",
        segments=(
            DomainSegment(1, 33, "N-term", "N-term (1–33)", "#cbd5e1", "#0f172a", "Fold initiation"),
            DomainSegment(34, 99, "DNA-Binding", "DNA-Binding (34–99)", "#1e3a5f", "#ffffff", "Operator DNA binding · α4 helix (62–68)"),
            DomainSegment(100, 160, "Dimerization", "Dimerization (100–160)", "#0f766e", "#ffffff", "Homodimer formation for pump repression"),
            DomainSegment(161, 165, "C-term", "", "#475569", "#ffffff", "C-terminal tail", is_chip=False),
        ),
        ticks=(1, 34, 100, 160, 165),
        accent=(62, 68, "#38bdf8", "α4 helix"),
    ),
    "atpE": GeneArchitecture(
        gene="atpE",
        locus="Rv1305",
        uniprot="P9WPS1",
        length=81,
        summary="ATP synthase F0 subunit c (bedaquiline target)",
        has_model_rules=False,
        source_label="Reference structure (UniProt P9WPS1)",
        segments=(
            DomainSegment(1, 4, "N-term", "", "#cbd5e1", "#0f172a", "N-terminal tail", is_chip=False),
            DomainSegment(5, 25, "TM Helix 1", "TM Helix 1 (5–25)", "#1e3a5f", "#ffffff", "Transmembrane helix 1"),
            DomainSegment(26, 56, "Cytoplasmic Loop", "Loop (26–56)", "#94a3b8", "#0f172a", "Polar loop (D28 cleft)"),
            DomainSegment(57, 77, "TM Helix 2", "TM Helix 2 (57–77)", "#0f766e", "#ffffff", "Transmembrane helix 2 (E61 site)"),
            DomainSegment(78, 81, "C-term", "", "#475569", "#ffffff", "C-terminal tail", is_chip=False),
        ),
        ticks=(1, 5, 26, 57, 78, 81),
        accent=(61, 61, "#38bdf8", "E61 proton site"),
    ),
    "pepQ": GeneArchitecture(
        gene="pepQ",
        locus="Rv2535c",
        uniprot="I6YDN6",
        length=372,
        summary="Cytoplasmic proline peptidase (low-level BDQ modifier)",
        has_model_rules=False,
        source_label="Reference structure (UniProt I6YDN6 / Pfam)",
        segments=(
            DomainSegment(1, 5, "N-term", "", "#cbd5e1", "#0f172a", "N-terminus", is_chip=False),
            DomainSegment(6, 134, "Creatinase N-term", "Creatinase N-term (6–134)", "#1e3a5f", "#ffffff", "Pfam PF01321 N-terminal domain"),
            DomainSegment(135, 141, "Linker", "", "#94a3b8", "#0f172a", "Inter-domain linker"),
            DomainSegment(142, 344, "Peptidase M24", "Peptidase M24 (142–344)", "#0f766e", "#ffffff", "Pfam PF00557 catalytic domain"),
            DomainSegment(345, 372, "C-term", "", "#475569", "#ffffff", "C-terminal tail", is_chip=False),
        ),
        ticks=(1, 6, 134, 142, 344, 372),
    ),
    "mmpL5": GeneArchitecture(
        gene="mmpL5",
        locus="Rv0676c",
        uniprot="P9WJV1",
        length=964,
        summary="RND efflux transporter (LoF triggers epistatic pump veto)",
        has_model_rules=False,
        source_label="Reference structure (UniProt P9WJV1 / IPR004869)",
        segments=(
            DomainSegment(1, 62, "N-term / TM1", "", "#cbd5e1", "#0f172a", "N-terminus & TM helix 1"),
            DomainSegment(63, 393, "RND TM Domain 1", "RND TM1 (63–393)", "#1e3a5f", "#ffffff", "Pfam PF03176 N-terminal TM bundle"),
            DomainSegment(394, 614, "Periplasmic Loop", "Porter Loop (394–614)", "#475569", "#ffffff", "Soluble periplasmic porter domain"),
            DomainSegment(615, 949, "RND TM Domain 2", "RND TM2 (615–949)", "#0f766e", "#ffffff", "Pfam PF03176 C-terminal TM bundle"),
            DomainSegment(950, 964, "C-term", "", "#64748b", "#ffffff", "Cytoplasmic tail", is_chip=False),
        ),
        ticks=(1, 63, 394, 615, 950, 964),
    ),
    "mmpS5": GeneArchitecture(
        gene="mmpS5",
        locus="Rv0677c",
        uniprot="P9WJS7",
        length=142,
        summary="MmpL5 accessory membrane protein (pump complex partner)",
        has_model_rules=False,
        source_label="Reference structure (UniProt P9WJS7 / PF05423)",
        segments=(
            DomainSegment(1, 6, "N-term", "", "#cbd5e1", "#0f172a", "Cytoplasmic N-terminus"),
            DomainSegment(7, 26, "TM Anchor", "TM (7–26)", "#1e3a5f", "#ffffff", "Single-pass membrane anchor"),
            DomainSegment(27, 142, "MmpS5 Periplasmic", "MmpS5 Periplasmic Domain (27–142)", "#0f766e", "#ffffff", "Pfam PF05423 MmpL5-binding domain"),
        ),
        ticks=(1, 7, 27, 142),
    ),
}


def _resolve_position(
    arch: GeneArchitecture, mutation: str, lof_annotation: dict | None,
) -> tuple[int, bool, str, str] | None:
    """Return (aa_pos, is_trunc, badge_label, trunc_noun) or None if non-coding."""
    flags = classify(mutation)
    if flags["IS_SYN"] or flags["IS_NULL"] or flags["is_prom"]:
        return None

    stop = lof_annotation.get("first_stop_aa") if lof_annotation else None
    if isinstance(stop, (int, float)) and stop > 0:
        aa = int(stop)
        if 1 <= aa <= arch.length:
            return aa, True, f"STOP \u00b7 aa {aa}", "Truncated"

    if flags["is_stop"]:
        m = re.search(r"(\d+)", mutation)
        if m:
            aa = int(m.group(1))
            if 1 <= aa <= arch.length:
                return aa, True, f"STOP \u00b7 aa {aa}", "Truncated"
        return None

    if flags["is_mis"] and flags["pos"] is not None:
        aa = int(flags["pos"])
        if 1 <= aa <= arch.length:
            return aa, False, f"VARIANT \u00b7 aa {aa}", "Substitution"
        return None

    if flags["is_indel"]:
        suf = mutation.split("@")[-1]
        m = re.match(r"^(-?\d+)", suf)
        if not m:
            return None
        nt_pos = int(m.group(1))
        if nt_pos <= 0:
            return None
        aa = min(arch.length, max(1, ((nt_pos - 1) // 3) + 1))
        # Check whether indel length is a multiple of 3 (in-frame vs frameshift)
        indel_m = re.search(r"_(?:ins|del)_([a-zA-Z]+)$", suf)
        is_inframe = bool(indel_m and len(indel_m.group(1)) % 3 == 0)
        if is_inframe:
            return aa, False, f"IN-FRAME \u00b7 aa {aa}", "In-frame indel"
        return aa, True, f"FRAMESHIFT \u00b7 aa {aa}", "Frameshift"

    return None


def _build_promoter_track(arch: GeneArchitecture, mutation: str, flags: dict) -> dict | None:
    """Build an upstream promoter / start-junction track (-40 to +20 nt) for negative-position variants."""
    suf = mutation.split("@")[-1]
    m = re.match(r"^[acgt]?(-?\d+)", suf, re.IGNORECASE)
    if not m:
        return None
    nt_pos = int(m.group(1))
    if nt_pos > 0:
        return None

    span_len = 1
    is_del = False
    if flags["is_indel"]:
        im = re.search(r"_(ins|del)_([a-zA-Z]+)$", suf)
        if im:
            is_del = im.group(1).lower() == "del"
            span_len = len(im.group(2)) if is_del else 1
    nt_end = nt_pos + span_len - 1 if is_del else nt_pos
    crosses_start = nt_end >= 1

    # Map window [-40, +20] (60 nt) to [24, 616]
    win_min, win_max = -40.0, 20.0
    scale = _TRACK_W / (win_max - win_min)

    def _x(pos: float) -> float:
        clamped = max(win_min, min(win_max, pos))
        return _TRACK_X0 + (clamped - win_min) * scale

    x_start_codon = _x(0)  # boundary between -1 and +1
    x_cds_3 = _x(3)
    pin_x = _x(nt_pos)
    span_x = _x(nt_pos)
    span_w = max(4.0, _x(nt_end + 1) - _x(nt_pos)) if is_del and span_len > 1 else 0.0
    badge_x = max(84.0, min(556.0, pin_x + span_w / 2.0))

    if crosses_start:
        badge_label = f"START CODON DEL \u00b7 {nt_pos}..+{nt_end}"
        right_pill = f"Crosses start codon ({nt_pos} to +{nt_end} nt)"
        pill_cls = "pill-trunc"
    else:
        badge_label = f"PROMOTER \u00b7 nt {nt_pos}"
        right_pill = f"Upstream 5\u2032 region (nt {nt_pos})"
        pill_cls = "pill-mis"

    ticks = []
    for t_nt, label in [(-40, "-40"), (-20, "-20"), (0, "+1 ATG/GTG"), (20, "+20 nt")]:
        tx = _x(t_nt)
        if abs(tx - pin_x) < 20.0 and t_nt not in (-40, 20):
            continue
        anchor = "middle" if t_nt not in (-40, 20) else ("start" if t_nt == -40 else "end")
        ticks.append({"label": label, "x": round(tx, 1), "anchor": anchor})

    chips = [
        {
            "label": "Upstream Promoter (\u221240 to \u22121)",
            "state": f"{'Deletion' if is_del else 'Variant'} at nt {nt_pos}",
            "css_class": "chip-hit",
        },
        {
            "label": "Translation Start (+1..+3)",
            "state": "Start codon disrupted" if crosses_start else "Intact",
            "css_class": "chip-lost" if crosses_start else "chip-intact",
        },
    ]

    return {
        "kind": "promoter",
        "gene": arch.gene,
        "locus": arch.locus,
        "length": arch.length,
        "source_label": "Upstream promoter & start junction",
        "right_pill": right_pill,
        "pill_cls": pill_cls,
        "badge_label": badge_label,
        "badge_x": round(badge_x, 1),
        "pin_x": round(pin_x, 1),
        "span_x": round(span_x, 1),
        "span_w": round(span_w, 1),
        "x_start_codon": round(x_start_codon, 1),
        "x_cds_3": round(x_cds_3, 1),
        "ticks": ticks,
        "chips": chips,
        "nt_pos": nt_pos,
    }


def build_domain_track(
    gene: str, mutation: str, lof_annotation: dict | None = None,
) -> dict | None:
    """Build SVG coordinates and domain impact chips for a gene variant card."""
    arch = ARCHITECTURES.get(gene)
    if arch is None:
        return None
    flags = classify(mutation)
    if flags["IS_SYN"] or flags["IS_NULL"]:
        return None
    if flags["is_prom"]:
        return _build_promoter_track(arch, mutation, flags)

    resolved = _resolve_position(arch, mutation, lof_annotation)
    if resolved is None:
        if flags["is_indel"]:
            return _build_promoter_track(arch, mutation, flags)
        return None

    aa_pos, is_trunc, badge_label, trunc_noun = resolved
    scale = _TRACK_W / float(arch.length)
    pin_x = _TRACK_X0 + (aa_pos * scale)
    badge_x = max(84.0, min(556.0, pin_x))
    retained_w = (pin_x - _TRACK_X0) if is_trunc else _TRACK_W
    retained_pct = round((aa_pos / float(arch.length)) * 100)

    segments = []
    for seg in arch.segments:
        x = _TRACK_X0 + ((seg.start - 1) * scale)
        w = (seg.end - seg.start + 1) * scale
        show_label = bool(seg.short_label) and w >= 72.0
        segments.append({
            "start": seg.start,
            "end": seg.end,
            "name": seg.name,
            "label": seg.short_label if show_label else "",
            "x": round(x, 1),
            "w": round(w, 1),
            "mid_x": round(x + w / 2.0, 1),
            "fill": seg.fill,
            "text_fill": seg.text_fill,
        })

    accent = None
    if arch.accent is not None:
        a_start, a_end, a_color, a_desc = arch.accent
        ax = _TRACK_X0 + ((a_start - 1) * scale)
        aw = max(4.0, (a_end - a_start + 1) * scale)
        accent = {"x": round(ax, 1), "w": round(aw, 1), "fill": a_color, "label": a_desc}

    # Filter ticks to prevent collisions with pin stem or adjacent ticks
    ticks = []
    last_tick_x = -999.0
    for idx, t_aa in enumerate(arch.ticks):
        if t_aa == 1:
            tx = _TRACK_X0
            anchor = "middle"
        elif t_aa == arch.length:
            tx = _TRACK_X0 + _TRACK_W
            anchor = "start"
        else:
            tx = _TRACK_X0 + ((t_aa - 1) * scale)
            anchor = "end" if idx == len(arch.ticks) - 2 and (_TRACK_X0 + _TRACK_W - tx) < 26 else "middle"
        if abs(tx - pin_x) < 22.0 or abs(tx - last_tick_x) < 22.0:
            continue
        last_tick_x = tx
        ticks.append({"aa": t_aa, "x": round(tx, 1), "anchor": anchor})

    chips = []
    for seg in arch.segments:
        if not seg.is_chip:
            continue
        if is_trunc:
            if aa_pos < seg.start:
                state = "Lost downstream"
                css_class = "chip-lost"
            elif seg.start <= aa_pos <= seg.end:
                state = f"{trunc_noun} at aa {aa_pos}"
                css_class = "chip-disrupted"
            else:
                state = "Preserved" if (gene == "Rv0678" and aa_pos > 160 and seg.end == 160) else "Intact"
                css_class = "chip-intact"
        else:
            if seg.start <= aa_pos <= seg.end:
                note = " \u00b7 \u03b14 helix" if (gene == "Rv0678" and 62 <= aa_pos <= 68) else ""
                state = f"{trunc_noun} (aa {aa_pos}{note})"
                css_class = "chip-hit"
            else:
                state = "Intact"
                css_class = "chip-intact"
        chips.append({
            "label": f"{seg.name} ({seg.start}\u2013{seg.end})" if len(arch.segments) != 4 or sum(1 for s in arch.segments if s.is_chip) != 3 else f"{seg.name} {seg.start}\u2013{seg.end}",
            "state": state,
            "css_class": css_class,
        })

    return {
        "kind": "protein",
        "gene": arch.gene,
        "locus": arch.locus,
        "uniprot": arch.uniprot,
        "length": arch.length,
        "summary": arch.summary,
        "has_model_rules": arch.has_model_rules,
        "source_label": arch.source_label,
        "aa_pos": aa_pos,
        "is_trunc": is_trunc,
        "badge_label": badge_label,
        "pin_x": round(pin_x, 1),
        "badge_x": round(badge_x, 1),
        "retained_w": round(retained_w, 1),
        "retained_pct": retained_pct,
        "segments": segments,
        "accent": accent,
        "ticks": ticks,
        "chips": chips,
    }


def attach_domain_tracks(record: dict) -> dict:
    """Return a shallow copy of a result record with domain_track populated per variant card."""
    evidence = record.get("evidence")
    if not evidence:
        return record
    enriched = []
    for item in evidence:
        card = dict(item)
        variant = card.get("variant", "")
        gene, _, mutation = variant.partition("@")
        if gene and mutation and "domain_track" not in card:
            card["domain_track"] = build_domain_track(gene, mutation, card.get("lof_annotation"))
        enriched.append(card)
    return {**record, "evidence": enriched}
