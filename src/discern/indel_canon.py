"""Sequence-based matching of indel names (GARC vs HGVS) so the SAME DNA change always gets the SAME name.

Why: Rv0678 has runs of one repeated base (e.g. GGGGGG at 193-198). An extra G anywhere in the run is
identical DNA, but tools name it differently (training/GARC = left end: 192_ins_g; HGVS/TB-Profiler = right
end: c.198dupG). Matching by name alone makes known variants look novel. Matching by the DNA they produce does not.

Convention (VERIFIED against the Fowler training names for Rv0678):
  "N_ins_x"  = insert x AFTER coding position N (1-based, coding strand)
  "N_del_x"  = delete the len(x) bases starting AT coding position N
Only coding-sequence indels are handled. Promoter (negative position) and delins are returned unchanged.
"""
from __future__ import annotations

import re

RV0678_CDS = (
    "GTGAGCGTCAACGACGGGGTCGATCAGATGGGCGCCGAGCCCGACATCATGGAATTCGTCGAACAGATGGGCGGCTATTT"
    "CGAGTCCAGGAGTTTGACTCGGTTGGCGGGTCGATTGTTGGGCTGGCTGCTGGTGTGTGATCCCGAGCGGCAGTCCTCGG"
    "AGGAACTGGCGACGGCGCTGGCGGCCAGCAGCGGGGGGATCAGCACCAATGCCCGGATGCTGATCCAATTTGGGTTCATT"
    "GAGCGGCTCGCGGTCGCCGGGGATCGGCGCACCTATTTCCGGTTGCGGCCCAACGCTTTCGCGGCTGGCGAGCGTGAACG"
    "CATCCGGGCAATGGCCGAACTGCAGGACCTGGCTGACGTGGGGCTGAGGGCGCTGGGCGACGCCCCGCCGCAGCGAAGCC"
    "GACGGCTGCGGGAGATGCGGGATCTGTTGGCATATATGGAGAACGTCGTCTCCGACGCCCTGGGGCGATACAGCCAGCGA"
    "ACCGGAGAGGACGACTGA"
)  # 498 nt, H37Rv Rv0678 coding strand; translates to 165 aa + stop (checked)

_INS = re.compile(r"^(\d+)_ins_([acgt]+)$")
_DEL = re.compile(r"^(\d+)_del_([acgt]+)$")


def apply_garc(cds: str, name: str) -> str | None:
    """Return the mutated coding sequence, or None if the name is not a plain coding indel."""
    m = _INS.match(name)
    if m:
        p, b = int(m.group(1)), m.group(2).upper()
        return cds[:p] + b + cds[p:] if 0 < p <= len(cds) else None
    m = _DEL.match(name)
    if m:
        p, b = int(m.group(1)), m.group(2).upper()
        return cds[: p - 1] + cds[p - 1 + len(b):] if 0 < p and p - 1 + len(b) <= len(cds) else None
    return None


def del_matches_reference(cds: str, name: str) -> bool:
    """Convention check: the bases named in a deletion must equal the reference at that position."""
    m = _DEL.match(name)
    return bool(m) and cds[int(m.group(1)) - 1: int(m.group(1)) - 1 + len(m.group(2))] == m.group(2).upper()


def hgvs_indel_to_garc(cds: str, change: str) -> str | None:
    """c.139dupG / c.139dup / c.139_140dup / c.141_142insC / c.209delA / c.209del / c.193_195del -> GARC name."""
    s = change.strip()
    s = s.removeprefix("c.")
    m = re.match(r"^(\d+)(?:_(\d+))?dup([ACGTacgt]*)$", s)
    if m:
        a, e = int(m.group(1)), int(m.group(2) or m.group(1))
        b = (m.group(3) or cds[a - 1: e]).lower()
        return f"{e}_ins_{b}"
    m = re.match(r"^(\d+)_(\d+)ins([ACGTacgt]+)$", s)
    if m and int(m.group(2)) == int(m.group(1)) + 1:
        return f"{int(m.group(1))}_ins_{m.group(3).lower()}"
    m = re.match(r"^(\d+)(?:_(\d+))?del([ACGTacgt]*)$", s)
    if m:
        a, e = int(m.group(1)), int(m.group(2) or m.group(1))
        ref = cds[a - 1: e].lower()
        if m.group(3) and m.group(3).lower() != ref:
            return None  # stated bases disagree with the reference: refuse to guess
        return f"{a}_del_{ref}"
    return None


def left_align(cds: str, name: str) -> str:
    """Canonical (left-most) name for a coding indel; unchanged if not a plain coding indel."""
    m = _INS.match(name)
    if m:
        p, b = int(m.group(1)), m.group(2).upper()
        while p > 0 and cds[p - 1] == b[-1]:
            b, p = b[-1] + b[:-1], p - 1
        return f"{p}_ins_{b.lower()}"
    m = _DEL.match(name)
    if m:
        p, n = int(m.group(1)), len(m.group(2))
        while p > 1 and cds[p - 2] == cds[p - 2 + n]:
            p -= 1
        return f"{p}_del_{cds[p - 1: p - 1 + n].lower()}"
    return name


class IndelIndex:
    """Map 'mutated sequence' -> training name, so any equivalent spelling finds the training variant."""

    def __init__(self, cds: str, training_names, counts: dict | None = None):
        self.cds, self.by_seq = cds, {}
        for n in sorted(set(training_names), key=lambda x: -(counts or {}).get(x, 0)):
            s = apply_garc(cds, n)
            if s is not None:
                self.by_seq.setdefault(s, n)

    def to_training_name(self, name: str) -> tuple[str, str]:
        """Return (name_to_use, how) where how in {exact, sequence_match, left_aligned_novel, unchanged}."""
        s = apply_garc(self.cds, name)
        if s is None:
            return name, "unchanged"
        hit = self.by_seq.get(s)
        if hit is not None:
            return hit, ("exact" if hit == name else "sequence_match")
        return left_align(self.cds, name), "left_aligned_novel"
