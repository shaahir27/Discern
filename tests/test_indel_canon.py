"""Tests for sequence-based matching of indel names (GARC vs HGVS)."""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest

from discern.indel_canon import (
    RV0678_CDS as CDS,
)
from discern.indel_canon import (
    IndelIndex,
    apply_garc,
    del_matches_reference,
    hgvs_indel_to_garc,
    left_align,
)


def test_rv0678_cds_is_clean_orf():
    bases = "TCAG"
    aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    tab = {a + b + c: aas[i] for i, (a, b, c) in enumerate((x, y, z) for x in bases for y in bases for z in bases)}
    prot = "".join(tab[CDS[i:i + 3]] for i in range(0, len(CDS), 3))
    assert len(CDS) == 498
    assert prot.endswith("*")
    assert "*" not in prot[:-1]
    assert CDS.startswith("GTG")


def test_training_indel_canon_consistency():
    path = Path("data/raw/tb-bdq-cat-main/data/MUTATIONS_BDQ.csv.gz")
    if not path.is_file():
        pytest.skip("MUTATIONS_BDQ.csv.gz not found")
    m = pd.read_csv(path)
    names = pd.concat([m[m.GENE == "Rv0678"].MUTATION, m[m.GENE == "Rv0678"].MINOR_MUTATION]).dropna().astype(str)
    plain = [n for n in names.unique() if re.match(r"^\d+_(ins|del)_[acgt]+$", n)]
    counts = names.value_counts().to_dict()

    dels = [n for n in plain if "_del_" in n]
    bad = [n for n in dels if not del_matches_reference(CDS, n)]
    assert not bad, f"deletion convention mismatch: {bad[:5]}"

    idx = IndelIndex(CDS, plain, counts)

    def respellings(name):
        out, seen = [name], {name}
        m_ins = re.match(r"^(\d+)_ins_([acgt]+)$", name)
        if m_ins:
            p, b = int(m_ins.group(1)), m_ins.group(2).upper()
            while p < len(CDS) and CDS[p] == b[0]:
                b, p = b[1:] + b[0], p + 1
                nm = f"{p}_ins_{b.lower()}"
                if nm not in seen:
                    out.append(nm)
                    seen.add(nm)
        m_del = re.match(r"^(\d+)_del_([acgt]+)$", name)
        if m_del:
            p, n = int(m_del.group(1)), len(m_del.group(2))
            while p - 1 + n < len(CDS) and CDS[p - 1] == CDS[p - 1 + n]:
                p += 1
                nm = f"{p}_del_{CDS[p - 1:p - 1 + n].lower()}"
                if nm not in seen:
                    out.append(nm)
                    seen.add(nm)
        return out

    spell_groups = {}
    for n in plain:
        spell_groups.setdefault(apply_garc(CDS, n), []).append(n)

    for n in plain:
        alts = respellings(n)
        best = max(spell_groups[apply_garc(CDS, n)], key=lambda x: counts.get(x, 0))
        for a in alts:
            got, how = idx.to_training_name(a)
            assert got == best, f"{a} -> {got} ({how}); expected {best}"
            assert apply_garc(CDS, a) == apply_garc(CDS, n)


def test_hgvs_known_spellings_to_training_variants():
    idx = IndelIndex(
        CDS,
        ["141_ins_c", "192_ins_g", "138_ins_g"],
        {"141_ins_c": 168, "192_ins_g": 99, "138_ins_g": 58},
    )
    cases = [
        ("c.144dupC", "141_ins_c"),
        ("c.198dupG", "192_ins_g"),
        ("c.139dupG", "138_ins_g"),
    ]
    for hgvs, want in cases:
        g = hgvs_indel_to_garc(CDS, hgvs)
        assert g is not None
        got, how = idx.to_training_name(g)
        assert got == want
        assert how in ("exact", "sequence_match")

    assert hgvs_indel_to_garc(CDS, "c.209del") == "209_del_a"
    assert left_align(CDS, "198_ins_g") == "192_ins_g"
