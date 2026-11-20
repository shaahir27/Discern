"""G2 / P-E16: Rv0678 loss-of-function mechanism biology (Snobre 2024 reasoning).

Classifies truncating Rv0678 variants by first-stop position using the reference CDS
(NC_000962.3:778,990-779,487, + strand, GTG start, 165 aa). Regions: DBD <=99, mid 100-160,
late >160. Late stops preserve the dimerization domain (function-preserving proxy for Snobre's
conserved alternative-reading-frame mechanism). A full phylogenetic ARF-conservation alignment is
out of scope; late-stop preservation is the dominant, exactly-computable function-preserving case.

HIGH-confidence LoF = stop / frameshift / whole-gene deletion, EXCLUDING single-nt homopolymer
indels (sequencing-artifact prone). The epistatic veto downgrades R->U only when a high-confidence
Rv0678 LoF co-occurs with mmpL5/S5 LoF and no atpE variant is present (never on atpE-driven R).
"""
from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AAS[i * 16 + j * 4 + k]
         for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}
DBD_MAX = 99
LATE_MIN = 161  # first-stop aa > 160 => late (function-preserving)


def load_cds(path: Path) -> str:
    return "".join(l.strip() for l in Path(path).read_text().splitlines() if not l.startswith(">"))


def translate(seq: str, start_to_m: bool = True) -> str:
    prot = []
    for i in range(0, len(seq) - 2, 3):
        aa = CODON.get(seq[i:i + 3], "X")
        if i == 0 and start_to_m and seq[i:i + 3] in ("GTG", "TTG"):
            aa = "M"
        prot.append(aa)
        if aa == "*":
            break
    return "".join(prot)


def first_stop_aa(seq: str) -> int:
    """1-based amino-acid position of the first stop codon (166 for full-length WT)."""
    prot = translate(seq, start_to_m=False)
    return prot.index("*") + 1 if "*" in prot else len(seq) // 3 + 1


def region(aa_pos: int) -> str:
    if aa_pos <= 1:
        return "full_del"
    if aa_pos <= DBD_MAX:
        return "DBD"
    if aa_pos <= 160:
        return "mid"
    return "late"


def _is_homopolymer(cds: str, pos: int, length: int, nt: str) -> bool:
    """Single-nt indel inside a run of >=3 identical bases (artifact-prone). pos is 1-based."""
    if abs(length) != 1:
        return False
    p = int(pos)
    if length > 0:  # insertion of nt after position p
        b = nt.lower()
        left = cds[max(0, p - 2):p]
        right = cds[p:p + 2]
        return (len(left) == 2 and left == b * 2) or (len(right) == 2 and right == b * 2) \
            or (p >= 1 and p < len(cds) and cds[p - 1] == b and cds[p] == b)
    b = cds[p - 1].lower() if 1 <= p <= len(cds) else ""
    i = p - 1
    run = 1
    j = i - 1
    while j >= 0 and cds[j].lower() == b:
        run += 1; j -= 1
    j = i + 1
    while j < len(cds) and cds[j].lower() == b:
        run += 1; j += 1
    return run >= 3


def apply_indel(cds: str, pos: float, length: float, nt: str):
    """Return mutated CDS, or a sentinel: 'full_del' (whole-gene), 'promoter' (pos<0), None (inapplicable)."""
    if pos is None or (isinstance(pos, float) and np.isnan(pos)):
        return "full_del" if length <= -400 else None
    p = int(pos)
    if p < 0:
        return "promoter"
    L = int(length)
    if L <= -400:
        return "full_del"
    if L > 0:  # insertion after nt p
        return cds[:p] + str(nt).lower() + cds[p:]
    if L < 0:  # deletion of |L| nt starting at nt p
        return cds[:p - 1] + cds[p - 1 + abs(L):]
    return None


def classify_truncating(cds: str, row) -> dict | None:
    """Classify a single Rv0678 stop/indel variant. None if not truncating (e.g. promoter indel)."""
    suffix = row["suffix"]
    if row["is_stop"]:
        fs = int(row["AMINO_ACID_NUMBER"])
        return {"type": "stop", "first_stop_aa": fs, "region": region(fs),
                "high_conf": True, "late_preserved": fs > 160}
    if row["is_indel"]:
        pos = row["GENE_POSITION"]
        length = row["INDEL_LENGTH"]
        nt = row["INDEL_NUCLEOTIDES"]
        suf = str(suffix)
        if pd.isna(length):
            ntl = len(str(nt)) if pd.notna(nt) else 1
            length = ntl if "_ins_" in suf else -ntl
        if pd.isna(pos):
            mp = re.match(r"^(-?\d+)", suf)
            pos = float(mp.group(1)) if mp else np.nan
        if pd.isna(nt):
            nt = ""
        new = apply_indel(cds, pos, length, nt)
        if new == "promoter" or new is None:
            return None
        if new == "full_del":
            return {"type": "full_del", "first_stop_aa": 1, "region": "full_del",
                    "high_conf": True, "late_preserved": False}
        fs = first_stop_aa(new)
        inframe = (int(length) % 3 == 0)
        homo = _is_homopolymer(cds, pos, length, nt)
        return {"type": "inframe_indel" if inframe else "frameshift",
                "first_stop_aa": fs, "region": region(fs),
                "high_conf": (not inframe) and (not homo),
                "late_preserved": fs > 160}
    return None


def variant_effects(m: pd.DataFrame, cds: str) -> pd.DataFrame:
    """Per distinct Rv0678 truncating variant: mechanism class + carrier/resistant counts."""
    rv = m[(m.GENE == "Rv0678") & (m.is_stop | m.is_indel) & (~m.IS_SYN) & (~m.IS_NULL)]
    recs = {}
    for _, row in rv.drop_duplicates("MUT").iterrows():
        c = classify_truncating(cds, row)
        if c is not None:
            recs[row["MUT"]] = c
    eff = pd.DataFrame([{"MUT": k, **v} for k, v in recs.items()])
    if eff.empty:
        return eff
    carriers = rv.groupby("MUT").UNIQUEID.nunique().rename("n_carriers")
    res = rv[rv.y == 1].groupby("MUT").UNIQUEID.nunique().rename("n_resistant") if "y" in rv else None
    eff = eff.merge(carriers, on="MUT", how="left")
    if res is not None:
        eff = eff.merge(res, on="MUT", how="left").fillna({"n_resistant": 0})
        eff["frac_resistant"] = eff.n_resistant / eff.n_carriers
    return eff.sort_values(["region", "first_stop_aa"]).reset_index(drop=True)


def veto_flags(
    retained: pd.DataFrame, ids: Sequence[str] | np.ndarray, high_conf: set[str],
    tolerance: float = 0.20,
) -> pd.DataFrame:
    """Evaluate the epistatic pump veto with allele-frequency linkage checking.

    Under WHO Catalogue 3rd Edition (2026) and TB-Profiler v6.6.5, an mmpL5/mmpS5 LoF
    must be genetically linked to the Rv0678 variant to inactivate efflux:
        FRS(mmpL5/S5) >= FRS(Rv0678) - tolerance
    A minor subclone pump LoF (e.g. 10%) cannot veto a dominant (e.g. 95%) Rv0678 LoF.
    """
    rv = retained[(retained.GENE == "Rv0678") & retained.MUT.isin(high_conf) & (retained.FRSf > 0)]
    pump = retained[retained.GENE.isin(["mmpL5", "mmpS5"]) & (retained.is_indel | retained.is_stop) & (retained.FRSf > 0)]
    atpe = retained[retained.GENE == "atpE"]
    flags = pd.DataFrame(index=ids)
    flags["rv_highconf_lof"] = flags.index.isin(rv.UNIQUEID)
    flags["pump_lof"] = flags.index.isin(pump.UNIQUEID)
    flags["atpE_present"] = flags.index.isin(atpe.UNIQUEID)

    if not rv.empty and not pump.empty:
        rv_frs = rv.groupby("UNIQUEID").FRSf.max()
        pump_frs = pump.groupby("UNIQUEID").FRSf.max()
        s_rv = pd.Series(flags.index, index=flags.index).map(rv_frs).fillna(0.0)
        s_pump = pd.Series(flags.index, index=flags.index).map(pump_frs).fillna(0.0)
        linked = (s_pump >= (s_rv - tolerance)) & (s_pump > 0)
    else:
        linked = pd.Series(False, index=flags.index)
    flags["frs_linked"] = linked.values
    flags["veto_fire"] = flags.rv_highconf_lof & flags.pump_lof & ~flags.atpE_present & flags.frs_linked
    return flags


def sample_mech_flags(m: pd.DataFrame, ids, cds: str, frs_min: float = 0.1):
    """Per-sample mechanism features (max FRS) + epistatic-veto flags."""
    mm = m[(~m.IS_SYN) & (~m.IS_NULL) & (m.FRSf >= frs_min)]
    rv = mm[mm.GENE == "Rv0678"]
    cls = {}
    for _, row in rv.drop_duplicates("MUT").iterrows():
        c = classify_truncating(cds, row)
        if c is not None:
            cls[row["MUT"]] = c

    def mx(mask_series):
        return mask_series.groupby("UNIQUEID").FRSf.max().reindex(ids, fill_value=0.0)

    rv2 = rv.assign(
        hc=rv.MUT.map(lambda v: cls.get(v, {}).get("high_conf", False)),
        dbd=rv.MUT.map(lambda v: cls.get(v, {}).get("region") == "DBD"),
        late=rv.MUT.map(lambda v: cls.get(v, {}).get("late_preserved", False)),
        trunc=rv.MUT.map(lambda v: v in cls),
    )
    feat = pd.DataFrame(index=ids)
    feat["rv_lof_any"] = mx(rv2[rv2.trunc])
    feat["rv_lof_highconf"] = mx(rv2[rv2.hc])
    feat["rv_lof_dbd"] = mx(rv2[rv2.hc & rv2.dbd])
    feat["rv_lof_late"] = mx(rv2[rv2.hc & rv2.late])
    high_conf = {mutation for mutation, effect in cls.items() if effect["high_conf"]}
    # E16/E19 used the original stop and indel screen when producing the archive.
    historical_candidates = high_conf | set(rv.loc[rv.is_indel | rv.is_stop, "MUT"])
    return feat, veto_flags(mm, ids, historical_candidates)


def apply_veto(calls: np.ndarray, veto_fire: np.ndarray) -> np.ndarray:
    """Downgrade R->U where the epistatic veto fires (never touches atpE-driven R; see flags)."""
    out = np.array(calls, dtype=object)
    out[(out == "R") & veto_fire] = "U"
    return out
