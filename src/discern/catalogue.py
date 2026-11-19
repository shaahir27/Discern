"""Tier A catalogue backbone (AGENTS.md section 7). Golden implementation."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import piezo
from catomatic.BinaryCatalogue import BinaryBuilder

SC = ["UNIQUEID", "DRUG", "SOURCE", "METHOD_3", "PHENOTYPE", "METHOD_MIC"]
EXCLUDE = ["mmpL5", "mmpS5"]


def cat_muts(m: pd.DataFrame) -> pd.DataFrame:
    return m[~m.IS_SYN & ~m.IS_NULL & ~m.GENE.isin(EXCLUDE)]


def build_catalogue(train_s, m, path, wild_path, frs=0.1, background=0.1, tails="two", test="Binomial"):
    x = cat_muts(m[m.UNIQUEID.isin(set(train_s.UNIQUEID))])
    xb = x[["UNIQUEID", "GENE", "MUT", "FRSf"]].rename(columns={"MUT": "MUTATION", "FRSf": "FRS"})
    b = BinaryBuilder(samples=train_s[SC], mutations=xb, frs=frs)
    b.build(test=test, background=background, tails=tails, strict_unlock=True)
    wild = json.load(open(wild_path))
    df = b.build_piezo("NC_000962.3", "BDQ-cv", "1.0", "BDQ", wild, json_dumps=True)
    df.to_csv(path, index=False)
    return piezo.ResistanceCatalogue(path)


def _aggregate(series) -> str:
    if (series == "R").any():
        return "R"
    if (series == "U").any():
        return "U"
    return "S"


def apply_catalogue(rc, test_s, m, frs_min=0.1):
    x = cat_muts(m[m.UNIQUEID.isin(set(test_s.UNIQUEID)) & (m.FRSf >= frs_min)])
    cache = {}
    for v in x.MUT.unique():
        try:
            r = rc.predict(v)
            cache[v] = r.get("BDQ", "S") if isinstance(r, dict) else str(r)
        except Exception:
            cache[v] = "S"
    per = x.assign(c=x.MUT.map(cache)).groupby("UNIQUEID").c.agg(_aggregate)
    return test_s.UNIQUEID.map(per).fillna("S").values


def who_as_is(rc, s, m, thr):
    mm = m[~m.IS_SYN & ~m.IS_NULL & (m.FRSf >= thr)]
    cache = {}
    for v in mm.MUT.unique():
        try:
            r = rc.predict(v)
            cache[v] = r.get("BDQ", "S") if isinstance(r, dict) else str(r)
        except Exception:
            cache[v] = "S"
    per = mm.assign(c=mm.MUT.map(cache)).groupby("UNIQUEID").c.agg(_aggregate)
    calls = s.UNIQUEID.map(per).fillna("S").values
    return calls, cache