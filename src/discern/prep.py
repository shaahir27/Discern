"""Data preparation (AGENTS.md section 5). Golden implementation."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .grammar import classify

SAMPLE_COLS = ["UNIQUEID", "DRUG", "SOURCE", "METHOD_3", "METHOD_MIC", "PHENOTYPE", "QUALITY", "SRC"]


def load_samples(raw: Path) -> pd.DataFrame:
    cry = pd.read_csv(raw / "data/DST_MEASUREMENTS_BDQ_CRYPTIC.csv.gz", low_memory=False)
    cry = cry[cry.QUALITY == "HIGH"].copy()
    cry["SRC"] = "CRyPTIC"
    nic = pd.read_csv(raw / "data/DST_MEASUREMENTS_BDQ_NICD.csv.gz", low_memory=False)
    nic["SRC"] = "NICD"
    s = pd.concat([cry[SAMPLE_COLS], nic[SAMPLE_COLS]], ignore_index=True)
    s["_r"] = (s.PHENOTYPE == "R").astype(int)
    s["_m"] = s.METHOD_MIC.notna().astype(int)
    s = (
        s.sort_values(["UNIQUEID", "_r", "_m"], ascending=[True, False, False], kind="stable")
        .drop_duplicates("UNIQUEID", keep="first")
        .drop(columns=["_r", "_m"])
    )
    g = pd.read_csv(raw / "data/GENOMES.csv.gz", low_memory=False)[["UNIQUEID", "SITEID", "SUBJID"]]
    g = g.drop_duplicates("UNIQUEID")
    s = s.merge(g, on="UNIQUEID", how="inner").reset_index(drop=True)
    s["y"] = (s.PHENOTYPE == "R").astype(int)
    return s


def load_mutations(raw: Path, ids) -> pd.DataFrame:
    m = pd.read_csv(raw / "data/MUTATIONS_BDQ.csv.gz", low_memory=False)
    m = m[m.UNIQUEID.isin(set(ids))].copy()
    minor = m.IS_MINOR_ALLELE.astype(bool)
    m["MUT"] = m.GENE.astype(str) + "@" + np.where(minor, m.MINOR_MUTATION.astype(str), m.MUTATION.astype(str))
    m["FRSf"] = np.where(minor, m.FRS, 1.0)
    m["suffix"] = m.MUT.str.split("@").str[-1]
    flags = m.suffix.map(classify)
    for key in ("IS_SYN", "IS_NULL", "is_stop", "is_indel", "is_prom", "is_mis", "pos"):
        m[key] = flags.map(lambda d: d[key])
    return m.reset_index(drop=True)


def prepare(raw: Path, processed: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    processed.mkdir(parents=True, exist_ok=True)
    s = load_samples(raw)
    m = load_mutations(raw, s.UNIQUEID)
    s.to_parquet(processed / "samples.parquet")
    m.to_parquet(processed / "mutations.parquet")
    return s, m


def e01_counts(s: pd.DataFrame, m: pd.DataFrame) -> dict:
    minor = m.IS_MINOR_ALLELE.astype(bool)
    nul = m[m.IS_NULL]
    return {
        "n_samples": int(len(s)),
        "n_resistant": int(s.y.sum()),
        "cryptic_n": int((s.SRC == "CRyPTIC").sum()),
        "cryptic_R": int(s[s.SRC == "CRyPTIC"].y.sum()),
        "nicd_n": int((s.SRC == "NICD").sum()),
        "nicd_R": int(s[s.SRC == "NICD"].y.sum()),
        "unique_subjid": int(s.SUBJID.nunique()),
        "mutation_rows": int(len(m)),
        "minor_allele_rows": int(minor.sum()),
        "median_FRSf_minor": float(m[minor].FRSf.median()),
        "null_rows": int(len(nul)),
        "null_rows_samples": int(nul.UNIQUEID.nunique()),
    }


E01_EXPECTED = {
    "n_samples": 11868, "n_resistant": 866, "cryptic_n": 10704, "cryptic_R": 83,
    "nicd_n": 1164, "nicd_R": 783, "unique_subjid": 11611, "mutation_rows": 25554,
    "minor_allele_rows": 1115, "null_rows": 11, "null_rows_samples": 4,
}