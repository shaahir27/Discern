"""Tier B features (AGENTS.md section 8). Golden implementation."""
from __future__ import annotations

import numpy as np
import pandas as pd

CORE = ["Rv0678", "atpE", "pepQ"]
AGG1 = ["rv_missense", "rv_mis_dbd", "rv_indel", "rv_stop", "rv_prom", "atpE_any", "pepQ_any"]
AGG2 = ["mmpL5_LoF", "rv_x_L5LoF"]


def sample_tables(m, ids, frs_min=0.1, dbd=(34, 99)):
    if m.empty:
        columns = [*AGG1, "mmpL5_LoF", "rv_any", "rv_x_L5LoF"]
        return pd.DataFrame(index=ids), pd.DataFrame(0.0, index=ids, columns=columns)

    mm = m[~m.IS_SYN & ~m.IS_NULL & (m.FRSf >= frs_min)]
    core = mm[mm.GENE.isin(CORE)].groupby(["UNIQUEID", "MUT"]).FRSf.max().unstack(fill_value=0.0)
    core = core.reindex(ids, fill_value=0.0)
    rv = mm[mm.GENE == "Rv0678"]

    def mx(df):
        return df.groupby("UNIQUEID").FRSf.max().reindex(ids, fill_value=0.0)

    a = pd.DataFrame(index=ids)
    a["rv_missense"] = mx(rv[rv.is_mis])
    a["rv_mis_dbd"] = mx(rv[rv.is_mis & rv.pos.between(dbd[0], dbd[1])])
    a["rv_indel"] = mx(rv[rv.is_indel])
    a["rv_stop"] = mx(rv[rv.is_stop])
    a["rv_prom"] = mx(rv[rv.is_prom])
    a["atpE_any"] = mx(mm[mm.GENE == "atpE"])
    a["pepQ_any"] = mx(mm[mm.GENE == "pepQ"])
    l5 = mm[(mm.GENE == "mmpL5") & (mm.is_indel | mm.is_stop)]
    a["mmpL5_LoF"] = mx(l5)
    a["rv_any"] = a[["rv_missense", "rv_indel", "rv_stop", "rv_prom"]].max(axis=1)
    a["rv_x_L5LoF"] = a.rv_any * a.mmpL5_LoF
    return core, a


def vocab_from(core, train_mask, min_n=3):
    cnt = (core.iloc[train_mask] > 0).sum(axis=0)
    return list(cnt[cnt >= min_n].index)


def design(core, a, vocab, fset="F2"):
    X = core.reindex(columns=vocab, fill_value=0.0)
    if fset in ("F1", "F2"):
        X = pd.concat([X, a[AGG1]], axis=1)
    if fset == "F2":
        X = pd.concat([X, a[AGG2]], axis=1)
    return X
