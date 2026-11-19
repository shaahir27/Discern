"""Split definitions S1-S4 (AGENTS.md section 6)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold


def s1_folds(s, n_splits=5):
    return GroupKFold(n_splits=n_splits).split(s, groups=s.SUBJID)


def loso_folds(s):
    """S2: leave-one-site-out. Yields (site, train_idx, test_idx)."""
    idx = np.arange(len(s))
    for site in sorted(s.SITEID.unique()):
        te = idx[s.SITEID.values == site]
        tr = idx[s.SITEID.values != site]
        yield site, tr, te


def cross_source_folds(s):
    """S3: train on one SRC, test on the other, both directions."""
    idx = np.arange(len(s))
    for a, b in (("CRyPTIC", "NICD"), ("NICD", "CRyPTIC")):
        tr = idx[s.SRC.values == a]
        te = idx[s.SRC.values == b]
        yield f"{a}_to_{b}", tr, te


def variant_set_keys(m, ids):
    """S4 group key: sorted '|'.join of mutation strings per sample (non-syn/non-null)."""
    mm = m[~m.IS_SYN & ~m.IS_NULL]
    keys = mm.groupby("UNIQUEID").MUT.apply(lambda x: "|".join(sorted(x)))
    return pd.Series(ids).map(keys).fillna("").values