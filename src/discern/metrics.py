"""Ternary R/S/U metrics shared across experiments."""
from __future__ import annotations

import numpy as np


def ternary_metrics(y, c) -> dict:
    y = np.asarray(y)
    c = np.asarray(c)
    firm = c != "U"
    R = c == "R"
    S = c == "S"
    nR = max((y == 1).sum(), 1)
    nS = max((y == 0).sum(), 1)
    nRfirm = max(((y == 1) & firm).sum(), 1)
    nSfirm = max(((y == 0) & firm).sum(), 1)
    return {
        "firm": float(firm.mean()),
        "sens_firm": float((R & (y == 1)).sum() / nRfirm),
        "spec_firm": float((S & (y == 0)).sum() / nSfirm),
        "sens_U_as_S": float((R & (y == 1)).sum() / nR),
        "spec_U_as_S": float(((~R) & (y == 0)).sum() / nS),
        "sens_U_as_R": float(((R | (c == "U")) & (y == 1)).sum() / nR),
        "spec_U_as_R": float((S & (y == 0)).sum() / nS),
        "n_U": int((c == "U").sum()),
        "n_U_resistant": int(((c == "U") & (y == 1)).sum()),
    }