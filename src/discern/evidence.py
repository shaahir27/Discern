"""Per-variant evidence cards from the production evidence store (carriers, resistant, Wilson CI,
by-source). NICD is enriched for resistance — cards carry a per-source note.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

NOTE = "NICD is enriched for resistance; read per-source numbers."


class EvidenceStore:
    def __init__(self, path: Path):
        self.df = pd.read_csv(path).set_index("variant")

    def has(self, variant: str) -> bool:
        return variant in self.df.index

    def card(self, variant: str) -> dict | None:
        if variant not in self.df.index:
            return None
        r = self.df.loc[variant]
        return {
            "variant": variant,
            "n_carriers": int(r.n_carriers), "n_resistant": int(r.n_resistant),
            "pct_resistant": float(r.pct_resistant),
            "ci95_low": float(r.ci95_low), "ci95_high": float(r.ci95_high),
            "by_source": {"CRyPTIC": {"n": int(r["CRyPTIC_n"]), "r": int(r["CRyPTIC_r"])},
                          "NICD": {"n": int(r["NICD_n"]), "r": int(r["NICD_r"])}},
            "note": NOTE,
        }