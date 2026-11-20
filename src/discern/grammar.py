"""GARC-style mutation grammar (AGENTS.md section 5)."""
from __future__ import annotations

import re

_PROM = re.compile(r"^[acgt]-\d+")
_MIS = re.compile(r"^[A-Z]\d+[A-Z]$")
_POS = re.compile(r"(\d+)")
_INDEL = re.compile(r"_ins_|_del_|indel|del_")
_MUTATION = re.compile(
    r"(?:[A-Z][0-9]+[A-Z!]|[acgtnxz]-?[0-9]+[acgtnxz]"
    r"|-?[0-9]+_(?:(?:ins|del|indel|mixed)(?:_(?:[acgtnxz]+|[0-9]+))?|minorindel)"
    r"|del_(?:0(?:\.[0-9]+)?|1(?:\.0+)?))"
)


def valid_mutation(suffix: str) -> bool:
    return len(suffix) <= 512 and _MUTATION.fullmatch(suffix) is not None


def suffix_of(mut: str) -> str:
    return mut.split("@")[-1]


def classify(suffix: str) -> dict:
    """Return parsing flags for a mutation suffix (the part after GENE@)."""
    is_syn = len(suffix) >= 1 and suffix[0] == suffix[-1]
    is_null = suffix.endswith("Z")
    ok = not is_syn and not is_null
    is_stop = ok and suffix.endswith("!")
    is_indel = ok and bool(_INDEL.search(suffix))
    is_prom = ok and bool(_PROM.match(suffix))
    is_mis = ok and bool(_MIS.match(suffix))
    pos = float(_POS.search(suffix).group(1)) if is_mis else float("nan")
    return {
        "IS_SYN": is_syn,
        "IS_NULL": is_null,
        "is_stop": is_stop,
        "is_indel": is_indel,
        "is_prom": is_prom,
        "is_mis": is_mis,
        "pos": pos,
    }
