"""Pydantic request schemas (AGENTS.md 12.3). Gene enum + grammar-based mutation validation
(reject free text) + FRS/prevalence/sample_id bounds.
"""
from __future__ import annotations

import re

from pydantic import BaseModel, ConfigDict, Field, ValidationInfo, field_validator

from ..decision import GENES

_SAMPLE_ID = re.compile(r"[A-Za-z0-9_-]{1,40}")
_MUTATION = re.compile(
    r"(?:[A-Z][0-9]+[A-Z!]|[acgtnxz]-?[0-9]+[acgtnxz]"
    r"|-?[0-9]+_(?:(?:ins|del|indel|mixed)(?:_(?:[acgtnxz]+|[0-9]+))?|minorindel)"
    r"|del_(?:0(?:\.[0-9]+)?|1(?:\.0+)?))"
)


def valid_mutation(suf: str) -> bool:
    return _MUTATION.fullmatch(suf) is not None


class Variant(BaseModel):
    model_config = ConfigDict(extra="forbid", hide_input_in_errors=True)

    gene: str
    mutation: str
    frs: float = Field(1.0, ge=0.0, le=1.0)

    @field_validator("gene")
    @classmethod
    def _gene(cls, v):
        if v not in GENES:
            raise ValueError(f"gene must be one of {sorted(GENES)}")
        return v

    @field_validator("mutation")
    @classmethod
    def _mutation(cls, v: str, info: ValidationInfo) -> str:
        prefix, separator, suffix = v.strip().partition("@")
        if separator and prefix != info.data.get("gene"):
            raise ValueError("mutation prefix must match the gene")
        suf = suffix if separator else prefix
        if not valid_mutation(suf):
            raise ValueError("mutation must be standard variant notation (free text rejected)")
        return suf


class PredictRequest(BaseModel):
    model_config = ConfigDict(extra="forbid", hide_input_in_errors=True)

    sample_id: str | None = None
    prevalence: float = Field(0.05, ge=0.001, le=0.5)
    variants: list[Variant] = Field(default_factory=list)

    @field_validator("sample_id")
    @classmethod
    def _sid(cls, v):
        if v is not None and not _SAMPLE_ID.fullmatch(v):
            raise ValueError("sample_id must match ^[A-Za-z0-9_-]{1,40}$")
        return v
