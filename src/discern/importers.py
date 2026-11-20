"""Native lab-export format auto-detection and HGVS-to-GARC normalization."""
from __future__ import annotations

import csv
import io
import json
import math
import re
from pathlib import Path

from .research import UploadError

AA3_TO_1 = {
    "ala": "A", "arg": "R", "asn": "N", "asp": "D", "cys": "C",
    "gln": "Q", "glu": "E", "gly": "G", "his": "H", "ile": "I",
    "leu": "L", "lys": "K", "met": "M", "phe": "F", "pro": "P",
    "ser": "S", "thr": "T", "trp": "W", "tyr": "Y", "val": "V",
    "ter": "!", "*": "!", "stop": "!",
}

BDQ_GENE_ALIAS = {
    "mmpr5": "Rv0678",
    "rv0678": "Rv0678",
    "atpe": "atpE",
    "pepq": "pepQ",
    "mmpl5": "mmpL5",
    "mmps5": "mmpS5",
}


def hgvs_to_garc(change: str) -> str:
    """Normalize HGVS syntax (p.Met146Thr, p.Glu13*, c.-11C>A, c.192dupG) into exact GARC syntax."""
    s = str(change).strip()
    if not s:
        return s

    # Strip protein prefix if present
    if s.startswith("p."):
        rest = s[2:]
        # Nonsense: e.g. Glu13* or Glu13Ter
        m_stop = re.match(r"^([A-Za-z]{3}|\*)(\d+)(?:\*|Ter|ter|\!)$", rest)
        if m_stop:
            aa = AA3_TO_1.get(m_stop.group(1).lower())
            return f"{aa}{m_stop.group(2)}!" if aa else s

        # Missense / Synonymous: e.g. Met146Thr
        m_mis = re.match(r"^([A-Za-z]{3})(\d+)([A-Za-z]{3}|\*)$", rest)
        if m_mis:
            from_aa = AA3_TO_1.get(m_mis.group(1).lower())
            pos = m_mis.group(2)
            to_aa = AA3_TO_1.get(m_mis.group(3).lower())
            return f"{from_aa}{pos}{to_aa}" if from_aa and to_aa else s

        # 1-letter format with p. prefix: e.g. p.M146T or p.E13*
        m_1letter = re.match(r"^([A-Za-z])(\d+)([A-Za-z\*!])$", rest)
        if m_1letter:
            to_char = "!" if m_1letter.group(3) in ("*", "!") else m_1letter.group(3).upper()
            return f"{m_1letter.group(1).upper()}{m_1letter.group(2)}{to_char}"

    # Strip nucleotide prefix if present
    if s.startswith("c."):
        rest = s[2:]
        # Promoter SNV: e.g. c.-11C>A or -11C>A
        m_prom = re.match(r"^-(\d+)([ACGTacgt])>([ACGTacgt])$", rest)
        if m_prom:
            ref = m_prom.group(2).lower()
            pos = m_prom.group(1)
            alt = m_prom.group(3).lower()
            return f"{ref}-{pos}{alt}"

        # Duplication: e.g. c.192dupG or c.192dup
        m_dup = re.match(r"^(\d+)(?:_\d+)?dup([A-Za-z]*)$", rest)
        if m_dup:
            nt = m_dup.group(2).lower() or "n"
            return f"{m_dup.group(1)}_ins_{nt}"

        # Insertion: e.g. c.192_193insG
        m_ins = re.match(r"^(\d+)_\d+ins([A-Za-z]+)$", rest)
        if m_ins:
            return f"{m_ins.group(1)}_ins_{m_ins.group(2).lower()}"

        # Deletion: e.g. c.193delG or c.193del or c.193_195del
        m_del = re.match(r"^(\d+)(?:_\d+)?del([A-Za-z]*)$", rest)
        if m_del:
            nt = m_del.group(2).lower() or "n"
            return f"{m_del.group(1)}_del_{nt}"

    # Already valid GARC format or unhandled string
    return s


def _fraction(value: object, line: int, *, percent: bool = False) -> float:
    text = str(value).strip()
    number = text.removesuffix("%")
    try:
        fraction = float(number) / (100 if percent or text.endswith("%") else 1)
    except (ValueError, OverflowError):
        fraction = math.nan
    if not 0 <= fraction <= 1:
        raise UploadError("invalid_fraction", "Enter a fraction or a percentage in the export's units.",
                          line=line, column="frs")
    return fraction


def detect_and_parse_upload(payload: bytes, filename: str = "") -> dict:
    """Identify supported lab exports and normalize their reported variant rows."""
    try:
        text = payload.decode("utf-8-sig").strip()
    except UnicodeError:
        raise UploadError("invalid_encoding", "Save the export as UTF-8 and try again.") from None
    default_id = re.sub(r"[^A-Za-z0-9_.-]", "_", Path(filename).stem.split(".")[0])[:64] or "sample"
    if not text:
        return {"format": "empty", "inferred_assay": "wgs", "rows": [], "sample_ids": []}
    if text.startswith("{"):
        try:
            data = json.loads(text)
        except (ValueError, RecursionError):
            raise UploadError("invalid_export", "Use a valid TB-Profiler JSON export.") from None
        groups = [data.get(name, []) for name in ("dr_variants", "other_variants")]
        if not any(name in data for name in ("dr_variants", "other_variants")) or any(
            not isinstance(group, list) for group in groups
        ):
            raise UploadError("invalid_export", "Use a TB-Profiler export with variant lists.")
        sid = str(data.get("id") or default_id).strip()
        rows = []
        for line, item in enumerate((item for group in groups for item in group), 1):
            if not isinstance(item, dict) or not isinstance(item.get("gene"), str):
                raise UploadError("invalid_export", "A variant record needs a gene name.", line=line, column="gene")
            gene = BDQ_GENE_ALIAS.get(item["gene"].strip().lower())
            if gene:
                rows.append({"sample_id": sid, "gene": gene, "mutation": hgvs_to_garc(item.get("change", "")),
                             "frs": _fraction(item.get("freq"), line)})
        return {"format": "tb_profiler_json", "inferred_assay": "wgs", "rows": rows, "sample_ids": [sid]}

    delimiter = "\t" if "\t" in text.splitlines()[0] else ","
    reader = csv.reader(io.StringIO(text, newline=""), delimiter=delimiter, strict=True)
    try:
        headers = [h.strip().lower() for h in next(reader, [])]
        if headers[:4] == ["sample_id", "gene", "mutation", "frs"] or headers[:1] == ["sample_id"] and len(headers) <= 2:
            return {"format": "discern_csv", "inferred_assay": None, "rows": [], "sample_ids": []}

        def column(*names: str) -> int | None:
            return next((i for i, header in enumerate(headers) if header in names), None)

        gene_idx = column("gene", "gene name", "locus")
        change_idx = column("variant", "change", "mutation", "aa change", "amino_acid_change")
        freq_idx = column("freq", "vaf", "af", "freq (%)", "variant frequency", "variant frequency (%)")
        sid_idx = column("sample", "sample_id", "sample id", "id")
        if None not in (gene_idx, change_idx, freq_idx):
            deeplex = (
                "%" in headers[freq_idx] or headers[freq_idx] == "variant frequency"
                or "gene name" in headers or "deeplex" in filename.lower()
            )
            ampore = delimiter == "\t" or sid_idx is not None or any(
                token in filename.lower() for token in ("ampore", "wf-tb")
            )
            if deeplex or ampore:
                rows = []
                sample_ids = set()
                for parts in reader:
                    if not parts or not any(part.strip() for part in parts):
                        continue
                    line = reader.line_num
                    if len(parts) != len(headers):
                        raise UploadError("invalid_export", "Use the same columns on all rows.", line=line)
                    sid = parts[sid_idx].strip() if sid_idx is not None else default_id
                    if not sid:
                        raise UploadError("invalid_export", "Fill in the sample identifier.", line=line, column="sample_id")
                    sample_ids.add(sid)
                    gene = BDQ_GENE_ALIAS.get(parts[gene_idx].strip().lower())
                    if gene:
                        rows.append({"sample_id": sid, "gene": gene,
                                     "mutation": hgvs_to_garc(parts[change_idx]),
                                     "frs": _fraction(parts[freq_idx], line, percent=deeplex)})
                return {
                    "format": "deeplex_csv" if deeplex else "ampore_tb_tsv",
                    "inferred_assay": "deeplex_myc_tb" if deeplex else "ampore_tb",
                    "rows": rows, "sample_ids": sorted(sample_ids),
                }
    except csv.Error:
        raise UploadError("malformed_csv", "Check the export quoting and separators.", line=reader.line_num or 1) from None
    return {"format": "discern_csv", "inferred_assay": None, "rows": [], "sample_ids": []}
