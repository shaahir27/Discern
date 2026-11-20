"""Native lab-export format auto-detection and HGVS-to-GARC normalization."""
from __future__ import annotations

import csv
import io
import json
import re
from pathlib import Path

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
            aa = AA3_TO_1.get(m_stop.group(1).lower(), m_stop.group(1)[0].upper())
            return f"{aa}{m_stop.group(2)}!"

        # Missense / Synonymous: e.g. Met146Thr
        m_mis = re.match(r"^([A-Za-z]{3})(\d+)([A-Za-z]{3}|\*)$", rest)
        if m_mis:
            from_aa = AA3_TO_1.get(m_mis.group(1).lower(), m_mis.group(1)[0].upper())
            pos = m_mis.group(2)
            to_aa = AA3_TO_1.get(m_mis.group(3).lower(), m_mis.group(3)[0].upper())
            return f"{from_aa}{pos}{to_aa}"

        # 1-letter format with p. prefix: e.g. p.M146T or p.E13*
        m_1letter = re.match(r"^([A-Za-z])(\d+)([A-Za-z\*!])$", rest)
        if m_1letter:
            to_char = "!" if m_1letter.group(3) in ("*", "!") else m_1letter.group(3).upper()
            return f"{m_1letter.group(1).upper()}{m_1letter.group(2)}{to_char}"

    # Strip nucleotide prefix if present
    if s.startswith("c."):
        rest = s[2:]
        # Promoter SNV: e.g. c.-11C>A or -11C>A
        m_prom = re.match(r"^-?(\d+)([A-Za-z])>([A-Za-z])$", rest)
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


def detect_and_parse_upload(payload: bytes, filename: str = "") -> dict:
    """Auto-detect format of uploaded payload (TB-Profiler JSON, AmPORE-TB TSV/CSV, Deeplex CSV, or Discern CSV)."""
    text = payload.decode("utf-8", errors="replace").strip()
    if not text:
        return {"format": "empty", "inferred_assay": "wgs", "rows": [], "sample_ids": []}

    # 1. TB-Profiler JSON detection
    if text.startswith("{") and ("dr_variants" in text or "other_variants" in text or "tbprofiler" in text.lower()):
        try:
            data = json.loads(text)
            sid = str(data.get("id") or Path(filename).stem.split(".")[0] or "sample")
            all_raw = data.get("dr_variants", []) + data.get("other_variants", [])
            rows = []
            for item in all_raw:
                raw_gene = str(item.get("gene", "")).strip()
                canonical_gene = BDQ_GENE_ALIAS.get(raw_gene.lower())
                if canonical_gene:
                    change = str(item.get("change", "")).strip()
                    garc_mut = hgvs_to_garc(change)
                    freq = float(item.get("freq", 1.0))
                    rows.append({
                        "sample_id": sid,
                        "gene": canonical_gene,
                        "mutation": garc_mut,
                        "frs": freq,
                    })
            return {
                "format": "tb_profiler_json",
                "inferred_assay": "wgs",
                "rows": rows,
                "sample_ids": [sid],
            }
        except (json.JSONDecodeError, UnicodeDecodeError, KeyError, TypeError, ValueError):
            pass

    # 2. Delimited text parsing
    dialect_delim = "\t" if "\t" in text.splitlines()[0] else ("," if "," in text.splitlines()[0] else None)
    if not dialect_delim:
        dialect_delim = ","

    reader = csv.reader(io.StringIO(text), delimiter=dialect_delim)
    try:
        header_row = next(reader)
    except StopIteration:
        return {"format": "empty", "inferred_assay": "wgs", "rows": [], "sample_ids": []}

    headers_lower = [h.strip().lower() for h in header_row]

    # Check for Deeplex Myc-TB (Gene name, AA change, Freq (%))
    is_deeplex = any("gene" in h for h in headers_lower) and any("aa" in h or "change" in h for h in headers_lower) and any("freq" in h for h in headers_lower)
    if is_deeplex and any("%" in h or "deeplex" in filename.lower() for h in headers_lower + [filename]):
        gene_idx = next(i for i, h in enumerate(headers_lower) if "gene" in h)
        change_idx = next(i for i, h in enumerate(headers_lower) if "change" in h or "aa" in h or "variant" in h)
        freq_idx = next(i for i, h in enumerate(headers_lower) if "freq" in h or "vaf" in h)
        sid_idx = next((i for i, h in enumerate(headers_lower) if "sample" in h or "id" in h), None)

        rows = []
        sids = set()
        default_sid = Path(filename).stem.split(".")[0] or "sample"
        for line_num, parts in enumerate(reader, start=2):
            if not parts or len(parts) <= max(gene_idx, change_idx, freq_idx):
                continue
            sid = parts[sid_idx].strip() if sid_idx is not None and len(parts) > sid_idx else default_sid
            sids.add(sid)
            raw_gene = parts[gene_idx].strip()
            canonical_gene = BDQ_GENE_ALIAS.get(raw_gene.lower())
            if canonical_gene:
                mut = hgvs_to_garc(parts[change_idx].strip())
                try:
                    fval = float(parts[freq_idx].strip().rstrip("%"))
                    if fval > 1.0:
                        fval = fval / 100.0
                except ValueError:
                    fval = 1.0
                rows.append({"sample_id": sid, "gene": canonical_gene, "mutation": mut, "frs": fval})

        return {
            "format": "deeplex_csv",
            "inferred_assay": "deeplex_myc_tb",
            "rows": rows,
            "sample_ids": sorted(sids),
        }

    # Check for AmPORE-TB / wf-tb-amr (sample, gene/locus, variant, vaf/af/freq)
    is_ampore = any("ampore" in filename.lower() or "wf-tb" in filename.lower() for _ in [1]) or (
        dialect_delim == "\t" and any("variant" in h for h in headers_lower) and any("vaf" in h or "af" in h or "freq" in h for h in headers_lower)
    )
    if is_ampore:
        gene_idx = next(i for i, h in enumerate(headers_lower) if "gene" in h or "locus" in h)
        change_idx = next(i for i, h in enumerate(headers_lower) if "variant" in h or "change" in h or "mutation" in h)
        freq_idx = next(i for i, h in enumerate(headers_lower) if "vaf" in h or "af" in h or "freq" in h)
        sid_idx = next((i for i, h in enumerate(headers_lower) if "sample" in h or "id" in h), None)

        rows = []
        sids = set()
        default_sid = Path(filename).stem.split(".")[0] or "sample"
        for line_num, parts in enumerate(reader, start=2):
            if not parts or len(parts) <= max(gene_idx, change_idx, freq_idx):
                continue
            sid = parts[sid_idx].strip() if sid_idx is not None and len(parts) > sid_idx else default_sid
            sids.add(sid)
            raw_gene = parts[gene_idx].strip()
            canonical_gene = BDQ_GENE_ALIAS.get(raw_gene.lower())
            if canonical_gene:
                mut = hgvs_to_garc(parts[change_idx].strip())
                try:
                    fval = float(parts[freq_idx].strip().rstrip("%"))
                    if fval > 1.0:
                        fval = fval / 100.0
                except ValueError:
                    fval = 1.0
                rows.append({"sample_id": sid, "gene": canonical_gene, "mutation": mut, "frs": fval})

        return {
            "format": "ampore_tb_tsv",
            "inferred_assay": "ampore_tb",
            "rows": rows,
            "sample_ids": sorted(sids),
        }

    # Default: Standard Discern CSV
    return {
        "format": "discern_csv",
        "inferred_assay": None,
        "rows": [],
        "sample_ids": [],
    }
