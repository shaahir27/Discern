"""Unit tests for Native Lab-Export Importers & Mechanism-Stratified Reflex Guidance."""
from __future__ import annotations

import json

from discern.decision import Predictor
from discern.importers import detect_and_parse_upload, hgvs_to_garc


def test_hgvs_to_garc_conversions():
    # Missense
    assert hgvs_to_garc("p.Met146Thr") == "M146T"
    assert hgvs_to_garc("p.Arg50Gln") == "R50Q"
    assert hgvs_to_garc("p.Ala63Pro") == "A63P"
    
    # Nonsense / Stop
    assert hgvs_to_garc("p.Glu13*") == "E13!"
    assert hgvs_to_garc("p.Glu13Ter") == "E13!"
    assert hgvs_to_garc("p.Gln131*") == "Q131!"

    # Promoter
    assert hgvs_to_garc("c.-11C>A") == "c-11a"
    assert hgvs_to_garc("c.-8T>C") == "t-8c"

    # Duplications & Insertions
    assert hgvs_to_garc("c.192dupG") == "192_ins_g"
    assert hgvs_to_garc("c.192_193insG") == "192_ins_g"
    assert hgvs_to_garc("c.141_142insC") == "141_ins_c"

    # Deletions
    assert hgvs_to_garc("c.193delG") == "193_del_g"
    assert hgvs_to_garc("c.606delC") == "606_del_c"

    # GARC Idempotency
    assert hgvs_to_garc("M146T") == "M146T"
    assert hgvs_to_garc("E13!") == "E13!"
    assert hgvs_to_garc("192_ins_g") == "192_ins_g"
    assert hgvs_to_garc("c-11a") == "c-11a"


def test_tb_profiler_json_import():
    tb_profiler_sample = {
        "id": "ERR1234567",
        "dr_variants": [
            {
                "gene": "mmpR5",
                "change": "p.Met146Thr",
                "freq": 0.82,
                "type": "missense",
            },
            {
                "gene": "rpoB",  # Not a BDQ core gene; should be filtered out
                "change": "p.Ser450Leu",
                "freq": 0.99,
            }
        ],
        "other_variants": [
            {
                "gene": "mmpL5",
                "change": "c.606delC",
                "freq": 0.85,
                "type": "frameshift",
            }
        ]
    }
    raw_bytes = json.dumps(tb_profiler_sample).encode("utf-8")
    parsed = detect_and_parse_upload(raw_bytes, "ERR1234567.results.json")
    assert parsed["format"] == "tb_profiler_json"
    assert parsed["inferred_assay"] == "wgs"
    assert len(parsed["rows"]) == 2  # Rv0678 and mmpL5; rpoB excluded

    # Verify rows normalized to GARC
    v1 = parsed["rows"][0]
    assert v1["sample_id"] == "ERR1234567"
    assert v1["gene"] == "Rv0678"
    assert v1["mutation"] == "M146T"
    assert v1["frs"] == 0.82

    v2 = parsed["rows"][1]
    assert v2["sample_id"] == "ERR1234567"
    assert v2["gene"] == "mmpL5"
    assert v2["mutation"] == "606_del_c"
    assert v2["frs"] == 0.85


def test_ampore_and_deeplex_import():
    # AmPORE-TB TSV header & row
    ampore_tsv = (
        b"sample\tgene\tvariant\tvaf\n"
        b"SMPL01\tmmpR5\tp.Gln131*\t0.92\n"
        b"SMPL01\tatpE\tp.Ala63Pro\t0.88\n"
    )
    p_ampore = detect_and_parse_upload(ampore_tsv, "ampore_run.tsv")
    assert p_ampore["format"] == "ampore_tb_tsv"
    assert p_ampore["inferred_assay"] == "ampore_tb"
    assert len(p_ampore["rows"]) == 2
    assert p_ampore["rows"][0]["gene"] == "Rv0678"
    assert p_ampore["rows"][0]["mutation"] == "Q131!"

    # Deeplex Myc-TB CSV with percentages (0-100)
    deeplex_csv = (
        b"Sample ID,Gene name,AA change,Freq (%)\n"
        b"DPX_001,Rv0678,p.Met146Thr,85.5\n"
    )
    p_deeplex = detect_and_parse_upload(deeplex_csv, "deeplex_export.csv")
    assert p_deeplex["format"] == "deeplex_csv"
    assert p_deeplex["inferred_assay"] == "deeplex_myc_tb"
    assert p_deeplex["rows"][0]["frs"] == 0.855
    assert p_deeplex["rows"][0]["mutation"] == "M146T"


def test_reason_codes_and_clinical_reflex():
    pred = Predictor()

    # 1. Target resistance (atpE) -> CFZ-sparing reflex
    res_atpe = pred.predict([{"gene": "atpE", "mutation": "A63P", "frs": 0.95}])
    assert "Target-site (atpE)" in " ".join(res_atpe["clinical_reflex"])
    assert "CFZ-sparing" in " ".join(res_atpe["clinical_reflex"])
    assert any(rc["role"] == "driving_R" for rc in res_atpe["tier_a"]["reason_codes"])

    # 2. Efflux resistance (Rv0678) -> CFZ cross-resistance + MGIT reflex
    res_rv = pred.predict([{"gene": "Rv0678", "mutation": "141_ins_c", "frs": 0.95}])
    assert "Efflux-mediated (Rv0678)" in " ".join(res_rv["clinical_reflex"])
    assert "clofazimine (CFZ) cross-resistance" in " ".join(res_rv["clinical_reflex"])
    assert any(rc["role"] == "driving_R" for rc in res_rv["tier_a"]["reason_codes"])

    # 3. Late truncation Rv0678 (stop at 162: Q162!) -> surfaces Snobre 2024 & Schön 2025
    res_late = pred.predict([{"gene": "Rv0678", "mutation": "Q162!", "frs": 0.95}])
    reflex_text = " ".join(res_late["clinical_reflex"])
    assert "Late Rv0678 truncation" in reflex_text
    assert "Snobre et al., AAC 2024" in reflex_text
    assert "Schön et al., AAC 2025" in reflex_text


def test_web_review_tb_profiler_upload():
    from fastapi.testclient import TestClient

    from discern.api.app import app

    with TestClient(app) as client:
        tb_profiler_sample = {
            "id": "TB_TEST_001",
            "dr_variants": [
                {
                    "gene": "mmpR5",
                    "change": "p.Met146Thr",
                    "freq": 0.82,
                    "type": "missense",
                }
            ],
            "other_variants": []
        }
        raw_payload = json.dumps(tb_profiler_sample)

        # Post as form
        resp = client.post("/review", data={"csv": raw_payload, "assay": "wgs"}, follow_redirects=False)
        assert resp.status_code == 303
        loc = resp.headers["location"]
        assert "sample_id=TB_TEST_001" in loc

        # Follow redirect
        res_page = client.get(loc)
        assert res_page.status_code == 200
        assert "TB_TEST_001" in res_page.text
        assert "Call drivers" in res_page.text
        assert "Rv0678 Domain Track" in res_page.text
        assert "Efflux-mediated (Rv0678)" in res_page.text
