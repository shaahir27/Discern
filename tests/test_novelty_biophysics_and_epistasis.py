"""Unit and regression tests for novel biophysics, epistatic veto, and subclonal heteroresistance."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from discern.decision import Predictor
from discern.grammar import classify
from discern.mechanism import apply_veto, veto_flags
from discern.structure import (
    grantham_distance,
    charge_delta,
    parse_pdb,
    sample_biophysical_features,
)


def test_grantham_and_charge_calculations():
    # Conservative substitution: I to V is 29
    assert grantham_distance("I", "V") == 29
    # Radical substitution: I to S is 142
    assert grantham_distance("I", "S") == 142
    # Self distance is 0
    assert grantham_distance("A", "A") == 0
    # Symmetric lookup
    assert grantham_distance("S", "I") == 142

    # Charge difference
    assert charge_delta("R", "D") == 2.0  # Basic (+1) to Acidic (-1)
    assert charge_delta("A", "A") == 0.0


def test_sample_biophysical_features_recognition_helix():
    # Dummy structure map
    struct = {
        64: {"min_dna_dist": 3.2, "min_dimer_dist": 15.0, "burial": 12},
        120: {"min_dna_dist": 22.0, "min_dimer_dist": 4.1, "burial": 18},
    }
    
    # Mutation in recognition helix alpha4 (codon 64: R64P)
    df_rec = pd.DataFrame([
        {"GENE": "Rv0678", "MUT": "Rv0678@R64P", "suffix": "R64P", "pos": 64, "is_mis": True}
    ])
    feat_rec = sample_biophysical_features(df_rec, struct)
    assert feat_rec["rec_helix"] == 1
    assert feat_rec["helix_breaker"] == 1
    assert feat_rec["max_grantham"] == grantham_distance("R", "P")  # 103

    # Mutation outside recognition helix (codon 120: L120V)
    df_out = pd.DataFrame([
        {"GENE": "Rv0678", "MUT": "Rv0678@L120V", "suffix": "L120V", "pos": 120, "is_mis": True}
    ])
    feat_out = sample_biophysical_features(df_out, struct)
    assert feat_out["rec_helix"] == 0
    assert feat_out["helix_breaker"] == 0
    assert feat_out["max_grantham"] == grantham_distance("L", "V")  # 32
    assert feat_out["min_dimer_dist"] == 4.1


def test_epistatic_veto_firing_and_target_protection():
    # Case 1: Rv0678 frameshift + mmpL5 frameshift (NO atpE) -> Epistatic veto MUST fire
    retained_veto = pd.DataFrame([
        {"UNIQUEID": "iso_veto", "GENE": "Rv0678", "MUT": "Rv0678@192_ins_g",
         "is_indel": True, "is_stop": False, "FRSf": 1.0},
        {"UNIQUEID": "iso_veto", "GENE": "mmpL5", "MUT": "mmpL5@141_ins_c",
         "is_indel": True, "is_stop": False, "FRSf": 1.0},
    ])
    flags_veto = veto_flags(retained_veto, ["iso_veto"], {"Rv0678@192_ins_g"})
    assert flags_veto.veto_fire.iloc[0] == True
    
    # Calls before veto: ["R"]
    calls = np.array(["R"])
    vetoed_calls = apply_veto(calls, flags_veto.veto_fire.values)
    assert vetoed_calls[0] == "U"  # Downgraded from R to U!

    # Case 2: atpE target resistance co-occurring -> Epistatic veto MUST NEVER touch atpE!
    retained_atpe = pd.DataFrame([
        {"UNIQUEID": "iso_atpe", "GENE": "Rv0678", "MUT": "Rv0678@192_ins_g",
         "is_indel": True, "is_stop": False, "FRSf": 1.0},
        {"UNIQUEID": "iso_atpe", "GENE": "mmpL5", "MUT": "mmpL5@141_ins_c",
         "is_indel": True, "is_stop": False, "FRSf": 1.0},
        {"UNIQUEID": "iso_atpe", "GENE": "atpE", "MUT": "atpE@A63P",
         "is_indel": False, "is_stop": False, "FRSf": 1.0},
    ])
    flags_atpe = veto_flags(retained_atpe, ["iso_atpe"], {"Rv0678@192_ins_g"})
    assert flags_atpe.veto_fire.iloc[0] == False  # Protected because atpE is present

    # Case 3: Minor unlinked mmpL5 subclone (e.g. 10%) with fixed dominant Rv0678 (95%)
    # Under WHO Catalogue v3 / TB-Profiler v6.6.5, this must NOT veto!
    retained_unlinked = pd.DataFrame([
        {"UNIQUEID": "iso_unlinked", "GENE": "Rv0678", "MUT": "Rv0678@192_ins_g",
         "is_indel": True, "is_stop": False, "FRSf": 0.95},
        {"UNIQUEID": "iso_unlinked", "GENE": "mmpL5", "MUT": "mmpL5@141_ins_c",
         "is_indel": True, "is_stop": False, "FRSf": 0.10},
    ])
    flags_unlinked = veto_flags(retained_unlinked, ["iso_unlinked"], {"Rv0678@192_ins_g"})
    assert flags_unlinked.frs_linked.iloc[0] == False
    assert flags_unlinked.veto_fire.iloc[0] == False  # Dominant Rv0678 clone remains active!

    # Case 4: Linked subclone (e.g. Rv0678 90% and mmpL5 85%) -> veto fires!
    retained_linked = pd.DataFrame([
        {"UNIQUEID": "iso_linked", "GENE": "Rv0678", "MUT": "Rv0678@192_ins_g",
         "is_indel": True, "is_stop": False, "FRSf": 0.90},
        {"UNIQUEID": "iso_linked", "GENE": "mmpL5", "MUT": "mmpL5@141_ins_c",
         "is_indel": True, "is_stop": False, "FRSf": 0.85},
    ])
    flags_linked = veto_flags(retained_linked, ["iso_linked"], {"Rv0678@192_ins_g"})
    assert flags_linked.frs_linked.iloc[0] == True
    assert flags_linked.veto_fire.iloc[0] == True


def test_subclonal_heteroresistance_warning(tmp_path):
    predictor = Predictor()
    
    # Subclonal variant (FRS = 0.35) -> should trigger HETERORESISTANCE_SUBCLONAL
    res_subclonal = predictor.predict([
        {"gene": "Rv0678", "mutation": "141_ins_c", "frs": 0.35}
    ])
    assert "HETERORESISTANCE_SUBCLONAL" in res_subclonal["warnings"]
    assert "MINOR_ALLELE_ONLY" in res_subclonal["warnings"]

    # Clonal variant (FRS = 0.95) -> should NOT trigger HETERORESISTANCE_SUBCLONAL
    res_clonal = predictor.predict([
        {"gene": "Rv0678", "mutation": "141_ins_c", "frs": 0.95}
    ])
    assert "HETERORESISTANCE_SUBCLONAL" not in res_clonal["warnings"]
    assert "MINOR_ALLELE_ONLY" not in res_clonal["warnings"]


def test_epistatic_veto_end_to_end():
    predictor = Predictor()
    
    # Rv0678 + mmpL5 LoF -> Epistatic veto fires, call downgraded from R to U, warning added
    res = predictor.predict([
        {"gene": "Rv0678", "mutation": "141_ins_c", "frs": 1.0},
        {"gene": "mmpL5", "mutation": "141_ins_c", "frs": 1.0},
    ])
    assert res["tier_a"].get("epistatic_veto_applied") is True
    assert res["call"] == "U"
    assert "EPISTATIC_PUMP_INACTIVATION" in res["warnings"]
    assert "PUMP_DAMAGE_MODIFIER" in res["warnings"]
