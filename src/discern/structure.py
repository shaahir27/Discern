"""Structure and biophysical features from Rv0678_WT.pdb (AGENTS.md section 11.1 & novel extension).
Plain-text PDB parser and physicochemical mutation matrices (Grantham, BLOSUM62).

Chains C/D (DNA) are a modeled complex; chains A/B form the functional homodimer.
Features include:
- Distance to DNA operator (modeled)
- Distance to homodimer interface (Chain B)
- Grantham physicochemical distance (composition, polarity, volume)
- Recognition helix alpha4 (residues 62-68, major groove DNA contact)
- Helix-breaker substitutions (proline / glycine disruptions)
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

# Grantham matrix (Grantham R. Science 1974; 185:862-4)
GRANTHAM: dict[tuple[str, str], int] = {
    ('A', 'A'): 0, ('A', 'C'): 195, ('A', 'D'): 126, ('A', 'E'): 107, ('A', 'F'): 113,
    ('A', 'G'): 60, ('A', 'H'): 86, ('A', 'I'): 94, ('A', 'K'): 106, ('A', 'L'): 96,
    ('A', 'M'): 84, ('A', 'N'): 111, ('A', 'P'): 27, ('A', 'Q'): 91, ('A', 'R'): 112,
    ('A', 'S'): 99, ('A', 'T'): 58, ('A', 'V'): 64, ('A', 'W'): 148, ('A', 'Y'): 112,
    ('C', 'C'): 0, ('C', 'D'): 154, ('C', 'E'): 170, ('C', 'F'): 205, ('C', 'G'): 159,
    ('C', 'H'): 174, ('C', 'I'): 198, ('C', 'K'): 202, ('C', 'L'): 198, ('C', 'M'): 196,
    ('C', 'N'): 139, ('C', 'P'): 169, ('C', 'Q'): 154, ('C', 'R'): 180, ('C', 'S'): 112,
    ('C', 'T'): 149, ('C', 'V'): 192, ('C', 'W'): 215, ('C', 'Y'): 194,
    ('D', 'D'): 0, ('D', 'E'): 45, ('D', 'F'): 177, ('D', 'G'): 81, ('D', 'H'): 81,
    ('D', 'I'): 168, ('D', 'K'): 101, ('D', 'L'): 172, ('D', 'M'): 160, ('D', 'N'): 23,
    ('D', 'P'): 108, ('D', 'Q'): 61, ('D', 'R'): 96, ('D', 'S'): 65, ('D', 'T'): 85,
    ('D', 'V'): 152, ('D', 'W'): 181, ('D', 'Y'): 160,
    ('E', 'E'): 0, ('E', 'F'): 140, ('E', 'G'): 98, ('E', 'H'): 40, ('E', 'I'): 134,
    ('E', 'K'): 56, ('E', 'L'): 138, ('E', 'M'): 121, ('E', 'N'): 42, ('E', 'P'): 93,
    ('E', 'Q'): 29, ('E', 'R'): 54, ('E', 'S'): 80, ('E', 'T'): 65, ('E', 'V'): 121,
    ('E', 'W'): 152, ('E', 'Y'): 122,
    ('F', 'F'): 0, ('F', 'G'): 153, ('F', 'H'): 100, ('F', 'I'): 21, ('F', 'K'): 102,
    ('F', 'L'): 22, ('F', 'M'): 28, ('F', 'N'): 158, ('F', 'P'): 132, ('F', 'Q'): 116,
    ('F', 'R'): 97, ('F', 'S'): 155, ('F', 'T'): 130, ('F', 'V'): 50, ('F', 'W'): 40, ('F', 'Y'): 32,
    ('G', 'G'): 0, ('G', 'H'): 87, ('G', 'I'): 135, ('G', 'K'): 127, ('G', 'L'): 138,
    ('G', 'M'): 127, ('G', 'N'): 80, ('G', 'P'): 42, ('G', 'Q'): 87, ('G', 'R'): 125,
    ('G', 'S'): 56, ('G', 'T'): 59, ('G', 'V'): 109, ('G', 'W'): 176, ('G', 'Y'): 147,
    ('H', 'H'): 0, ('H', 'I'): 94, ('H', 'K'): 32, ('H', 'L'): 99, ('H', 'M'): 87,
    ('H', 'N'): 68, ('H', 'P'): 77, ('H', 'Q'): 24, ('H', 'R'): 29, ('H', 'S'): 89,
    ('H', 'T'): 77, ('H', 'V'): 84, ('H', 'W'): 115, ('H', 'Y'): 83,
    ('I', 'I'): 0, ('I', 'K'): 102, ('I', 'L'): 5, ('I', 'M'): 10, ('I', 'N'): 149,
    ('I', 'P'): 95, ('I', 'Q'): 109, ('I', 'R'): 97, ('I', 'S'): 142, ('I', 'T'): 89,
    ('I', 'V'): 29, ('I', 'W'): 61, ('I', 'Y'): 33,
    ('K', 'K'): 0, ('K', 'L'): 107, ('K', 'M'): 95, ('K', 'N'): 94, ('K', 'P'): 103,
    ('K', 'Q'): 53, ('K', 'R'): 26, ('K', 'S'): 121, ('K', 'T'): 101, ('K', 'V'): 97,
    ('K', 'W'): 110, ('K', 'Y'): 85,
    ('L', 'L'): 0, ('L', 'M'): 15, ('L', 'N'): 153, ('L', 'P'): 98, ('L', 'Q'): 113,
    ('L', 'R'): 102, ('L', 'S'): 145, ('L', 'T'): 92, ('L', 'V'): 32, ('L', 'W'): 61, ('L', 'Y'): 36,
    ('M', 'M'): 0, ('M', 'N'): 142, ('M', 'P'): 87, ('M', 'Q'): 101, ('M', 'R'): 91,
    ('M', 'S'): 135, ('M', 'T'): 81, ('M', 'V'): 21, ('M', 'W'): 67, ('M', 'Y'): 36,
    ('N', 'N'): 0, ('N', 'P'): 91, ('N', 'Q'): 46, ('N', 'R'): 86, ('N', 'S'): 46,
    ('N', 'T'): 65, ('N', 'V'): 133, ('N', 'W'): 174, ('N', 'Y'): 143,
    ('P', 'P'): 0, ('P', 'Q'): 76, ('P', 'R'): 103, ('P', 'S'): 74, ('P', 'T'): 38,
    ('P', 'V'): 68, ('P', 'W'): 147, ('P', 'Y'): 110,
    ('Q', 'Q'): 0, ('Q', 'R'): 43, ('Q', 'S'): 68, ('Q', 'T'): 42, ('Q', 'V'): 96,
    ('Q', 'W'): 130, ('Q', 'Y'): 99,
    ('R', 'R'): 0, ('R', 'S'): 110, ('R', 'T'): 91, ('R', 'V'): 96, ('R', 'W'): 101, ('R', 'Y'): 77,
    ('S', 'S'): 0, ('S', 'T'): 17, ('S', 'V'): 124, ('S', 'W'): 177, ('S', 'Y'): 144,
    ('T', 'T'): 0, ('T', 'V'): 69, ('T', 'W'): 128, ('T', 'Y'): 92,
    ('V', 'V'): 0, ('V', 'W'): 88, ('V', 'Y'): 55,
    ('W', 'W'): 0, ('W', 'Y'): 37,
    ('Y', 'Y'): 0,
}

CHARGE_MAP: dict[str, float] = {'R': 1.0, 'K': 1.0, 'H': 0.5, 'D': -1.0, 'E': -1.0}
REC_HELIX_RANGE: tuple[int, int] = (62, 68)  # Recognition helix alpha4 (NSIRRRTA)


def grantham_distance(wt_aa: str, mut_aa: str) -> int:
    """Return Grantham chemical distance (0-215). Defaults to 100 if unknown."""
    return GRANTHAM.get((wt_aa, mut_aa), GRANTHAM.get((mut_aa, wt_aa), 100))


def charge_delta(wt_aa: str, mut_aa: str) -> float:
    """Absolute charge difference between amino acids."""
    return abs(CHARGE_MAP.get(mut_aa, 0.0) - CHARGE_MAP.get(wt_aa, 0.0))


def parse_pdb(path: Path) -> dict[int, dict[str, float | int]]:
    """Parse Rv0678 PDB structure, computing distances to DNA and Chain B dimer interface."""
    protein_atoms = {}  # (chain, resseq) -> list of (x,y,z)
    ca = {}  # (chain, resseq) -> (x,y,z)
    dna_atoms = []
    atoms_b = []
    
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.startswith("ATOM"):
            continue
        chain = line[21]
        try:
            resseq = int(line[22:26])
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        name = line[12:16].strip()
        if chain in ("C", "D"):
            dna_atoms.append((x, y, z))
        elif chain == "A":
            protein_atoms.setdefault((chain, resseq), []).append((x, y, z))
            if name == "CA":
                ca[(chain, resseq)] = (x, y, z)
        elif chain == "B":
            atoms_b.append((x, y, z))
            if name == "CA":
                ca[(chain, resseq)] = (x, y, z)

    dna = np.array(dna_atoms) if dna_atoms else np.empty((0, 3))
    chain_b = np.array(atoms_b) if atoms_b else np.empty((0, 3))
    ca_arr = np.array(list(ca.values())) if ca else np.empty((0, 3))
    struct = {}

    for (chain, resseq), atoms in protein_atoms.items():
        if chain != "A":
            continue
        a = np.array(atoms)
        min_dna = float(np.min(np.linalg.norm(a[:, None, :] - dna[None, :, :], axis=2))) if len(dna) else 30.0
        min_dimer = float(np.min(np.linalg.norm(a[:, None, :] - chain_b[None, :, :], axis=2))) if len(chain_b) else 30.0
        burial = 0
        if (chain, resseq) in ca and len(ca_arr):
            c = np.array(ca[(chain, resseq)])
            burial = int((np.linalg.norm(ca_arr - c, axis=1) <= 10.0).sum())
        struct[resseq] = {
            "min_dna_dist": min_dna,
            "min_dimer_dist": min_dimer,
            "burial": burial,
        }
    return struct


def sample_features(positions, struct, dbd=(34, 99)):
    """Extract sample features matching AGENTS.md section 11.1 golden interface."""
    pos = [int(p) for p in positions if not np.isnan(p)]
    in_struct = [p for p in pos if p in struct]
    outside = int(any(p not in struct for p in pos))
    min_dist = min((struct[p]["min_dna_dist"] for p in in_struct), default=30.0)
    max_burial = max((struct[p]["burial"] for p in in_struct), default=0)
    dbd_flag = int(any(dbd[0] <= p <= dbd[1] for p in pos))
    return {
        "min_dna_dist": float(min_dist),
        "max_burial": float(max_burial),
        "dbd_flag": dbd_flag,
        "n_pos": len(pos),
        "outside_structure": outside,
    }


def sample_biophysical_features(
    mutations_df, struct, dbd=(34, 99), rec_helix_range=(62, 68),
) -> dict[str, float | int]:
    """Enhanced biophysical feature extraction incorporating Grantham distance and secondary structure."""
    pos = [int(p) for p in mutations_df.pos if not np.isnan(p)] if "pos" in mutations_df else []
    in_struct = [p for p in pos if p in struct]
    outside = int(any(p not in struct for p in pos))
    min_dist = min((struct[p]["min_dna_dist"] for p in in_struct), default=30.0)
    min_dimer = min((struct[p].get("min_dimer_dist", 30.0) for p in in_struct), default=30.0)
    max_burial = max((struct[p]["burial"] for p in in_struct), default=0)
    dbd_flag = int(any(dbd[0] <= p <= dbd[1] for p in pos))
    rec_helix = int(any(rec_helix_range[0] <= p <= rec_helix_range[1] for p in pos))

    grantham_dists = []
    charge_diffs = []
    helix_breakers = []

    for _, r in mutations_df.iterrows():
        suf = str(r.get("suffix", "")).split("@")[-1]
        if len(suf) >= 3 and suf[0].isalpha() and suf[-1].isalpha():
            wt, mut = suf[0], suf[-1]
            p = r.get("pos", np.nan)
            grantham_dists.append(grantham_distance(wt, mut))
            charge_diffs.append(charge_delta(wt, mut))
            is_breaker = int(mut == 'P' or (mut == 'G' and pd.notna(p) and dbd[0] <= p <= dbd[1]))
            helix_breakers.append(is_breaker)

    max_grantham = max(grantham_dists, default=0)
    max_charge = max(charge_diffs, default=0.0)
    has_helix_breaker = max(helix_breakers, default=0)

    return {
        "n_pos": len(pos),
        "dbd_flag": dbd_flag,
        "outside_structure": outside,
        "min_dna_dist": float(min_dist),
        "min_dimer_dist": float(min_dimer),
        "max_burial": float(max_burial),
        "rec_helix": rec_helix,
        "max_grantham": float(max_grantham),
        "max_charge": float(max_charge),
        "helix_breaker": has_helix_breaker,
    }
