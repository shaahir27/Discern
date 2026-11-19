"""Load the fixed external KM cohort in Discern's evaluation format."""

from pathlib import Path

import numpy as np
import pandas as pd
from pyarrow import ArrowException

from discern.grammar import classify

BDQ_GENES = ["Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5"]
MECHANISM_COLUMNS = ["GENE_POSITION", "INDEL_LENGTH", "INDEL_NUCLEOTIDES", "AMINO_ACID_NUMBER"]


def load_km_cohort(ext: Path, proc: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    mutations = ext / "MUTATIONS.parquet"
    if not mutations.exists():
        raise FileNotFoundError(f"Missing {mutations} (Zenodo 15680920)")
    dst = pd.read_parquet(ext / "DST_MEASUREMENTS.parquet").reset_index()
    km = dst[(dst.SOURCE == "KMVALIDATION2025") & (dst.DRUG == "BDQ")].drop_duplicates("UNIQUEID").copy()
    km["y"] = (km.PHENOTYPE == "R").astype(int)
    genomes = pd.read_parquet(ext / "GENOMES.parquet").reset_index().set_index("UNIQUEID")
    km["lineage"] = km.UNIQUEID.map(genomes.LINEAGE)
    ids = set(km.UNIQUEID)
    train = pd.read_parquet(proc / "samples.parquet", columns=["UNIQUEID"])
    overlap = len(ids & set(train.UNIQUEID))
    assert overlap == 0, f"training overlap detected: {overlap}"

    try:
        raw = pd.read_parquet(mutations, filters=[("GENE", "in", BDQ_GENES)]).reset_index()
    except (ArrowException, NotImplementedError, TypeError, ValueError):
        raw = pd.read_parquet(mutations).reset_index()
        raw = raw[raw.GENE.isin(BDQ_GENES)]
    raw = raw[raw.UNIQUEID.isin(ids)]
    minor_col = "IS_MINOR" if "IS_MINOR" in raw.columns else "IS_MINOR_ALLELE"
    minor = raw[minor_col].astype(bool)
    suffix = np.where(minor, raw.MINOR_MUTATION.astype(str), raw.MUTATION.astype(str))
    muts = pd.DataFrame({
        "UNIQUEID": raw.UNIQUEID.values, "GENE": raw.GENE.values,
        "MUT": raw.GENE.astype(str).values + "@" + suffix,
        "FRSf": np.where(minor, raw.FRS.values, 1.0), "suffix": suffix,
    })
    for key in MECHANISM_COLUMNS:
        muts[key] = raw[key].to_numpy()
    for key, values in pd.DataFrame(muts.suffix.map(classify).tolist()).items():
        muts[key] = values
    return km, muts
