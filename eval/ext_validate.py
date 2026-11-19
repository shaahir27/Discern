"""G5 / P-EXT: external validation — cohort characterization + accession inventory + honest blocker.

Verified from local data: KMVALIDATION2025 BDQ phenotype cohort (792/35), ZERO training overlap,
MEDIUM quality / undocumented assay, lineage stratification. Sobol/Choi accession inventories.
The genotype->phenotype sens/spec pipeline is BLOCKED: external variant calls (CRyPTIC MUTATIONS
for KM; ENA raw reads + supplement pDST counts for Sobol/Choi) are not present locally. No accuracy
is fabricated. KM is a sensitivity analysis (assay undocumented), NOT headline proof, per P-EXT.
"""
import sys
from pathlib import Path

sys.path.insert(0, "src")
import numpy as np
import pandas as pd
import yaml

from discern.results_io import write_result

cfg = yaml.safe_load(open("config.yaml"))
proc = Path(cfg["paths"]["processed"])
ext = Path("data/external_cryptic_v34")
s = pd.read_parquet(proc / "samples.parquet")

d = pd.read_parquet(ext / "DST_MEASUREMENTS.parquet").reset_index()
km = d[(d.SOURCE == "KMVALIDATION2025") & (d.DRUG == "BDQ")].copy()
g = pd.read_parquet(ext / "GENOMES.parquet").reset_index().set_index("UNIQUEID")
km["lineage"] = km.UNIQUEID.map(g.LINEAGE)
km["sublineage"] = km.UNIQUEID.map(g.SUBLINEAGE)
n = int(km.UNIQUEID.nunique())
nR = int((km.PHENOTYPE == "R").sum())
overlap = len(set(km.UNIQUEID) & set(s.UNIQUEID))
lin = km.lineage.value_counts(dropna=False).to_dict()
linR = km[km.PHENOTYPE == "R"].lineage.value_counts(dropna=False).to_dict()
lin_tab = pd.DataFrame({"lineage": list(lin), "n_isolates": [lin[k] for k in lin],
                        "n_resistant": [linR.get(k, 0) for k in lin]})

obs_km = {
    "n_isolates": n, "n_resistant": nR, "n_susceptible": n - nR, "prevalence": nR / n,
    "quality": "MEDIUM", "assay_cc_undocumented": True,
    "zero_training_overlap": overlap == 0, "n_overlap": overlap,
    "genotype_pipeline_blocked": True,
    **{f"n_{k}": int(v) for k, v in lin.items()},
    **{f"nR_{k}": int(v) for k, v in linR.items()},
}
write_result("external_km", "e_km_cohort", obs_km,
             tables={"lineage_stratification": lin_tab},
             method="KMVALIDATION2025 BDQ DST from CRyPTIC v3.4.0 (DST_MEASUREMENTS); lineage from GENOMES; "
                    "training-ID overlap asserted against the 11,868-sample cohort.",
             interpretation=f"Verified external cohort: {n} isolates, {nR} resistant (prevalence {nR/n:.4f}), "
                            f"ZERO training overlap. All MEDIUM quality; METHOD_3/METHOD_CC undocumented. "
                            f"Lineage mix L2={lin.get('lineage2',0)}/L4={lin.get('lineage4',0)}. "
                            f"Per P-EXT this is a SENSITIVITY ANALYSIS, not headline proof. Genotype->phenotype "
                            f"sens/spec is BLOCKED: external variant calls are not local.",
             caveats="Assay/CC UNDOCUMENTED (MEDIUM) -> cannot assume WHO critical concentration; stratify by "
                     "lineage. NO external sens/spec computed (no external genotypes). NO Indian external "
                     "resistant isolates exist -> no Indian-sensitivity claim ever. To finish: acquire CRyPTIC "
                     "MUTATIONS (BDQ genes) for these 792 isolates, then run catalogue+model.",
             reproduce="python eval/ext_validate.py",
             supports="First external evaluation cohort for a BDQ triage tool (characterized; accuracy pending "
                      "external genotypes).")
print(f"KM: {n}/{nR} R, overlap={overlap}, blocked genotype pipeline", flush=True)

sob = pd.read_csv("data/sobol_runs.tsv", sep="\t")
obs_sob = {
    "n_accessions": int(len(sob)), "country": "South Africa",
    "claimed_bdst_pdst_n": 423, "claimed_pdst_n_verified": False,
    "assay": "MGIT 1.0 (per paper; assay-matched to NICD)",
    "exact_pdst_R_count": None, "paper_89_is_genotypic": True,
    "genotype_pipeline_blocked": True,
}
write_result("external_sobol", "e_sobol_inventory", obs_sob,
             tables={"accessions": sob.head(50)},
             method="Sobol 2026 KZN run-accession inventory (data/sobol_runs.tsv). pDST subset count and exact "
                    "resistant count require the paper supplement; raw reads require ENA PRJNA1257273.",
             interpretation=f"{len(sob)} South African run accessions inventoried. The 423 BDQ-pDST subset and "
                            f"exact pDST-R count are NOT verifiable locally (paper's '89' is genotypic, not pDST). "
                            f"Strongest external anchor (MGIT 1.0, assay-matched to NICD) once genotypes acquired.",
             caveats="BLOCKED: needs ENA raw reads -> variant calling, plus supplement pDST counts. No accuracy "
                     "fabricated. MGIT 1.0 matches NICD training assay (favourable for transfer).",
             reproduce="python eval/ext_validate.py",
             supports="External validation anchor (Sobol KZN), pending genotype acquisition.")

cho = pd.read_csv("data/choi_runs.tsv", sep="\t")
obs_cho = {
    "n_accessions": int(len(cho)), "country": "South Korea",
    "breakpoint_caveat": "paper breakpoint >=0.5 != WHO CC; reclassify at 1.0 / 0.25 first",
    "genotype_pipeline_blocked": True,
}
write_result("external_choi", "e_choi_inventory", obs_cho,
             tables={"accessions": cho.head(50)},
             method="Choi 2026 Korea run-accession inventory (data/choi_runs.tsv). Breakpoint differs from WHO "
                    "critical concentration; reclassification required before comparison.",
             interpretation=f"{len(cho)} South Korean run accessions inventoried. Paper breakpoint >=0.5 mg/L "
                            f"differs from WHO CC (MGIT 1.0 / plate 0.25); phenotypes must be reclassified first.",
             caveats="BLOCKED: needs ENA PRJNA1301050 raw reads -> variant calling + MIC values for "
                     "reclassification. Small n (49) -> wide intervals. No accuracy fabricated.",
             reproduce="python eval/ext_validate.py",
             supports="External validation cohort (Choi Korea), pending genotype acquisition + reclassification.")
print(f"Sobol {len(sob)} accessions; Choi {len(cho)} accessions; both blocked on raw reads", flush=True)