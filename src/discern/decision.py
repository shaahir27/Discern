"""Tier A catalogue decision with the optional D011 veto and Tier B probability and ranking.
Offline; loads only local production artefacts. No request bodies are logged here.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import piezo
import yaml

from .catalogue import _aggregate
from .evidence import EvidenceStore
from .features import design, sample_tables
from .grammar import classify
from .mechanism import apply_veto, veto_flags
from .model import score
from .panels import audit_sample_panel
from .prevalence import adjust, archived_predictive_values

GENES = {"Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5"}
CORE = {"Rv0678", "atpE", "pepQ"}
LABEL = {"R": "Resistant", "S": "Susceptible", "U": "Send for lab test"}
BANNER = ("Research use only - not a diagnostic test. Final call comes from the catalogue; "
          "the model only ranks unresolved samples and gives a probability.")


class CatalogueEvaluationError(RuntimeError):
    """The catalogue could not evaluate an eligible variant."""


class Predictor:
    def __init__(
        self, models_dir: Path = Path("models/lr_model"), results_dir: Path = Path("results"),
        *, config_path: Path = Path("config.yaml"),
    ):
        models_dir, results_dir = Path(models_dir), Path(results_dir)
        self.art = joblib.load(models_dir / "bdq_lr.joblib")
        self.rc = piezo.ResistanceCatalogue(str(models_dir / "catalogue_bdq.csv"))
        self.ev = EvidenceStore(models_dir / "evidence_store.csv")
        self.meta = json.loads((models_dir / "metadata.json").read_text(encoding="utf-8"))
        self.catver = json.loads((models_dir / "catalogue_version.json").read_text(encoding="utf-8"))
        self.pi_tr = float(self.art["pi_tr"])
        self.th = self.art["thresholds"]
        self.vocab = self.art["vocab"]
        e03 = json.loads((results_dir / "approaches/catalogue/e03_s1_catalogue_cv/results.json").read_text())["observed"]
        self.cat_metrics = e03
        self.cat_sens, self.cat_spec = float(e03["sens_firm"]), float(e03["spec_firm"])
        config_path = Path(config_path)
        config = {}
        if config_path.is_file():
            config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
        self.veto_high_conf: set[str] | None = None
        self.lof_effects: dict[str, dict] = {}
        effects_path = models_dir.parent / "lr_mechanism/lof_effects.csv"
        if effects_path.is_file():
            eff_df = pd.read_csv(effects_path)
            for _, erow in eff_df.iterrows():
                self.lof_effects[erow["MUT"]] = erow.to_dict()
        if config.get("enhanced", {}).get("epistatic_veto", True):
            if effects_path.is_file():
                effects = pd.read_csv(effects_path, usecols=["MUT", "high_conf"])
                self.veto_high_conf = set(effects.loc[effects.high_conf, "MUT"])
            else:
                self.veto_high_conf = set()
            e16_path = results_dir / "models/lr_mechanism/e16_lof_ablation/results.json"
            if e16_path.is_file():
                e16 = json.loads(e16_path.read_text(encoding="utf-8"))["observed"]
                self.cat_sens, self.cat_spec = e16["cat_s1_sens_firm_on"], e16["cat_s1_spec_firm_on"]
                metrics_csv = e16_path.with_name("s1_veto_metrics.csv")
                if metrics_csv.is_file():
                    self.cat_metrics = pd.read_csv(metrics_csv).set_index("state").loc["ON"].to_dict()
        replay = pd.read_parquet(results_dir / "oof_predictions.parquet")
        self.u_p = replay.loc[replay.cat_call == "U", "p"].values

    def _mutation_rows(self, variants, sample_id):
        rows = []
        for v in variants:
            gene, mut = v["gene"], v["mutation"]
            suf = mut.split("@")[-1]
            frs = float(v.get("frs", 1.0))
            flags = classify(suf)
            rows.append({"UNIQUEID": sample_id, "GENE": gene, "MUT": f"{gene}@{suf}",
                         "FRSf": frs, "suffix": suf, **flags})
        if not rows:
            return pd.DataFrame(columns=["UNIQUEID", "GENE", "MUT", "FRSf", "suffix",
                                         "IS_SYN", "IS_NULL", "is_stop", "is_indel", "is_prom", "is_mis", "pos"])
        return pd.DataFrame(rows)

    def _call_variant(self, mut: str, frs: float) -> str:
        if frs < self.art["frs_min"]:
            return "S"
        try:
            r = self.rc.predict(mut)
        except Exception:  # noqa: BLE001 — contain failures from the catalogue library
            raise CatalogueEvaluationError("Catalogue evaluation failed.") from None
        call = r.get("BDQ") if isinstance(r, dict) else r
        if not isinstance(call, str) or call not in LABEL:
            raise CatalogueEvaluationError("Catalogue returned an unsupported call.")
        return call

    def predict(self, variants, prevalence: float = 0.05, sample_id: str | None = None, assay: str = "wgs") -> dict:
        sid = sample_id or "input"
        mdf = self._mutation_rows(variants, sid)
        core, a = sample_tables(mdf, [sid], self.art["frs_min"], tuple(self.art["dbd"]))
        X = design(core, a, self.vocab, "F2").values

        variant_calls = []
        catalogue_failed = False
        for _, r in mdf.iterrows():
            c = "S"
            if r.GENE in CORE and not r.IS_SYN and not r.IS_NULL and r.FRSf >= self.art["frs_min"]:
                try:
                    c = self._call_variant(r.MUT, r.FRSf)
                except CatalogueEvaluationError:
                    c = "U"
                    catalogue_failed = True
            variant_calls.append({"variant": r.MUT, "gene": r.GENE, "frs": float(r.FRSf), "call": c})
        tier_a = {"variant_calls": variant_calls}
        calls = np.array([vc["call"] for vc in variant_calls], dtype=object)
        veto_applied = False
        if not catalogue_failed and self.veto_high_conf:
            veto_targets = (calls == "R") & (
                mdf.GENE.eq("Rv0678") & mdf.MUT.isin(self.veto_high_conf)
            ).to_numpy()
            if veto_targets.any():
                retained = mdf[(~mdf.IS_SYN) & (~mdf.IS_NULL) & (mdf.FRSf >= self.art["frs_min"])]
                flags = veto_flags(retained, [sid], self.veto_high_conf)
                if flags.veto_fire.iloc[0]:
                    calls = apply_veto(calls, veto_targets)
                    tier_a["epistatic_veto_applied"] = True
                    veto_applied = True
        call = "U" if catalogue_failed else _aggregate(pd.Series(calls))

        reason_codes = []
        for vc, v_call in zip(variant_calls, calls, strict=True):
            if vc["call"] != v_call:
                role = "vetoed_R_to_U"
            elif call == "R" and v_call == "R":
                role = "driving_R"
            elif call == "U" and v_call == "U":
                role = "driving_U"
            else:
                role = "benign_or_below_frs"
            reason_codes.append({"variant": vc["variant"], "call": vc["call"], "role": role})
        tier_a["reason_codes"] = reason_codes

        clinical_reflex = []
        has_atpe_r_or_u = any(
            vc["gene"] == "atpE" and vc["call"] in ("R", "U") for vc in variant_calls
        )
        if has_atpe_r_or_u:
            clinical_reflex.append(
                "Target-site (atpE): typically confers high-level BDQ resistance that is CFZ-sparing (efflux-independent) and unaffected by mmpL5/mmpS5 loss-of-function."
            )

        has_rv_r_or_u = any(
            vc["gene"] == "Rv0678" and vc["call"] in ("R", "U") for vc in variant_calls
        )
        if has_rv_r_or_u:
            clinical_reflex.append(
                "Efflux-mediated (Rv0678): typically confers low-to-moderate/borderline BDQ MIC elevation with high likelihood of clofazimine (CFZ) cross-resistance. Reflex to MGIT960 BDQ MIC and monitor intermediate-growth tubes."
            )

        has_late_trunc = False
        for _, r in mdf.iterrows():
            if r.GENE == "Rv0678" and not r.IS_SYN and not r.IS_NULL and r.FRSf >= self.art["frs_min"]:
                lof_info = self.lof_effects.get(r.MUT)
                if lof_info and (lof_info.get("first_stop_aa", 0) > 160 or lof_info.get("late_preserved")):
                    has_late_trunc = True
                    break
                if r.is_stop:
                    m_stop_pos = re.search(r"(\d+)", r.suffix)
                    if m_stop_pos and int(m_stop_pos.group(1)) > 160:
                        has_late_trunc = True
                        break

        if has_late_trunc:
            clinical_reflex.append(
                "Late Rv0678 truncation (first stop > aa 160): protein dimerization domain may be preserved (Snobre et al., AAC 2024), though susceptibility evidence is actively debated due to inter-lab MIC shifts (Schön et al., AAC 2025). Confirm with phenotypic DST."
            )

        if veto_applied:
            clinical_reflex.append(
                "Efflux pump inactivation (mmpL5/mmpS5 LoF): associated with BDQ hypersusceptibility; negates Rv0678-only efflux resistance when genetically linked (WHO Catalogue v3)."
            )

        raw, p = score(self.art["final"], self.art["platt"], X)
        raw, p = float(raw[0]), float(p[0])
        p_adj = float(adjust(p, self.pi_tr, prevalence))
        rank_pct = float((self.u_p < p).mean() * 100) if len(self.u_p) else 0.0

        veto_fired = bool(tier_a.get("epistatic_veto_applied", False))
        warnings = self._warnings(mdf, a, call, p, veto_fired=veto_fired)
        if catalogue_failed:
            warnings.append("CATALOGUE_EVALUATION_FAILED")

        panel_audit = audit_sample_panel(variants, assay)
        for pw in panel_audit.get("warnings", []):
            if pw not in warnings:
                warnings.append(pw)

        evidence = []
        for mut in mdf.MUT.unique():
            card = self.ev.card(mut)
            if card is None:
                card = {"variant": mut, "n_carriers": 0, "n_resistant": 0, "pct_resistant": None,
                        "ci95_low": None, "ci95_high": None, "by_source": {}, "note": "Novel variant - not in evidence store."}
            else:
                card = dict(card)
            card["lof_annotation"] = self.lof_effects.get(mut)
            matching_rows = mdf[mdf.MUT == mut]
            if not matching_rows.empty:
                r0 = matching_rows.iloc[0]
                pos_val = r0["pos"] if not np.isnan(r0["pos"]) else None
                if pos_val is None:
                    m_digits = re.search(r"(\d+)", r0["suffix"])
                    if m_digits:
                        pos_val = int(m_digits.group(1))
                card["pos"] = pos_val
            evidence.append(card)
        ppv, npv = archived_predictive_values(self.cat_metrics, prevalence)
        return {
            "sample_id": sid, "call": call, "call_label": LABEL[call],
            "tier_a": tier_a,
            "panel_audit": panel_audit,
            "model": {"score_raw": raw, "probability_train_prevalence": p,
                      "probability_at_prevalence": p_adj, "queue_rank_percentile": rank_pct},
            "predictive_values": {"prevalence": prevalence, "ppv_if_R": ppv, "npv_if_S": npv,
                                  "catalogue_sens": self.cat_sens, "catalogue_spec": self.cat_spec},
            "evidence": evidence, "warnings": warnings,
            "clinical_reflex": clinical_reflex,
            "versions": {"package": self.meta.get("package_versions", {}), "catalogue": self.catver.get("version"),
                         "model_git_commit": self.meta.get("git_commit"),
                         "data_manifest_sha256": self.meta.get("data_manifest_sha256")},
            "banner": BANNER,
        }

    def _warnings(self, mdf, a, call, p, veto_fired: bool = False) -> list[str]:
        w = []
        contrib = mdf[(mdf.GENE.isin(CORE)) & (~mdf.IS_SYN) & (~mdf.IS_NULL) & (mdf.FRSf >= self.art["frs_min"])]
        if len(contrib) and (contrib.FRSf < 0.9).all():
            w.append("MINOR_ALLELE_ONLY")
        if any(self.art["frs_min"] <= float(r.FRSf) < 0.75 for _, r in contrib.iterrows()):
            w.append("HETERORESISTANCE_SUBCLONAL")
        if any(not self.ev.has(mut) for mut in contrib.MUT):
            w.append("NOVEL_VARIANT")
        if any(self.ev.has(mut) and self.ev.df.loc[mut, "n_carriers"] < 10 for mut in contrib.MUT):
            w.append("LOW_EVIDENCE")
        if mdf.IS_NULL.any():
            w.append("NULL_CALL")
        has_rv = (mdf.GENE == "Rv0678").any()
        if float(a["mmpL5_LoF"].iloc[0]) > 0 and has_rv:
            w.append("PUMP_DAMAGE_MODIFIER")
            if veto_fired or not (mdf.GENE == "atpE").any():
                w.append("EPISTATIC_PUMP_INACTIVATION")
        if call == "S" and len(contrib) == 0:
            w.append("NO_VARIANT_INVISIBLE_RESISTANCE_POSSIBLE")
        if (call == "S" and p >= self.th["hi"]) or (call == "R" and p <= self.th["lo"]):
            w.append("MODEL_DISAGREES")
        return w
