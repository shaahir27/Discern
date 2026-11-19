"""E22: Rigorous Full-Cohort (11,868 samples) Evaluation of Biophysical Features on S1 and KM External Cohort."""
import sys
import time
from pathlib import Path

sys.path.insert(0, "src")
import numpy as np
import pandas as pd
import yaml
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score

from discern.features import design, sample_tables, vocab_from
from discern.model import fit_nested, score, sens_at_spec, thresholds
from discern.structure import parse_pdb, sample_biophysical_features
from discern.metrics import ternary_metrics

cfg = yaml.safe_load(open("config.yaml"))
raw, proc = Path(cfg["paths"]["raw"]), Path(cfg["paths"]["processed"])

print("=" * 70)
print("E22: FULL-COHORT (11,868 ISOLATES) BIOPHYSICS EVALUATION")
print("=" * 70)

# Load data
t0 = time.time()
s = pd.read_parquet(proc / "samples.parquet")
m = pd.read_parquet(proc / "mutations.parquet")
struct = parse_pdb(raw / "data/Rv0678_WT.pdb")
y = s.y.values
ids = s.UNIQUEID.values
g = s.SUBJID.values
folds = np.load(proc / "s1_folds.npy")

print(f"Loaded {len(s)} patient samples ({int(y.sum())} resistant, {int((y==0).sum())} susceptible).")

# 1. Extract biophysical features for ALL samples
print("Extracting biophysical features for all 11,868 samples...")
rv_muts = m[(m.GENE == "Rv0678") & (~m.IS_SYN) & (~m.IS_NULL) & (m.FRSf >= 0.1)]
by_uid = {uid: df for uid, df in rv_muts.groupby("UNIQUEID")}
empty_df = pd.DataFrame()

bio_list = []
for uid in ids:
    u_df = by_uid.get(uid, empty_df)
    bio_list.append(sample_biophysical_features(u_df, struct, tuple(cfg["features"]["dbd_range"])))

BIO_COLS = ["max_grantham", "rec_helix", "helix_breaker", "min_dna_dist", "max_burial", "max_charge"]
df_bio = pd.DataFrame(bio_list, index=ids)[BIO_COLS]
print(f"Biophysical feature matrix shape: {df_bio.shape}")

# 2. Extract baseline F2 features
core, a = sample_tables(m, ids, cfg["features"]["frs_min"], tuple(cfg["features"]["dbd_range"]))

# 3. Run 5-fold Patient-grouped CV (S1) for Baseline vs Enhanced
print("\nRunning 5-fold Patient-Grouped CV (S1)...")
P_base = np.zeros(len(s))
P_enh = np.zeros(len(s))

for k in range(5):
    te = (folds == k)
    tr = ~te
    vocab = vocab_from(core, tr, cfg["features"]["vocab_min_count"])
    
    # Baseline F2
    X_base = design(core, a, vocab, "F2")
    final_b, platt_b, poof_b = fit_nested(X_base.values[tr], y[tr], g[tr], cfg["model"]["C"], cfg["model"]["max_iter"])
    _, p_b = score(final_b, platt_b, X_base.values[te])
    P_base[te] = p_b
    
    # Enhanced F2 + Biophysics
    X_enh = pd.concat([X_base, df_bio], axis=1)
    final_e, platt_e, poof_e = fit_nested(X_enh.values[tr], y[tr], g[tr], cfg["model"]["C"], cfg["model"]["max_iter"])
    _, p_e = score(final_e, platt_e, X_enh.values[te])
    P_enh[te] = p_e
    
    print(f"  Fold {k}: Base features={X_base.shape[1]}, Enhanced features={X_enh.shape[1]}")

# Metrics computation
auc_b = float(roc_auc_score(y, P_base))
auc_e = float(roc_auc_score(y, P_enh))
auprc_b = float(average_precision_score(y, P_base))
auprc_e = float(average_precision_score(y, P_enh))
s985_b = float(sens_at_spec(y, P_base, 0.985))
s985_e = float(sens_at_spec(y, P_enh, 0.985))

print("\n" + "=" * 70)
print("S1 CROSS-VALIDATION RESULTS (11,868 SAMPLES):")
print("=" * 70)
print(f"Metric                    | Baseline F2 | F2 + Biophysics | Delta")
print(f"--------------------------|-------------|-----------------|-----------")
print(f"AUC                       | {auc_b:.5f}     | {auc_e:.5f}         | {auc_e - auc_b:+.5f}")
print(f"AUPRC                     | {auprc_b:.5f}     | {auprc_e:.5f}         | {auprc_e - auprc_b:+.5f}")
print(f"Sensitivity at Spec 0.985 | {s985_b:.5f}     | {s985_e:.5f}         | {s985_e - s985_b:+.5f}")
print("=" * 70)

# 4. Now evaluate on external KM cohort (792 samples)
ext_path = Path("data/external_cryptic_v34")
if ext_path.exists():
    try:
        from eval.km_cohort import load_km_cohort
        print("\nEvaluating generalization on External KM Cohort (792 samples)...")
        s_km, m_km = load_km_cohort(ext_path, proc)
        y_km = s_km.y.values
        km_ids = s_km.UNIQUEID.values
        
        # Train both models on ALL 11,868 training samples
        vocab_all = vocab_from(core, np.ones(len(s), dtype=bool), cfg["features"]["vocab_min_count"])
        X_train_base = design(core, a, vocab_all, "F2")
        X_train_enh = pd.concat([X_train_base, df_bio], axis=1)
        
        # Model fits
        m_base = LogisticRegression(C=cfg["model"]["C"], max_iter=cfg["model"]["max_iter"])
        m_base.fit(X_train_base.values, y)
        
        m_enh = LogisticRegression(C=cfg["model"]["C"], max_iter=cfg["model"]["max_iter"])
        m_enh.fit(X_train_enh.values, y)
        
        # KM features
        core_km, a_km = sample_tables(m_km, km_ids, cfg["features"]["frs_min"], tuple(cfg["features"]["dbd_range"]))
        X_km_base = design(core_km, a_km, vocab_all, "F2")
        
        rv_km = m_km[(m_km.GENE == "Rv0678") & (~m_km.IS_SYN) & (~m_km.IS_NULL) & (m_km.FRSf >= 0.1)]
        by_uid_km = {uid: df for uid, df in rv_km.groupby("UNIQUEID")}
        bio_km = []
        for uid in km_ids:
            u_df = by_uid_km.get(uid, empty_df)
            bio_km.append(sample_biophysical_features(u_df, struct, tuple(cfg["features"]["dbd_range"])))
        df_bio_km = pd.DataFrame(bio_km, index=km_ids)[BIO_COLS]
        X_km_enh = pd.concat([X_km_base, df_bio_km], axis=1)
        
        p_km_base = m_base.predict_proba(X_km_base.values)[:, 1]
        p_km_enh = m_enh.predict_proba(X_km_enh.values)[:, 1]
        
        auc_km_b = float(roc_auc_score(y_km, p_km_base))
        auc_km_e = float(roc_auc_score(y_km, p_km_enh))
        
        print("=" * 70)
        print("KM EXTERNAL COHORT RESULTS (792 UNSEEN SAMPLES):")
        print("=" * 70)
        print(f"Baseline KM AUC:     {auc_km_b:.5f}")
        print(f"Biophysics KM AUC:   {auc_km_e:.5f}")
        print(f"Delta on External:   {auc_km_e - auc_km_b:+.5f}")
        print("=" * 70)
    except Exception as exc:
        print(f"Could not load KM cohort: {exc}")

print(f"\nExecution finished in {time.time() - t0:.1f}s.")
