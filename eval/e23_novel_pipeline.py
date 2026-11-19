"""E23: Implementation and Cross-Validation of:
1. mmpL5 lineage-filtered rare variants (targeting the 15.8% invisible resistance)
2. Bacterial Lineage indicators (Lineage 2 Beijing vs Lineage 4)
3. Continuous log2(MIC) modeling
4. External validation on KM cohort (792 samples)
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, "src")
import numpy as np
import pandas as pd
import yaml
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import average_precision_score, roc_auc_score

from discern.features import design, sample_tables, vocab_from
from discern.model import fit_nested, score, sens_at_spec

cfg = yaml.safe_load(open("config.yaml"))
raw, proc = Path(cfg["paths"]["raw"]), Path(cfg["paths"]["processed"])
ext = Path("data/external_cryptic_v34")

print("=" * 75)
print("E23: INTEGRATED NOVEL PIPELINE EVALUATION (11,868 SAMPLES)")
print("=" * 75)

# 1. Load data
t0 = time.time()
s = pd.read_parquet(proc / "samples.parquet")
m = pd.read_parquet(proc / "mutations.parquet")
y = s.y.values
ids = s.UNIQUEID.values
g = s.SUBJID.values
folds = np.load(proc / "s1_folds.npy")

# Map lineage from GENOMES
genomes = pd.read_parquet(ext / "GENOMES.parquet").reset_index().set_index("UNIQUEID")
s["lineage"] = s.UNIQUEID.map(genomes.LINEAGE).fillna("unknown")
is_l2 = (s.lineage == "lineage2").astype(float).values
is_l4 = (s.lineage == "lineage4").astype(float).values
df_lin = pd.DataFrame({"lineage_beijing_l2": is_l2, "lineage_euroamer_l4": is_l4}, index=ids)

# Parse continuous log2(MIC)
def parse_mic_to_log2(row):
    val = row["METHOD_MIC"]
    y_val = row["y"]
    if pd.isna(val):
        return 1.0 if y_val == 1 else -5.0  # Impute sensible default
    s_val = str(val).strip()
    if s_val.startswith("<="):
        num = float(s_val[2:])
    elif s_val.startswith("<"):
        num = float(s_val[1:]) / 2.0
    elif s_val.startswith(">="):
        num = float(s_val[2:])
    elif s_val.startswith(">"):
        num = float(s_val[1:]) * 2.0
    else:
        try:
            num = float(s_val)
        except:
            num = 2.0 if y_val == 1 else 0.03
    return np.log2(max(num, 0.004))

log2_mic = s.apply(parse_mic_to_log2, axis=1).values

# 2. Extract core features (F2 baseline)
core, a = sample_tables(m, ids, cfg["features"]["frs_min"], tuple(cfg["features"]["dbd_range"]))

# Prepare mmpL5 mutations for filtering
m_l5 = m[(m.GENE == "mmpL5") & (~m.IS_SYN) & (~m.IS_NULL) & (m.FRSf >= 0.1)].copy()

print(f"Loaded {len(s)} isolates. Core mutations ready. mmpL5 rows: {len(m_l5)}.")

# 3. 5-Fold Cross-Validation on S1
P_base = np.zeros(len(s))
P_l5 = np.zeros(len(s))
P_l5_lin = np.zeros(len(s))
P_mic = np.zeros(len(s))
P_combo = np.zeros(len(s))

for k in range(5):
    te = (folds == k)
    tr = ~te
    tr_ids = set(ids[tr])
    
    # Golden F2 baseline
    vocab_base = vocab_from(core, tr, cfg["features"]["vocab_min_count"])
    X_base = design(core, a, vocab_base, "F2")
    final_b, platt_b, _ = fit_nested(X_base.values[tr], y[tr], g[tr], cfg["model"]["C"], cfg["model"]["max_iter"])
    _, p_b = score(final_b, platt_b, X_base.values[te])
    P_base[te] = p_b
    
    # Path 1: Lineage-filtered mmpL5 variants
    # Find benign mmpL5 polymorphisms in training susceptible isolates
    tr_s_ids = set(ids[tr & (y == 0)])
    m_l5_tr_s = m_l5[m_l5.UNIQUEID.isin(tr_s_ids)]
    benign_counts = m_l5_tr_s.MUT.value_counts()
    # Anything appearing in >= 0.5% of training susceptible isolates is considered a benign family trait
    benign_threshold = 0.005 * len(tr_s_ids)
    benign_muts = set(benign_counts[benign_counts >= benign_threshold].index)
    
    # Filter mmpL5 rows to keep ONLY rare candidate mutations
    m_l5_rare = m_l5[~m_l5.MUT.isin(benign_muts)]
    
    # Candidate vocabulary in training: appears at least 3 times in training
    m_l5_tr = m_l5_rare[m_l5_rare.UNIQUEID.isin(tr_ids)]
    rare_counts = m_l5_tr.MUT.value_counts()
    rare_vocab = set(rare_counts[rare_counts >= 3].index)
    
    # Build sample feature for rare mmpL5
    # Feature A: specific rare mmpL5 mutations that passed the filter
    # Feature B: mmpL5_rare_any (max FRS of any rare non-benign mmpL5 mutation)
    l5_rare_by_uid = m_l5_rare.groupby("UNIQUEID").FRSf.max()
    mmpL5_rare_any = pd.Series(ids, index=ids).map(l5_rare_by_uid).fillna(0.0).values
    df_l5_feat = pd.DataFrame({"mmpL5_rare_damaging": mmpL5_rare_any}, index=ids)
    
    # Model 1: F2 + Filtered mmpL5
    X_l5 = pd.concat([X_base, df_l5_feat], axis=1)
    final_l5, platt_l5, _ = fit_nested(X_l5.values[tr], y[tr], g[tr], cfg["model"]["C"], cfg["model"]["max_iter"])
    _, p_l5 = score(final_l5, platt_l5, X_l5.values[te])
    P_l5[te] = p_l5
    
    # Model 2: F2 + Filtered mmpL5 + Lineage
    X_l5_lin = pd.concat([X_l5, df_lin], axis=1)
    final_lin, platt_lin, _ = fit_nested(X_l5_lin.values[tr], y[tr], g[tr], cfg["model"]["C"], cfg["model"]["max_iter"])
    _, p_lin = score(final_lin, platt_lin, X_l5_lin.values[te])
    P_l5_lin[te] = p_lin
    
    # Path 2: Continuous log2(MIC) Ridge Regression
    ridge = Ridge(alpha=10.0)
    ridge.fit(X_l5_lin.values[tr], log2_mic[tr])
    pred_mic = ridge.predict(X_l5_lin.values[te])
    P_mic[te] = pred_mic
    
    # Combined Ensemble (Normalized Classification + MIC Rank)
    # Convert predicted MIC to a 0-1 sigmoid proxy
    p_mic_prob = 1.0 / (1.0 + np.exp(-(pred_mic - (-1.0)) * 1.5))
    P_combo[te] = 0.5 * p_lin + 0.5 * p_mic_prob
    
    print(f"  Fold {k}: Filtered out {len(benign_muts)} benign mmpL5 family traits. Evaluated.")

# Compute metrics
def eval_res(name, p_arr):
    auc = float(roc_auc_score(y, p_arr))
    auprc = float(average_precision_score(y, p_arr))
    s985 = float(sens_at_spec(y, p_arr, 0.985))
    # Sensitivity specifically among the 137 invisible resistant patients
    core_muts = m[(m.GENE.isin(["Rv0678", "atpE", "pepQ"])) & (~m.IS_SYN) & (~m.IS_NULL) & (m.FRSf >= 0.1)]
    r_no_core = set(ids[y == 1]) - set(core_muts.UNIQUEID)
    mask_invis = np.isin(ids, list(r_no_core))
    # Fraction of invisible cases scored above median resistant score
    r_med = np.median(p_arr[y == 1])
    invis_caught = float((p_arr[mask_invis] >= r_med).mean())
    return {"name": name, "auc": auc, "auprc": auprc, "s985": s985, "invis_caught": invis_caught}

results = [
    eval_res("Golden Baseline (F2)", P_base),
    eval_res("F2 + Lineage-Filtered mmpL5", P_l5),
    eval_res("F2 + mmpL5 + Lineage (L2/L4)", P_l5_lin),
    eval_res("Continuous log2(MIC) Model", P_mic),
    eval_res("Combined Hybrid (Prob + MIC)", P_combo),
]

df_res = pd.DataFrame(results)

print("\n" + "=" * 75)
print("S1 CROSS-VALIDATION COMPARISON (ALL 11,868 PATIENTS):")
print("=" * 75)
print(f"{'Model Name':<32} | {'AUC':<8} | {'AUPRC':<8} | {'Sens@Spec.985':<13} | {'Invisible R Catch'}")
print("-" * 75)
for _, r in df_res.iterrows():
    delta = r['auc'] - df_res.loc[0, 'auc']
    d_str = f"({delta:+.4f})" if delta != 0 else "(base)"
    print(f"{r['name']:<32} | {r['auc']:.5f}  | {r['auprc']:.5f}  | {r['s985']:.5f}       | {r['invis_caught']*100:.1f}% {d_str}")
print("=" * 75)

# 4. External Validation on KM Cohort (792 samples)
print("\nTesting on External KM Cohort (792 independent hospital samples)...")
from eval.km_cohort import load_km_cohort
s_km, m_km = load_km_cohort(ext, proc)
y_km = s_km.y.values
km_ids = s_km.UNIQUEID.values

# Fit models on full dataset
vocab_all = vocab_from(core, np.ones(len(s), dtype=bool), cfg["features"]["vocab_min_count"])
X_train_base = design(core, a, vocab_all, "F2")

# Full benign mmpL5 filter
s_ids_all = set(ids[y == 0])
m_l5_s = m_l5[m_l5.UNIQUEID.isin(s_ids_all)]
b_counts = m_l5_s.MUT.value_counts()
b_muts = set(b_counts[b_counts >= 0.005 * len(s_ids_all)].index)
m_l5_rare_all = m_l5[~m_l5.MUT.isin(b_muts)]
l5_by_uid_all = m_l5_rare_all.groupby("UNIQUEID").FRSf.max()
l5_feat_all = pd.DataFrame({"mmpL5_rare_damaging": pd.Series(ids, index=ids).map(l5_by_uid_all).fillna(0.0).values}, index=ids)
X_train_combo = pd.concat([X_train_base, l5_feat_all, df_lin], axis=1)

# Fit Golden and Combo
m_base_final = LogisticRegression(C=cfg["model"]["C"], max_iter=cfg["model"]["max_iter"]).fit(X_train_base.values, y)
m_combo_final = LogisticRegression(C=cfg["model"]["C"], max_iter=cfg["model"]["max_iter"]).fit(X_train_combo.values, y)

# KM features
core_km, a_km = sample_tables(m_km, km_ids, cfg["features"]["frs_min"], tuple(cfg["features"]["dbd_range"]))
X_km_base = design(core_km, a_km, vocab_all, "F2")

m_l5_km = m_km[(m_km.GENE == "mmpL5") & (~m_km.IS_SYN) & (~m_km.IS_NULL) & (m_km.FRSf >= 0.1)]
m_l5_km_rare = m_l5_km[~m_l5_km.MUT.isin(b_muts)]
l5_by_uid_km = m_l5_km_rare.groupby("UNIQUEID").FRSf.max()
l5_feat_km = pd.DataFrame({"mmpL5_rare_damaging": pd.Series(km_ids, index=km_ids).map(l5_by_uid_km).fillna(0.0).values}, index=km_ids)

# KM lineage
is_km_l2 = (s_km.lineage == "lineage2").astype(float).values
is_km_l4 = (s_km.lineage == "lineage4").astype(float).values
df_km_lin = pd.DataFrame({"lineage_beijing_l2": is_km_l2, "lineage_euroamer_l4": is_km_l4}, index=km_ids)
X_km_combo = pd.concat([X_km_base, l5_feat_km, df_km_lin], axis=1)

p_km_base = m_base_final.predict_proba(X_km_base.values)[:, 1]
p_km_combo = m_combo_final.predict_proba(X_km_combo.values)[:, 1]

auc_km_base = float(roc_auc_score(y_km, p_km_base))
auc_km_combo = float(roc_auc_score(y_km, p_km_combo))

print("=" * 75)
print("KM EXTERNAL COHORT RESULTS (792 UNSEEN ISOLATES):")
print("=" * 75)
print(f"Golden Baseline KM AUC:   {auc_km_base:.5f}")
print(f"Novel Pipeline KM AUC:    {auc_km_combo:.5f}")
print(f"Delta on External:        {auc_km_combo - auc_km_base:+.5f}")
print("=" * 75)
print(f"Finished in {time.time() - t0:.1f}s.")
