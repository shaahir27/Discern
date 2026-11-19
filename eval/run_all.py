"""Run the Discern evaluation pipeline (all experiments, or a named subset) in dependency order.

Usage (from repo root, with the discern interpreter):
    python eval/run_all.py                       # run everything in order
    python eval/run_all.py e04_lr_f2.py e05b_model_comparison.py   # run just these

Each script is run as a subprocess with the SAME interpreter that launched this file, from the
repo root, so relative paths (config.yaml, data/, results/) resolve correctly. Stops on first failure.
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

ORDER = [
    "e01_prep.py",
    "e02_who.py",
    "e03_catalogue_cv.py",
    "e04_lr_f2.py",
    "e05_ablation.py",
    "e05b_model_comparison.py",
    "e06_cross_source.py",
    "e07_loso.py",
    "e08_novel_variant.py",
    "e09_hybrid.py",
    "e10_prevalence_shift.py",
    "e11_imbalance.py",
    "e12_priority.py",
    "e13_ppv_npv.py",
]


def main():
    scripts = sys.argv[1:] or ORDER
    failed = []
    for name in scripts:
        path = ROOT / "eval" / name
        if not path.exists():
            print(f"!! missing script: {path}")
            failed.append(name)
            break
        print(f"\n=== RUN {name} ===", flush=True)
        r = subprocess.run([sys.executable, str(path)], cwd=ROOT)
        if r.returncode != 0:
            print(f"!! FAILED {name} (exit {r.returncode})", flush=True)
            failed.append(name)
            break
    print("\n=== RUN_ALL SUMMARY ===")
    print(f"ran {len(scripts) - len(failed)}/{len(scripts)} successfully")
    if failed:
        print("failed:", failed)
        sys.exit(1)
    print("ALL OK - now run: python eval/verify.py")


if __name__ == "__main__":
    main()