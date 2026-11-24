#!/usr/bin/env python3
"""
run_all.py — End-to-end execution of the public + aggregated GLA analysis pipeline.
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STEPS = [
    "scripts/01_curation_models.py",
    "scripts/02_female_carriers.py",
    "scripts/03_ancestry_stratification.py",
    "scripts/04_hypomorphic_gradient.py",
    "scripts/05_biobank_brava_ot.py",
    "scripts/06_biomnis_validation.py",
    "scripts/07_tenon_clinical_series.py",
    "scripts/08_figures.py",
]


def main():
    for step in STEPS:
        print(f"\n>>> Running {step}")
        res = subprocess.run([sys.executable, str(ROOT / step)], check=True)
        if res.returncode != 0:
            sys.exit(res.returncode)
    print("\nAll pipeline steps completed successfully.")


if __name__ == "__main__":
    main()
