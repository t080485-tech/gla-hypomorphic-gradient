#!/usr/bin/env python3
"""
07_tenon_clinical_series.py — [CONTROLLED ACCESS] Clinical & biochemical
analysis of the referred Hôpital Tenon EXOME_TNN series (Group T1, n=8; Group T3, n=59).
Pseudonymized individual-level workbooks are NOT distributed in the public repository;
this script runs when local controlled-access files are provided.
"""
import argparse
import sys
from pathlib import Path
import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config as cfg


def main():
    cfg.ensure_dirs()
    parser = argparse.ArgumentParser()
    parser.add_argument("--t1-xlsx", type=Path, default=cfg.ROOT / "data/raw/Fabry_T1_patients.local.xlsx")
    parser.add_argument("--t3-xlsx", type=Path, default=cfg.ROOT / "data/raw/Fabry_T3_suspicion_sans_variant.local.xlsx")
    args = parser.parse_args()

    if not args.t1_xlsx.exists() or not args.t3_xlsx.exists():
        print("[CONTROLLED ACCESS] Local Tenon Excel workbooks not present — skipping step 07.")
        print("See Table 5, Table 6, and Supplementary Tables S6 & S10 for published pseudonymized summaries.")
        return

    t1 = pd.read_excel(args.t1_xlsx, sheet_name="synthese")
    t3 = pd.read_excel(args.t3_xlsx, sheet_name="synthese")

    f_t1 = int((t1["sexe_redcap"] == "Femme").sum())
    m_t1 = len(t1) - f_t1
    f_t3 = int((t3["sexe_redcap"] == "Femme").sum())
    m_t3 = len(t3) - f_t3
    odds, p = fisher_exact([[f_t1, m_t1], [f_t3, m_t3]])
    print(f"Group T1 (n={len(t1)}): {f_t1}F/{m_t1}M ({100*f_t1/len(t1):.1f}% F)")
    print(f"Group T3 (n={len(t3)}): {f_t3}F/{m_t3}M ({100*f_t3/len(t3):.1f}% F)")
    print(f"Female enrichment T1 vs T3: OR = {odds:.2f}, Fisher p = {p:.4f}")


if __name__ == "__main__":
    main()
