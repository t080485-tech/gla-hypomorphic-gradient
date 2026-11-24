#!/usr/bin/env python3
"""
02_female_carriers.py — Empirical and theoretical female carrier prevalence
in gnomAD v4.1.2 (X-linked; N_females derived from AN_max = N_males + 2*N_females).
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config as cfg


def load_model(name: str) -> pd.DataFrame:
    d = pd.read_csv(cfg.TABLES_DIR / f"Model{name}_GLA.csv", low_memory=False)
    for c in ["Allele Count", "Allele Number", "Allele Frequency", "Homozygote Count", "Hemizygote Count"]:
        d[c] = pd.to_numeric(d[c], errors="coerce").fillna(0)
    return d


def main():
    cfg.ensure_dirs()
    print(f"Denominators: N_males = {cfg.N_MALES:,} | AN_max = {cfg.AN_MAX:,} -> N_females = {cfg.N_FEMALES:,.0f}")

    rows = []
    for name in ["A", "B", "C"]:
        m = load_model(name)
        ac = float(m["Allele Count"].sum())
        hemi = float(m["Hemizygote Count"].sum())
        hom = float(m["Homozygote Count"].sum())
        fem_carriers_obs = (ac - hemi) - hom
        sum_af = float(m["Allele Frequency"].sum())
        fem_emp = fem_carriers_obs / cfg.N_FEMALES
        fem_theo = 2.0 * sum_af
        print(f"Model {name}: {fem_carriers_obs:.0f} observed female carriers (1 in {1/fem_emp:,.0f}) | "
              f"2*sumAF = 1 in {1/fem_theo:,.0f}")
        rows.append({
            "model": name,
            "n_variants": len(m),
            "AC": ac,
            "hemizygotes_male": hemi,
            "homozygotes_female": hom,
            "female_carriers_obs": fem_carriers_obs,
            "female_prev_emp_1inX": 1.0 / fem_emp if fem_emp > 0 else np.nan,
            "female_prev_theo_1inX": 1.0 / fem_theo if fem_theo > 0 else np.nan,
            "male_prev_1inX": cfg.N_MALES / hemi if hemi > 0 else np.nan,
        })

    # Classic pLoF vs late-onset split in females (Model B)
    mB = load_model("B")
    for label, sub in [
        ("pLoF (classic)", mB[mB["VEP Annotation"].isin(cfg.PLOF_VEP)]),
        ("Missense/intronic (late-onset)", mB[~mB["VEP Annotation"].isin(cfg.PLOF_VEP)]),
    ]:
        fc = (sub["Allele Count"].sum() - sub["Hemizygote Count"].sum()) - sub["Homozygote Count"].sum()
        print(f"  Model B female {label}: {len(sub)} variants, {fc:.0f} carriers (1 in {cfg.N_FEMALES/fc:,.0f})")

    pd.DataFrame(rows).to_csv(cfg.TABLES_DIR / "GLA_female_prevalence.csv", index=False)
    print("Saved output/tables/GLA_female_prevalence.csv.")


if __name__ == "__main__":
    main()
