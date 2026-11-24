#!/usr/bin/env python3
"""
03_ancestry_stratification.py — Ancestry-stratified GLA prevalence under Model B
and ancestry distribution of the five reclassified benign alleles.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config as cfg


def main():
    cfg.ensure_dirs()
    df = pd.read_csv(cfg.GNOMAD_CSV, low_memory=False)
    modB = pd.read_csv(cfg.TABLES_DIR / "ModelB_GLA.csv", low_memory=False)

    denoms = {}
    for anc in cfg.ANCESTRIES:
        an_max = pd.to_numeric(df[f"Allele Number {anc}"], errors="coerce").max()
        n_m = cfg.GLOBAL_MALE_FRAC * an_max
        n_f = (an_max - n_m) / 2.0
        denoms[anc] = {"AN_max": an_max, "N_males_est": n_m, "N_fem_est": n_f}

    rows = []
    print(f"\n{'Ancestry':28s} | {'AC':>5s} | {'sumAF':>10s} | {'Hemi':>4s} | {'Male 1/X':>10s} | {'Fem 1/X':>10s}")
    print("-" * 82)
    for anc in cfg.ANCESTRIES:
        ac_s = pd.to_numeric(modB[f"Allele Count {anc}"], errors="coerce").fillna(0)
        an_s = pd.to_numeric(modB[f"Allele Number {anc}"], errors="coerce")
        hem_s = pd.to_numeric(modB[f"Hemizygote Count {anc}"], errors="coerce").fillna(0)
        ac = float(ac_s.sum())
        sum_af = float((ac_s / an_s).sum())
        hemi = float(hem_s.sum())
        d = denoms[anc]
        male_prev = (d["N_males_est"] / hemi) if hemi > 0 else np.nan
        fem_prev = (1.0 / (2.0 * sum_af)) if sum_af > 0 else np.nan
        rows.append({
            "ancestry": anc,
            "AN_max": d["AN_max"],
            "N_males_est": d["N_males_est"],
            "AC": ac,
            "sumAF": sum_af,
            "hemi": hemi,
            "male_prev_1inx": male_prev,
            "fem_carrier_1inx": fem_prev,
        })
        m_str = f"1/{male_prev:,.0f}" if pd.notna(male_prev) else "0"
        f_str = f"1/{fem_prev:,.0f}" if pd.notna(fem_prev) else "0"
        print(f"{anc:28s} | {ac:5.0f} | {sum_af:10.6f} | {hemi:4.0f} | {m_str:>10s} | {f_str:>10s}")

    pd.DataFrame(rows).to_csv(cfg.TABLES_DIR / "GLA_ancestry_ModelB.csv", index=False)

    # Reclassified benign alleles across ancestries
    gidx = df.set_index("gnomAD ID")
    benign_ids = {
        "X-101398432-C-A": "p.Asp313Tyr",
        "X-101403828-G-A": "p.Arg118Cys",
        "X-101401752-C-T": "p.Ala143Thr",
        "X-101401803-T-C": "p.Ser126Gly",
        "X-101403984-C-G": "p.Glu66Gln",
    }
    b_rows = []
    for gid, prot in benign_ids.items():
        r = gidx.loc[gid]
        row = {"gnomAD ID": gid, "Variant": prot, "AF_global": r["Allele Frequency"], "Hemizygotes": r["Hemizygote Count"]}
        for anc in cfg.ANCESTRIES:
            ac = pd.to_numeric(r[f"Allele Count {anc}"], errors="coerce")
            an = pd.to_numeric(r[f"Allele Number {anc}"], errors="coerce")
            row[f"AF_{anc}"] = (ac / an) if (pd.notna(ac) and pd.notna(an) and an > 0) else np.nan
        b_rows.append(row)
    pd.DataFrame(b_rows).to_csv(cfg.TABLES_DIR / "GLA_benign_ancestry.csv", index=False)
    print("\nSaved output/tables/GLA_ancestry_ModelB.csv and GLA_benign_ancestry.csv.")


if __name__ == "__main__":
    main()
