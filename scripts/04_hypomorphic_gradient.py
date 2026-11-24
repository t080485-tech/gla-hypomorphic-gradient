#!/usr/bin/env python3
"""
04_hypomorphic_gradient.py — Correlation between residual alpha-Gal A enzyme
activity (Lukas et al. 2020; Togawa et al. 2012) and gnomAD v4.1.2 allele
frequency, and calibration of the empirical pathogenicity frequency ceiling.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config as cfg


def act_class(a: float) -> str:
    if a < 1:
        return "Null (<1%)"
    if a < 10:
        return "Severe (1-10%)"
    if a < 25:
        return "Threshold (10-25%)"
    if a < 50:
        return "Moderate hypomorph (25-50%)"
    return "Pseudodeficiency (>50%)"


def main():
    cfg.ensure_dirs()
    df = pd.read_csv(cfg.GNOMAD_CSV, low_memory=False)
    for c in ["Allele Count", "Allele Number", "Allele Frequency", "Homozygote Count", "Hemizygote Count", "revel_max", "cadd"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    full = df.set_index("gnomAD ID")

    rows = []
    for gid, (prot, act_ih, act_glp, src) in cfg.RESIDUAL_ACTIVITY.items():
        if gid not in full.index:
            continue
        r = full.loc[gid]
        acts = [a for a in (act_ih, act_glp) if a is not None]
        act_mean = float(np.mean(acts))
        rows.append({
            "gnomAD ID": gid,
            "protein": prot,
            "act_inhouse": act_ih,
            "act_glp": act_glp,
            "act_mean": act_mean,
            "source_activity": src,
            "ClinVar": r["ClinVar Germline Classification"],
            "VEP": r["VEP Annotation"],
            "AF": r["Allele Frequency"],
            "AC": r["Allele Count"],
            "hemi": r["Hemizygote Count"],
            "hom": r["Homozygote Count"],
            "revel": r.get("revel_max"),
            "cadd": r.get("cadd"),
        })

    grad = pd.DataFrame(rows)
    grad["class"] = grad["act_mean"].apply(act_class)

    obs = grad[grad["AC"] > 0].copy()
    obs["logAF"] = np.log10(obs["AF"])
    rho, pval = stats.spearmanr(obs["act_mean"], obs["logAF"])
    print(f"Spearman residual activity vs allele frequency (n = {len(obs)} observed): rho = {rho:.3f}, p = {pval:.4e}")

    plp_obs = df[(df["ClinVar Germline Classification"].isin(cfg.PLP_LABELS)) &
                 (df["VEP Annotation"].isin(cfg.CODING_VEP)) &
                 (df["Allele Count"] > 0)]
    max_plp = plp_obs.loc[plp_obs["Allele Frequency"].idxmax()]
    print(f"Empirical pathogenicity ceiling: {max_plp['Protein Consequence']} AF = {max_plp['Allele Frequency']:.4e} "
          f"({max_plp['gnomAD ID']})")

    grad.to_csv(cfg.TABLES_DIR / "GLA_hypomorphic_gradient.csv", index=False)
    print("Saved output/tables/GLA_hypomorphic_gradient.csv.")


if __name__ == "__main__":
    main()
