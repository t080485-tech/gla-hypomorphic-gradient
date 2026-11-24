#!/usr/bin/env python3
"""
01_curation_models.py — GLA prevalence estimation under three curation models
(Model A: ClinVar-strict; Model B: Expert Strategy; Model C: in silico)
using gnomAD v4.1.2 (GRCh38/hg38, ENST00000218516.4 / NM_000169.3).
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config as cfg


def quality_pass(row: pd.Series) -> bool:
    fe, fg = row["Filters - exomes"], row["Filters - genomes"]
    if fe == "PASS":
        return True
    if pd.isna(fe) or "AC0" in str(fe):
        return fg == "PASS"
    return False


def summarize_model(mod: pd.DataFrame, name: str) -> dict:
    sum_af = float(mod["Allele Frequency"].sum())
    hemi = float(mod["Hemizygote Count"].sum())
    fem_theo = 2.0 * sum_af
    male_prev = cfg.N_MALES / hemi if hemi > 0 else np.nan
    fem_prev_theo = 1.0 / fem_theo if fem_theo > 0 else np.nan
    print(f"\n{'=' * 80}")
    print(f"MODEL {name} — {len(mod)} variants | sumAF = {sum_af:.6f}")
    if hemi > 0:
        print(f"  Male hemizygotes   : {hemi:.0f} / {cfg.N_MALES:,} = 1 in {male_prev:,.0f}")
    else:
        print(f"  Male hemizygotes   : 0 / {cfg.N_MALES:,}")
    print(f"  Female carriers    : 2*sumAF = {fem_theo:.6f} = 1 in {fem_prev_theo:,.0f}")
    return {
        "model": name,
        "n_variants": len(mod),
        "sumAF": sum_af,
        "hemizygotes": hemi,
        "male_prev_1inX": male_prev,
        "female_carrier_theo_1inX": fem_prev_theo,
    }


def main():
    cfg.ensure_dirs()
    df = pd.read_csv(cfg.GNOMAD_CSV, low_memory=False)
    for c in ["Allele Count", "Allele Number", "Allele Frequency",
              "Homozygote Count", "Hemizygote Count", "cadd", "revel_max"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df["qpass"] = df.apply(quality_pass, axis=1)
    df["obs"] = df["Allele Count"] > 0
    df["evidence_status"] = df["gnomAD ID"].map(lambda g: cfg.EVIDENCE.get(g, ("", "", ""))[0])
    df["evidence_note"] = df["gnomAD ID"].map(lambda g: cfg.EVIDENCE.get(g, ("", "", ""))[1])
    df["evidence_refs"] = df["gnomAD ID"].map(lambda g: cfg.EVIDENCE.get(g, ("", "", ""))[2])

    # Model A: ClinVar P/LP strict
    modA = df[
        df["ClinVar Germline Classification"].isin(cfg.PLP_LABELS)
        & df["VEP Annotation"].isin(cfg.CODING_VEP)
        & df["obs"]
        & df["qpass"]
    ].copy()

    # Model B: Expert Strategy = A + HC pLoF (PVS1, MAF <= 0.01%) + validated late-onset
    plof = df[df["VEP Annotation"].isin(cfg.PLOF_VEP)].copy()
    plof_hc = plof[~plof["Flags"].isin(["lof_flag", "lc_lof"])]
    plof_hc_unc = plof_hc[
        (~plof_hc["ClinVar Germline Classification"].isin(cfg.PLP_LABELS))
        & plof_hc["obs"]
        & plof_hc["qpass"]
        & (plof_hc["Allele Frequency"] <= cfg.MAF_MAX)
    ]
    late_val = df[(df["evidence_status"] == "LATE_ONSET_VALIDATED") & df["obs"] & df["qpass"]]
    modB = pd.concat([modA, plof_hc_unc, late_val]).drop_duplicates(subset="gnomAD ID")

    # Model C: in silico comparator (pLoF + missense REVEL >= 0.5 or CADD >= 20, excl. Benign/LB)
    miss = df[(df["VEP Annotation"] == "missense_variant") & df["obs"] & df["qpass"]].copy()
    miss_insilico = miss[(miss["revel_max"] >= cfg.REVEL_MIN) | (miss["cadd"] >= cfg.CADD_MIN)]
    miss_insilico = miss_insilico[~miss_insilico["ClinVar Germline Classification"].isin(cfg.BENIGN_LABELS)]
    modC = pd.concat([plof[plof["obs"] & plof["qpass"]], miss_insilico]).drop_duplicates(subset="gnomAD ID")

    resA = summarize_model(modA, "A (ClinVar strict)")
    resB = summarize_model(modB, "B (Expert Strategy)")
    resC = summarize_model(modC, "C (in silico)")

    # Classic pLoF vs late-onset split within Model B
    classic = modB[modB["VEP Annotation"].isin(cfg.PLOF_VEP)]
    late = modB[~modB["VEP Annotation"].isin(cfg.PLOF_VEP)]
    print(f"\nModel B split — pLoF (classic): {len(classic)} variants, hemi={classic['Hemizygote Count'].sum():.0f}")
    print(f"Model B split — Missense/intronic (late-onset): {len(late)} variants, hemi={late['Hemizygote Count'].sum():.0f}")

    # Export
    modA.to_csv(cfg.TABLES_DIR / "ModelA_GLA.csv", index=False)
    modB.to_csv(cfg.TABLES_DIR / "ModelB_GLA.csv", index=False)
    modC.to_csv(cfg.TABLES_DIR / "ModelC_GLA.csv", index=False)
    pd.DataFrame([resA, resB, resC]).to_csv(cfg.TABLES_DIR / "GLA_prevalence_summary.csv", index=False)

    ev = df[df["evidence_status"] != ""][
        ["gnomAD ID", "Position", "Protein Consequence", "VEP Annotation",
         "ClinVar Germline Classification", "Allele Count", "Allele Frequency",
         "Hemizygote Count", "evidence_status", "evidence_note", "evidence_refs"]
    ]
    ev.to_csv(cfg.TABLES_DIR / "GLA_evidence_table.csv", index=False)
    print("\nSaved tables to output/tables/.")


if __name__ == "__main__":
    main()
