#!/usr/bin/env python3
"""
05_biobank_brava_ot.py — Summary of BRaVa rare-variant association results
(~1.2M individuals across 10 biobanks) and OpenTargets GLA disease scores.
"""
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config as cfg


def main():
    cfg.ensure_dirs()
    brava = pd.read_csv(cfg.BRAVA_CSV)
    ot = pd.read_csv(cfg.OPENTARGETS_CSV)

    gene_lvl = brava[brava["Level"] == "Gene"]
    var_lvl = brava[brava["Level"] == "Variant"]
    print(f"BRaVa summary: {len(gene_lvl)} gene-level entries, {len(var_lvl)} variant-level entries")
    print(f"Min gene-level SKAT-O p-value: {gene_lvl['p_SKATO'].min():.4e} (Bonferroni = 1.39e-7)")
    print(f"Top OpenTargets association: {ot.iloc[0]['Disease']} (score = {ot.iloc[0]['Overall score']:.3f})")

    brava.to_csv(cfg.TABLES_DIR / "Table_S11_BRaVa_GLA.csv", index=False)
    ot.to_csv(cfg.TABLES_DIR / "Table_S12_OpenTargets_GLA.csv", index=False)
    print("Saved output/tables/Table_S11_BRaVa_GLA.csv and Table_S12_OpenTargets_GLA.csv.")


if __name__ == "__main__":
    main()
