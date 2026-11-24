#!/usr/bin/env python3
"""
06_biomnis_validation.py — Clinical exome cohort validation (Eurofins Biomnis,
Illumina n=36,184 + MGI n=1,933 = 38,117 consecutive patients).

Modes:
  --mode aggregated (default): reproduces all enrichment, Poisson, and gradient
      statistics from the public, anonymized variant-level aggregated tables in
      data/aggregated/.
  --mode raw: parses controlled-access Biomnis Illumina + MGI VCF-derived TSV
      extractions (never committed to GitHub), performs prescriber-based index
      deduplication, checks zero overlap across platforms, and exports strictly
      anonymized variant-level counts.
"""
import argparse
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, poisson, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config as cfg


def gla_ann(row: pd.Series) -> pd.Series:
    genes = str(row["ANN[*].GENE"]).split(",")
    effects = str(row["ANN[*].EFFECT"]).split(",")
    feats = str(row["ANN[*].FEATUREID"]).split(",")
    hgvsc = str(row["ANN[*].HGVS_C"]).split(",")
    hgvsp = str(row["ANN[*].HGVS_P"]).split(",")
    for i, g in enumerate(genes):
        if g.strip() == "GLA" and (i < len(feats) and "NM_000169" in feats[i]):
            return pd.Series({
                "gla_effect": effects[i].strip() if i < len(effects) else "",
                "gla_hgvsc": hgvsc[i].strip() if i < len(hgvsc) else "",
                "gla_hgvsp": hgvsp[i].strip() if i < len(hgvsp) else "",
            })
    return pd.Series({"gla_effect": "", "gla_hgvsc": "", "gla_hgvsp": ""})


def extract_carrier_records(df: pd.DataFrame, source_name: str) -> pd.DataFrame:
    records = []
    for _, r in df.iterrows():
        key = (int(r["POS"]), r["REF_x"], r["ALT"])
        for st in ["NEPHRO", "CARDIO", "AUTRE", "NA"]:
            for gt in ["0/1", "1/1"]:
                cols = [c for c in df.columns if c == f"{st}_G_{gt}" or c.startswith(f"{st}_G_{gt}_")]
                for c in cols:
                    val = r[c]
                    if pd.notna(val) and str(val).strip() not in ("", ".", "nan"):
                        for chunk in str(val).split(";"):
                            chunk = chunk.strip()
                            if not chunk:
                                continue
                            m = re.match(r"([^:(]+):\(([^)]+)\)", chunk)
                            if m:
                                presc = m.group(1).strip()
                                for d in m.group(2).split(","):
                                    records.append({
                                        "POS": key[0], "REF_x": key[1], "ALT": key[2],
                                        "source": source_name, "stratum": st, "gt": gt,
                                        "presc": presc, "dossier": d.strip(),
                                    })
    return pd.DataFrame(records)


def run_raw(ill_path: Path, mgi_path: Path):
    ill = pd.read_csv(ill_path, sep="\t", low_memory=False)
    mgi = pd.read_csv(mgi_path, sep="\t", low_memory=False)
    modA = pd.read_csv(cfg.TABLES_DIR / "ModelA_GLA.csv", low_memory=False)
    modB = pd.read_csv(cfg.TABLES_DIR / "ModelB_GLA.csv", low_memory=False)
    modC = pd.read_csv(cfg.TABLES_DIR / "ModelC_GLA.csv", low_memory=False)
    grad = pd.read_csv(cfg.TABLES_DIR / "GLA_hypomorphic_gradient.csv")

    rec_ill = extract_carrier_records(ill, "Illumina")
    rec_mgi = extract_carrier_records(mgi, "MGI")
    assert len(set(rec_ill["dossier"]) & set(rec_mgi["dossier"])) == 0, "Overlap detected!"

    rec_all = pd.concat([rec_ill, rec_mgi], ignore_index=True)
    prio = {"NEPHRO": 0, "CARDIO": 1, "AUTRE": 2, "NA": 3}
    rec_all["prio"] = rec_all["stratum"].map(prio)
    rec_dedup = (
        rec_all.sort_values(["POS", "REF_x", "ALT", "presc", "prio", "dossier"])
        .drop_duplicates(subset=["POS", "REF_x", "ALT", "presc"], keep="first")
    )

    def agg(rec_df, suffix):
        g = rec_df.groupby(["POS", "REF_x", "ALT", "stratum"]).size().unstack(fill_value=0)
        for st in ["NEPHRO", "CARDIO", "AUTRE", "NA"]:
            if st not in g.columns:
                g[st] = 0
        g = g[["NEPHRO", "CARDIO", "AUTRE", "NA"]].copy()
        g["TOT"] = g.sum(axis=1)
        g.columns = [f"{c}_{suffix}" for c in g.columns]
        return g.reset_index()

    for df in (ill, mgi):
        ann = df.apply(gla_ann, axis=1)
        for c in ann.columns:
            df[c] = ann[c]

    meta_cols = ["POS", "REF_x", "ALT", "ID_x", "gla_effect", "gla_hgvsc", "gla_hgvsp",
                 "clinvar_CLNSIG", "gnomADex_AC", "gnomADex_AF", "gnomADex_revel_max",
                 "gnomADex_spliceai_ds_max", "dbNSFP_CADD_phred", "dbNSFP_AlphaMissense_pred"]
    meta = pd.concat([ill[meta_cols], mgi[meta_cols]], ignore_index=True).drop_duplicates(subset=["POS", "REF_x", "ALT"])

    master = (
        meta.merge(agg(rec_all, "raw"), on=["POS", "REF_x", "ALT"], how="left")
        .merge(agg(rec_dedup, "dedup"), on=["POS", "REF_x", "ALT"], how="left")
        .merge(agg(rec_ill, "ill"), on=["POS", "REF_x", "ALT"], how="left")
        .merge(agg(rec_mgi, "mgi"), on=["POS", "REF_x", "ALT"], how="left")
    )
    for c in master.columns:
        if any(c.endswith(s) for s in ["_raw", "_dedup", "_ill", "_mgi"]):
            master[c] = master[c].fillna(0).astype(int)

    is_coding = master["gla_effect"].str.contains("missense|stop_gained|frameshift|splice_acceptor|splice_donor|inframe", na=False)
    is_ivs4 = (master["POS"] == 101399747) & (master["ALT"] == "T")
    master_gla = master[(master["gla_hgvsc"] != "") & (is_coding | is_ivs4) & (master["TOT_raw"] > 0)].copy()

    keysA = set(zip(modA["Position"].astype(int), modA["Reference"], modA["Alternate"]))
    keysB = set(zip(modB["Position"].astype(int), modB["Reference"], modB["Alternate"]))
    keysC = set(zip(modC["Position"].astype(int), modC["Reference"], modC["Alternate"]))
    plp_set = {"Pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic"}
    blb_set = {"Benign", "Likely_benign", "Benign/Likely_benign", "Benign/Likely_benign|other"}

    master_gla["in_ModelA_gnomAD"] = [(int(r.POS), r.REF_x, r.ALT) in keysA for _, r in master_gla.iterrows()]
    master_gla["in_ModelB_gnomAD"] = [(int(r.POS), r.REF_x, r.ALT) in keysB for _, r in master_gla.iterrows()]
    master_gla["in_ModelC_gnomAD"] = [(int(r.POS), r.REF_x, r.ALT) in keysC for _, r in master_gla.iterrows()]
    master_gla["in_ExpertStrategy_Biomnis"] = master_gla["in_ModelB_gnomAD"] | master_gla["clinvar_CLNSIG"].isin(plp_set)
    master_gla["in_ModelC_Biomnis"] = master_gla["in_ModelC_gnomAD"] | (
        (~master_gla["clinvar_CLNSIG"].isin(blb_set)) & (~is_ivs4.loc[master_gla.index]) &
        ((pd.to_numeric(master_gla["gnomADex_revel_max"], errors="coerce") >= cfg.REVEL_MIN) |
         (pd.to_numeric(master_gla["dbNSFP_CADD_phred"], errors="coerce") >= cfg.CADD_MIN) |
         master_gla["gla_effect"].str.contains("stop_gained|frameshift|splice_acceptor|splice_donor|inframe", na=False))
    )

    AA3 = {'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G',
           'His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S',
           'Thr':'T','Trp':'W','Tyr':'Y','Val':'V','Ter':'*'}

    def three_to_one(hgvsp):
        if not isinstance(hgvsp, str) or not hgvsp.startswith("p."):
            return None
        m = re.match(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", hgvsp[2:].strip("()"))
        return f"p.{AA3[m.group(1)]}{m.group(2)}{AA3[m.group(3)]}" if (m and m.group(1) in AA3 and m.group(3) in AA3) else None

    grad["protein_short"] = grad["protein"].str.strip()
    master_gla["protein_short"] = master_gla["gla_hgvsp"].apply(three_to_one)
    mg = master_gla.merge(grad[["protein_short", "act_mean", "class"]], on="protein_short", how="left")
    es = master_gla[master_gla["in_ExpertStrategy_Biomnis"]].copy()

    es.to_csv(cfg.BIOMNIS_ES_CSV, index=False)
    mg.to_csv(cfg.BIOMNIS_ALL_CSV, index=False)
    return es, mg


def summarize_aggregated(es: pd.DataFrame, all_v: pd.DataFrame):
    N = cfg.BIOMNIS_DENOMS
    rate_50_50 = 0.5 * (57.0 / cfg.N_MALES) + 0.5 * (105.0 / cfg.N_FEMALES)
    expected = rate_50_50 * N["TOTAL"]

    print(f"Biomnis combined cohort: N = {N['TOTAL']:,} (NEPHRO {N['NEPHRO']:,}, CARDIO {N['CARDIO']:,}, "
          f"AUTRE {N['AUTRE']:,}, NA {N['NA']:,})")
    print(f"Observed variants: {len(all_v)} (54 coding/canonical splice + IVS4+919G>A)")

    for sfx in ["raw", "dedup"]:
        kn = int(es[f"NEPHRO_{sfx}"].sum())
        kc = int(es[f"CARDIO_{sfx}"].sum())
        ka = int(es[f"AUTRE_{sfx}"].sum())
        kt = int(es[f"TOT_{sfx}"].sum())
        o_c, p_c = fisher_exact([[kc, N["CARDIO"] - kc], [ka, N["AUTRE"] - ka]])
        o_n, p_n = fisher_exact([[kn, N["NEPHRO"] - kn], [ka, N["AUTRE"] - ka]])
        p_pois = poisson.sf(kt - 1, expected)
        print(f"  Expert Strategy [{sfx:5s}]: TOT={kt} (1/{N['TOTAL']/kt:,.0f}, vs exp {expected:.1f}, Poisson p={p_pois:.4f}) | "
              f"CARDIO={kc} (1/{N['CARDIO']/kc:,.0f}, OR={o_c:.2f}, p={p_c:.4f}) | "
              f"NEPHRO={kn} (1/{N['NEPHRO']/kn:,.0f}, OR={o_n:.2f}, p={p_n:.4f})")

    det = all_v[all_v["act_mean"].notna()]
    rho_raw, p_raw = spearmanr(det["act_mean"], det["TOT_raw"])
    rho_ded, p_ded = spearmanr(det["act_mean"], det["TOT_dedup"])
    print(f"  Activity gradient Spearman (n={len(det)} variants): raw rho={rho_raw:.3f} (p={p_raw:.4f}) | "
          f"dedup rho={rho_ded:.3f} (p={p_ded:.4f})")

    es.to_csv(cfg.TABLES_DIR / "Biomnis_ExpertStrategy_carriers.csv", index=False)
    all_v.to_csv(cfg.TABLES_DIR / "Biomnis_all_coding_GLA.csv", index=False)


def main():
    cfg.ensure_dirs()
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["aggregated", "raw"], default="aggregated")
    parser.add_argument("--illumina-tsv", type=Path, default=cfg.ROOT / "data/raw/GLA_TV2.2026-11-16.GENOTYPES.tsv")
    parser.add_argument("--mgi-tsv", type=Path, default=cfg.ROOT / "data/raw/GLA_MGI_TV2.2026-11-20.GENOTYPES.tsv")
    args = parser.parse_args()

    if args.mode == "raw":
        if not args.illumina_tsv.exists() or not args.mgi_tsv.exists():
            print("Controlled-access raw TSV files not found locally; falling back to aggregated counts.")
            es = pd.read_csv(cfg.BIOMNIS_ES_CSV)
            all_v = pd.read_csv(cfg.BIOMNIS_ALL_CSV)
        else:
            es, all_v = run_raw(args.illumina_tsv, args.mgi_tsv)
    else:
        es = pd.read_csv(cfg.BIOMNIS_ES_CSV)
        all_v = pd.read_csv(cfg.BIOMNIS_ALL_CSV)

    summarize_aggregated(es, all_v)


if __name__ == "__main__":
    main()
