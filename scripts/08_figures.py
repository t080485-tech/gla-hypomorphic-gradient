#!/usr/bin/env python3
"""
08_figures.py — Generate all seven publication figures (600 dpi) from public
gnomAD v4.1.2 tables, BRaVa summaries, Tenon aggregate coordinates, and
variant-level Biomnis counts.
"""
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config as cfg

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 100,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
})


def make_fig1(grad: pd.DataFrame, out_dir: Path):
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    colors = {
        "Null (<1%)": "#67000d",
        "Severe (1-10%)": "#cb181d",
        "Threshold (10-25%)": "#f16913",
        "Moderate hypomorph (25-50%)": "#238b45",
        "Pseudodeficiency (>50%)": "#08519c",
    }
    for cls, col in colors.items():
        sub = grad[(grad["class"] == cls) & (grad["AC"] > 0)]
        ax.scatter(sub["act_mean"], sub["AF"] * 100, c=col, s=70, edgecolors="black",
                   linewidths=0.6, zorder=4, label=cls)
        sub0 = grad[(grad["class"] == cls) & (grad["AC"] == 0)]
        if len(sub0):
            ax.scatter(sub0["act_mean"], [0.0011] * len(sub0), c=col, s=42, marker="x", zorder=3)

    off = {
        "p.R112H": (8, 6), "p.Y216C": (-30, 8), "p.I232T": (-34, 10), "p.R356W": (8, -4),
        "p.R301G": (8, 6), "p.R118C": (-32, -16), "p.A143T": (8, 8), "p.R363H": (-38, -4),
        "p.N215S": (8, 10), "p.P214L": (8, -14), "p.R356Q": (8, 6), "p.I198T": (8, 6),
        "p.E66Q": (8, 6), "p.D83N": (-32, 8), "p.D175E": (8, 6), "p.S126G": (8, 8),
        "p.D313Y": (-14, 10), "p.A368T": (8, -14),
    }
    for _, r in grad[grad["AC"] > 0].iterrows():
        dx, dy = off.get(r["protein"], (7, 5))
        ax.annotate(r["protein"], (r["act_mean"], r["AF"] * 100), textcoords="offset points",
                    xytext=(dx, dy), fontsize=7.2, style="italic")

    ax.axvspan(20, 25, color="grey", alpha=0.15, zorder=1)
    ax.text(22.5, 0.75, "substrate accumulation\nthreshold (20-25%)", ha="center",
            fontsize=7.5, color="dimgrey", style="italic")
    ax.axhline(0.01, color="grey", ls=":", lw=0.9)
    ax.text(104, 0.0125, "MAF 0.01%", ha="right", fontsize=7.5, color="dimgrey")
    ax.axhline(0.00316, color="#cb181d", ls="--", lw=1.1)
    ax.text(104, 0.0037, "pathogenicity ceiling 3.2x10$^{-5}$", ha="right", fontsize=7.5,
            color="#cb181d", fontweight="bold")
    ax.set_yscale("log")
    ax.set_xlim(-4, 110)
    ax.set_ylim(0.0007, 1.6)
    ax.set_xlabel("Residual $\\alpha$-Gal A activity in vitro (% wild type, mean of two assays)")
    ax.set_ylabel("Allele frequency in gnomAD v4.1.2 (%)")
    ax.set_title("Residual enzyme activity predicts allele frequency\n"
                 "Spearman $\\rho$ = 0.750, p = 1x10$^{-4}$ (n = 21 observed variants)", fontsize=10)
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    fig.savefig(out_dir / "Fig1_hypomorphic_gradient.png")
    plt.close()


def make_fig2(out_dir: Path):
    data = [
        ("Clinical registries\n(historical)", 40000, "#9e9e9e"),
        ("Clinical registries\n(upper bound)", 170000, "#9e9e9e"),
        ("NBS Taiwan 2009\n(Hwu et al.)", 1250, "#f4a582"),
        ("NBS Italy 2006\n(Spada et al.)", 4600, "#f4a582"),
        ("NBS Japan\n(Sawada et al.)", 11854, "#f4a582"),
        ("NBS revised after\nreclassification (Germain)", 6883, "#2166ac"),
        ("UK Biobank 469k\n(Giammanco et al.)", 6938, "#2166ac"),
        ("gnomAD v4.1.2\nModel A (this study)", 6949, "#2166ac"),
        ("gnomAD v4.1.2\nExpert Strategy (this study)", 6340, "#b2182b"),
        ("In silico curation\n(this study)", 408, "#7f0000"),
    ]
    fig, ax = plt.subplots(figsize=(7.0, 5.4))
    y = np.arange(len(data))[::-1]
    for yi, (_, v, c) in zip(y, data):
        ax.hlines(yi, 0, np.log10(v), color=c, lw=2.4, alpha=0.9)
        ax.scatter(np.log10(v), yi, s=100, c=c, edgecolors="black", linewidths=0.6, zorder=4)
        ax.text(np.log10(v) + 0.07, yi, f"1/{v:,}", va="center", fontsize=8.5, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels([d[0] for d in data], fontsize=8)
    ax.axvspan(np.log10(6340), np.log10(6949), color="#2166ac", alpha=0.12, zorder=1)
    ax.text(np.log10(6640), 6.9, "convergence", fontsize=8, color="#2166ac", ha="center",
            style="italic", rotation=90)
    ax.set_xlim(0, 5.5)
    ax.set_xticks([0, 1, 2, 3, 4, 5])
    ax.set_xticklabels(["1/1", "1/10", "1/100", "1/1,000", "1/10,000", "1/100,000"], fontsize=8.5)
    ax.set_xlabel("Prevalence in males (1 in X, log scale)")
    ax.set_title("Twenty years of prevalence divergence, resolved by expert curation\n"
                 "three independent platforms converge at 1 in 6,300-6,950 males", fontsize=10)
    fig.savefig(out_dir / "Fig2_prevalence_waterfall.png")
    plt.close()


def make_fig3(anc: pd.DataFrame, gla: pd.DataFrame, out_dir: Path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.6), gridspec_kw={"width_ratios": [1.25, 1]})
    anc_v = anc.dropna(subset=["male_prev_1inx"]).sort_values("male_prev_1inx")
    short = {
        "African/African American": "African/African\nAmerican",
        "Admixed American": "Admixed\nAmerican",
        "Ashkenazi Jewish": "Ashkenazi\nJewish",
        "East Asian": "East Asian*",
        "European (non-Finnish)": "European\n(non-Finnish)",
        "South Asian": "South Asian",
        "Remaining": "Remaining",
    }
    labels = [short.get(a, a) for a in anc_v["ancestry"]]
    vals = anc_v["male_prev_1inx"].values
    cols = ["#b2182b" if v < 4000 else ("#ef8a62" if v < 12000 else "#2166ac") for v in vals]
    bars = ax1.barh(labels, vals, color=cols, edgecolor="black", linewidth=0.5)
    for b, v in zip(bars, vals):
        ax1.text(b.get_width() * 1.07, b.get_y() + b.get_height() / 2, f"1/{v:,.0f}",
                 va="center", fontsize=8, fontweight="bold")
    ax1.axvline(6340, color="grey", ls="--", lw=0.9)
    ax1.text(7300, 6.45, "global 1/6,340", fontsize=7, color="dimgrey")
    ax1.set_xscale("log")
    ax1.set_xlim(700, 90000)
    ax1.set_xlabel("Male prevalence (1 in X, log scale)", fontsize=9)
    ax1.set_title("A  Male prevalence by ancestry", fontsize=10, loc="left")
    ax1.tick_params(axis="y", labelsize=7.5)
    ax1.text(0.98, 0.03, "no carriers: Finnish, Middle Eastern, Amish",
             transform=ax1.transAxes, ha="right", fontsize=7, style="italic", color="dimgrey")

    ancs = ["European (non-Finnish)", "South Asian", "Admixed American", "African/African American", "East Asian"]
    variants = [
        ("X-101398432-C-A", "p.D313Y"), ("X-101403828-G-A", "p.R118C"),
        ("X-101401752-C-T", "p.A143T"), ("X-101401803-T-C", "p.S126G"),
        ("X-101403984-C-G", "p.E66Q"),
    ]
    x = np.arange(len(ancs))
    w = 0.155
    cols5 = ["#455a64", "#7b1fa2", "#0277bd", "#00897b", "#c62828"]
    gidx = gla.set_index("gnomAD ID")
    for i, (gid, lab) in enumerate(variants):
        afs = []
        for a in ancs:
            ac = pd.to_numeric(gidx.loc[gid, f"Allele Count {a}"], errors="coerce")
            an = pd.to_numeric(gidx.loc[gid, f"Allele Number {a}"], errors="coerce")
            afs.append((ac / an * 100) if (pd.notna(ac) and pd.notna(an) and an > 0) else 0)
        ax2.bar(x + (i - 2) * w, afs, w, label=lab, color=cols5[i], edgecolor="black", linewidth=0.35)
    ax2.set_xticks(x)
    ax2.set_xticklabels(["European\n(non-Finnish)", "South\nAsian", "Admixed\nAmerican",
                         "African/\nAmerican", "East\nAsian"], fontsize=7.5)
    ax2.set_ylabel("Allele frequency (%)", fontsize=9)
    ax2.set_title("B  Ancestry structure of the in silico trap\n(five reclassified benign alleles)",
                  fontsize=10, loc="left")
    ax2.legend(fontsize=7.2, frameon=False)
    fig.savefig(out_dir / "Fig3_ancestry.png")
    plt.close()


def make_fig4(miss: pd.DataFrame, out_dir: Path):
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    cats = ["Pathogenic/LP", "Conflicting", "VUS", "Benign/LB"]
    cat_colors = {"Pathogenic/LP": "#b2182b", "Conflicting": "#f4a582", "VUS": "#bdbdbd", "Benign/LB": "#2166ac"}
    np.random.seed(7)
    for i, c in enumerate(cats):
        sub = miss[miss["cat"] == c]
        xs = i + np.random.uniform(-0.2, 0.2, len(sub))
        ys = sub["Allele Frequency"] * 100
        ss = 18 + np.minimum(sub["Hemizygote Count"].fillna(0), 80) * 1.6
        ax.scatter(xs, ys, s=ss, c=cat_colors[c], alpha=0.6, edgecolors="black", linewidths=0.4, zorder=3)

    ax.annotate("p.N215S — the ceiling", (0, 0.00316), xytext=(0.15, 0.0055),
                fontsize=7.5, style="italic", arrowprops=dict(arrowstyle="-", lw=0.6, color="grey"))
    ax.annotate("p.R118C · p.A143T · p.S126G\n(reclassified benign)", (1, 0.058),
                xytext=(1.15, 0.16), fontsize=7.5, style="italic",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="grey"))
    ax.annotate("p.I198T\n(validated late-onset)", (1, 0.0143), xytext=(0.35, 0.030),
                fontsize=7.5, style="italic", arrowprops=dict(arrowstyle="-", lw=0.6, color="grey"))
    ax.annotate("p.E66Q (benign)", (1, 0.0132), xytext=(1.28, 0.0075),
                fontsize=7.5, style="italic", arrowprops=dict(arrowstyle="-", lw=0.6, color="grey"))
    ax.annotate("p.D313Y (pseudodeficiency)", (3, 0.3779), xytext=(2.35, 0.65),
                fontsize=7.5, style="italic", arrowprops=dict(arrowstyle="-", lw=0.6, color="grey"))
    ax.axhline(0.00316, color="#b2182b", ls="--", lw=1.2)
    ax.text(3.48, 0.0042, "ceiling 3.2x10$^{-5}$", fontsize=8, color="#b2182b", ha="right", fontweight="bold")
    ax.axhspan(0.00316, 0.01, color="grey", alpha=0.10)
    ax.text(2.0, 0.0055, "no man's land", fontsize=7.5, color="dimgrey", style="italic", ha="center")
    ax.set_yscale("log")
    ax.set_ylim(0.0006, 1.8)
    ax.set_xticks(range(4))
    ax.set_xticklabels([f"{c}\n(n={(miss['cat']==c).sum()})" for c in cats], fontsize=8.5)
    ax.set_xlabel("ClinVar classification of observed missense variants")
    ax.set_ylabel("Allele frequency in gnomAD v4.1.2 (%)")
    ax.set_title("The frequency ceiling of pathogenicity in GLA\n"
                 "no validated pathogenic variant exceeds AF 3.2x10$^{-5}$; "
                 "every variant above 10$^{-4}$ is benign", fontsize=10)
    fig.savefig(out_dir / "Fig4_frequency_ceiling.png")
    plt.close()


def make_fig5(miss: pd.DataFrame, out_dir: Path):
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(10.5, 5.0), gridspec_kw={"width_ratios": [1.2, 1]})
    plp_miss = miss[miss["cat"] == "Pathogenic/LP"]
    np.random.seed(7)
    axA.scatter(plp_miss["Allele Frequency"] * 100, np.random.uniform(-0.16, 0.16, len(plp_miss)),
                s=26, c="#bdbdbd", edgecolors="black", linewidths=0.3, zorder=2,
                label=f"gnomAD P/LP missense (n={len(plp_miss)})")
    axA.axvspan(2e-6, 3e-5, color="#f4a582", alpha=0.18, zorder=1)
    axA.text(7.5e-6, 0.55, "private variants — not observed\nin gnomAD (allele count = 0)",
             fontsize=7.3, ha="center", style="italic", color="#8c510a")
    tenon_vars = [
        ("p.N215S (T1-07)", 3.157e-5, 0.34),
        ("p.R363H (T1-08)", 2.566e-5, -0.40),
        ("p.G395A (T1-06)", 1.324e-5, 0.12),
        ("p.R112C (T1-05)", None, 0.40),
        ("p.G261D (T1-01)", None, 0.24),
        ("p.I219L (T1-02)", None, 0.06),
        ("p.F273S (T1-04)", None, -0.10),
        ("p.A143V (T1-03)", None, -0.26),
    ]
    priv_x = 1.05e-5
    for lab, af, yy in tenon_vars:
        xval = (af * 100) if af is not None else priv_x
        axA.scatter(xval, yy, s=130, c="#b2182b", marker="D", edgecolors="black",
                    linewidths=0.7, zorder=5, alpha=1.0 if af is not None else 0.85)
        axA.annotate(lab, (xval, yy), textcoords="offset points",
                     xytext=(8, 2) if af is not None else (9, -2), fontsize=7.4, fontweight="bold")
    axA.axvline(3.157e-3, color="#b2182b", ls="--", lw=1.1)
    axA.text(3.4e-3, -0.58, "ceiling", fontsize=7.5, color="#b2182b", rotation=90)
    axA.set_xscale("log")
    axA.set_xlim(2e-6, 3e-2)
    axA.set_ylim(-0.68, 0.72)
    axA.set_yticks([])
    axA.set_xlabel("Allele frequency in gnomAD v4.1.2 (%)", fontsize=9)
    axA.set_title("A  All eight Tenon genotypes on the population map", fontsize=10, loc="left")
    axA.spines["left"].set_visible(False)
    axA.legend(fontsize=7.2, frameon=False, loc="lower left", bbox_to_anchor=(0.005, 0.01))

    pts = [("04", 27.0, 13.8), ("05", 139.0, 4.6), ("07", 29.0, 2.7), ("08", 8.0, 7.04)]
    axB.axvspan(0, 20, color="#fddbc7", alpha=0.5, zorder=1)
    axB.axhspan(0, 0.8, color="#fddbc7", alpha=0.5, zorder=1)
    axB.axvline(20, color="#d95f02", lw=0.9, ls="--")
    axB.axhline(0.8, color="#d95f02", lw=0.9, ls="--")
    for pid, e, l in pts:
        axB.scatter(e, l, s=110, c="#2166ac", edgecolors="black", linewidths=0.7, zorder=4)
        axB.annotate(f"T1-{pid}", (e, l), textcoords="offset points", xytext=(8, 4),
                     fontsize=8, fontweight="bold")
    axB.scatter([3], [1.637], s=110, c="#bdbdbd", edgecolors="black", linewidths=0.7, zorder=4)
    axB.annotate("T1-03 (A143V)\nborderline 0.7-1.6", (3, 1.637), textcoords="offset points",
                 xytext=(8, 4), fontsize=7.5)
    axB.set_xscale("log")
    axB.set_yscale("log")
    axB.set_xlim(2, 300)
    axB.set_ylim(0.3, 30)
    axB.set_xlabel("Leukocyte $\\alpha$-Gal A activity (nmol/h/mg prot.; normal $\\geq$20)", fontsize=9)
    axB.set_ylabel("Plasma lyso-Gb3 (nmol/L; normal <0.8)", fontsize=9)
    axB.set_title("B  Enzyme misses, lyso-Gb3 reveals\n(seven females, Tenon cohort)", fontsize=10, loc="left")
    axB.text(45, 0.45, "normal enzyme", fontsize=7.5, color="#d95f02")
    axB.text(2.4, 13, "lyso-Gb3 elevated\nin 4 of 5 measured", fontsize=8, color="#2166ac", fontweight="bold")
    fig.savefig(out_dir / "Fig5_Tenon_cohort.png")
    plt.close()


def make_fig6(out_dir: Path):
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 6.4), gridspec_kw={"width_ratios": [1.55, 1]})
    gene_rows = [
        ("Gout · EAS · dam. missense",     0.2239,  0.06905, 0.001175, "#8c510a"),
        ("Gout · AFR · other missense",    0.07738, 0.02348, 0.002188, "#01665e"),
        ("HF · AFR · other missense",      0.07783, 0.02438, 0.002138, "#01665e"),
        ("Gout · SAS · other missense",    0.1518,  0.0608,  0.01259,  "#762a83"),
        ("CRF · SAS · pLoF",               0.2754,  0.1055,  0.00912,  "#762a83"),
        ("HF · SAS · pLoF",                0.2667,  0.1039,  0.01023,  "#762a83"),
        ("AFib · SAS · pLoF",              0.2592,  0.1025,  0.01148,  "#762a83"),
        ("CRF · EAS · pLoF|dam. missense", 0.1311,  0.05713, 0.02188,  "#8c510a"),
    ]
    y = np.arange(len(gene_rows))[::-1]
    for yi, (_, beta, se, p, col) in zip(y, gene_rows):
        lo, hi = beta - 1.96 * se, beta + 1.96 * se
        axA.plot([lo, hi], [yi, yi], color=col, lw=1.8, solid_capstyle="butt", zorder=3)
        axA.scatter(beta, yi, s=52, c=col, edgecolors="black", linewidths=0.5, zorder=4)
        axA.text(hi + 0.012, yi, f"p={p:.4f}" if p >= 0.001 else f"p={p:.3f}",
                 va="center", fontsize=7.2, color="dimgrey")
    axA.axvline(0, color="black", lw=0.8, zorder=1)
    axA.set_yticks(y)
    axA.set_yticklabels([r[0] for r in gene_rows], fontsize=8)
    axA.set_xlabel("Beta (Burden, inverse-variance meta-analysis)", fontsize=9)
    axA.set_title("A  Gene-level suggestive associations (p < 0.05)\n"
                  "none reaches the Bonferroni threshold (1.4 x 10$^{-7}$)", fontsize=10, loc="left")
    axA.set_xlim(-0.06, 0.62)
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=8,
                      markeredgecolor="black", markeredgewidth=0.5, label=l)
               for c, l in [("#762a83", "South Asian"), ("#01665e", "African"), ("#8c510a", "East Asian")]]
    axA.legend(handles=handles, fontsize=7.5, frameon=False, loc="lower right")

    var_rows = [
        ("p.Asn215Ser\nHeart failure", 2.03, 0.70, "#b2182b"),
        ("p.Asn215Ser\nValvular disease", 1.62, 0.771, "#b2182b"),
        ("p.Ala368Thr\nAtrial fibrillation (AFR)", 2.6, 0.881, "#455a64"),
        ("p.Arg118Cys\nHypertension", 0.203, 0.0888, "#2166ac"),
    ]
    y2 = np.arange(len(var_rows))[::-1]
    for yi, (_, beta, se, col) in zip(y2, var_rows):
        orv = np.exp(beta)
        lo, hi = np.exp(beta - 1.96 * se), np.exp(beta + 1.96 * se)
        axB.plot([lo, hi], [yi, yi], color=col, lw=1.8, solid_capstyle="butt", zorder=3)
        axB.scatter(orv, yi, s=62, c=col, edgecolors="black", linewidths=0.5, zorder=4)
        axB.text(hi * 1.12, yi, f"OR {orv:.1f}", va="center", fontsize=7.5, color="dimgrey")
    axB.axvline(1, color="black", lw=0.8, ls="--", zorder=1)
    axB.set_yticks(y2)
    axB.set_yticklabels([r[0] for r in var_rows], fontsize=8)
    axB.set_xscale("log")
    axB.set_xlabel("Odds ratio (95% CI, log scale)", fontsize=9)
    axB.set_title("B  Variant-level associations\n"
                  "the ceiling variant p.Asn215Ser carries the cardiac signal", fontsize=10, loc="left")
    axB.set_xlim(0.8, 60)
    axB.text(1.05, len(var_rows) - 0.45, "no effect", fontsize=7, color="dimgrey", rotation=90, va="top")
    plt.tight_layout()
    fig.savefig(out_dir / "Fig6_BRaVa_forest.png")
    plt.close()


def make_fig7(bio_all: pd.DataFrame, bio_es: pd.DataFrame, out_dir: Path):
    N = cfg.BIOMNIS_DENOMS
    classes_ord = ["Null (<1%)", "Severe (1-10%)", "Threshold (10-25%)",
                   "Moderate hypomorph (25-50%)", "Pseudodeficiency (>50%)"]
    classes_lbl = ["Null\n(<1%)", "Severe\n(1-10%)", "Threshold\n(10-25%)",
                   "Moderate\n(25-50%)", "Pseudodef.\n(>50%)"]

    nephro = np.array([int(bio_all[bio_all["class"] == c]["NEPHRO_raw"].sum()) for c in classes_ord])
    cardio = np.array([int(bio_all[bio_all["class"] == c]["CARDIO_raw"].sum()) for c in classes_ord])
    autre = np.array([int(bio_all[bio_all["class"] == c]["AUTRE_raw"].sum()) for c in classes_ord])
    na_arr = np.array([int(bio_all[bio_all["class"] == c]["NA_raw"].sum()) for c in classes_ord])
    totals = nephro + cardio + autre + na_arr

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(13.4, 4.8), gridspec_kw={"width_ratios": [1.25, 0.85, 1.05]})
    x = np.arange(len(classes_lbl))
    ax1.bar(x, nephro, 0.62, label=f"NEPHRO (n = {N['NEPHRO']:,})", color="#b2182b", edgecolor="black", linewidth=0.5)
    ax1.bar(x, cardio, 0.62, bottom=nephro, label=f"CARDIO (n = {N['CARDIO']:,})", color="#f4a582", edgecolor="black", linewidth=0.5)
    ax1.bar(x, autre + na_arr, 0.62, bottom=nephro + cardio, label=f"AUTRE / NA (n = {N['AUTRE']+N['NA']:,})",
            color="#2166ac", edgecolor="black", linewidth=0.5)
    for xi, t in zip(x, totals):
        ax1.text(xi, t * 1.18, f"{t}\n(1/{N['TOTAL']//t:,})", ha="center", fontsize=7.4, fontweight="bold")
    ax1.set_yscale("log")
    ax1.set_ylim(0.5, 1200)
    ax1.set_xticks(x)
    ax1.set_xticklabels(classes_lbl, fontsize=8)
    ax1.set_ylabel(f"Carriers in Biomnis (N = {N['TOTAL']:,})", fontsize=9)
    ax1.set_title("A  Activity gradient replicated (264-fold)\nSpearman $\\rho$ = 0.560, p = 0.047 (n = 13 variants)",
                  fontsize=9.5, loc="left")
    ax1.legend(fontsize=7.2, frameon=False, loc="upper left")

    k_ded = int(bio_es["TOT_dedup"].sum())
    k_raw = int(bio_es["TOT_raw"].sum())
    ax2.bar(["Expected\n(gnomAD)", "Observed\n(index)", "Observed\n(total)"], [7.7, k_ded, k_raw],
            0.56, color=["#bdbdbd", "#2166ac", "#b2182b"], edgecolor="black", linewidth=0.5)
    ax2.text(0, 8.1, "7.7", ha="center", fontsize=9.5, fontweight="bold")
    ax2.text(1, k_ded + 0.4, f"{k_ded}\n(p = 0.026)", ha="center", fontsize=8.3, fontweight="bold")
    ax2.text(2, k_raw + 0.4, f"{k_raw}\n(p = 0.006)", ha="center", fontsize=8.3, fontweight="bold")
    ax2.set_ylabel("Expert Strategy carriers (9 variants)", fontsize=9)
    ax2.set_ylim(0, 20.5)
    ax2.tick_params(axis="x", labelsize=8)
    ax2.set_title("B  Enrichment over gnomAD\n1.8- to 2.1-fold (Poisson)", fontsize=9.5, loc="left")

    ka = int(bio_es["AUTRE_dedup"].sum())
    kn = int(bio_es["NEPHRO_dedup"].sum())
    kc = int(bio_es["CARDIO_dedup"].sum())
    prev_per_10k = [ka / N["AUTRE"] * 10000, kn / N["NEPHRO"] * 10000, kc / N["CARDIO"] * 10000]
    ax3.bar(np.arange(3), prev_per_10k, 0.55, color=["#2166ac", "#ef8a62", "#b2182b"], edgecolor="black", linewidth=0.5)
    notes_c = [
        f"{ka} (1/{N['AUTRE']//ka:,})\nreference",
        f"{kn} (1/{N['NEPHRO']//kn:,})\nOR 1.66 (p=0.36)",
        f"{kc} (1/{N['CARDIO']//kc:,})\nOR 10.8 (p=0.021)",
    ]
    for xi, h, nt in zip(np.arange(3), prev_per_10k, notes_c):
        ax3.text(xi, h + 1.1, nt, ha="center", fontsize=7.5, fontweight="bold")
    ax3.axhline(1 / 4943 * 10000, color="grey", ls="--", lw=1.1, zorder=2, label="gnomAD baseline (1/4,943)")
    ax3.legend(fontsize=7.3, frameon=False, loc="upper left")
    ax3.set_xticks(np.arange(3))
    ax3.set_xticklabels([f"AUTRE\n(n = {N['AUTRE']:,})", f"NEPHRO\n(n = {N['NEPHRO']:,})", f"CARDIO\n(n = {N['CARDIO']:,})"], fontsize=8)
    ax3.set_ylabel("Index carriers per 10,000 patients", fontsize=9)
    ax3.set_ylim(0, 38)
    ax3.set_title("C  Phenotypic enrichment by stratum\ncardiac stratum enriched 10.8-fold", fontsize=9.5, loc="left")

    plt.tight_layout()
    fig.savefig(out_dir / "Fig7_Biomnis_gradient.png")
    plt.close()


def main():
    cfg.ensure_dirs()
    gla = pd.read_csv(cfg.GNOMAD_CSV, low_memory=False)
    for c in ["Allele Count", "Allele Frequency", "Hemizygote Count"]:
        gla[c] = pd.to_numeric(gla[c], errors="coerce")
    grad = pd.read_csv(cfg.TABLES_DIR / "GLA_hypomorphic_gradient.csv")
    anc = pd.read_csv(cfg.TABLES_DIR / "GLA_ancestry_ModelB.csv")
    bio_all = pd.read_csv(cfg.BIOMNIS_ALL_CSV)
    bio_es = pd.read_csv(cfg.BIOMNIS_ES_CSV)

    miss = gla[(gla["VEP Annotation"] == "missense_variant") & (gla["Allele Count"] > 0)].copy()

    def cat(r):
        cl = r["ClinVar Germline Classification"]
        if cl in cfg.PLP_LABELS:
            return "Pathogenic/LP"
        if cl in cfg.BENIGN_LABELS:
            return "Benign/LB"
        if cl == "Conflicting classifications of pathogenicity":
            return "Conflicting"
        if cl == "Uncertain significance":
            return "VUS"
        return "No ClinVar"

    miss["cat"] = miss.apply(cat, axis=1)

    make_fig1(grad, cfg.FIGURES_DIR)
    make_fig2(cfg.FIGURES_DIR)
    make_fig3(anc, gla, cfg.FIGURES_DIR)
    make_fig4(miss, cfg.FIGURES_DIR)
    make_fig5(miss, cfg.FIGURES_DIR)
    make_fig6(cfg.FIGURES_DIR)
    make_fig7(bio_all, bio_es, cfg.FIGURES_DIR)
    print(f"Generated all 7 figures (600 dpi) in {cfg.FIGURES_DIR}.")


if __name__ == "__main__":
    main()
