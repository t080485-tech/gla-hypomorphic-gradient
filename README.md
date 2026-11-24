# GLA Hypomorphic Gradient & Fabry Disease Prevalence Analysis

Code repository for:

**"The hypomorphic gradient of *GLA* revealed by population genetics: expert curation resolves the Fabry disease prevalence debate, exposes the frequency ceiling of pathogenicity, and redefines the diagnostic strategy for late-onset disease in chronic kidney disease"**

Thomas Robert, MD, PhD; Fanny Ponce, PharmD; Xavier Vanhoye, PharmD; Laure Raymond, PharmD; Laurent Mesnard, MD, PhD (2026)

---

## Overview

This repository provides the reproducible analytical pipeline supporting the manuscript. It evaluates the genetic prevalence of Fabry disease (*GLA*, GRCh38/hg38, `ENSG00000102393`, canonical transcript `ENST00000218516.4` / `NM_000169.3`) across population-scale sequencing, biobanks, and clinical exome cohorts.

Key analyses:

- **Three-model variant curation** in gnomAD v4.1.2 (807,186 individuals; 361,371 males; 425,246 females): ClinVar-strict (Model A), Expert Strategy (Model B), and *in silico* comparator (Model C)
- **Sex-stratified & X-linked carrier estimation** (empirical female carriers and `2 × ΣAF` sensitivity; classic pLoF vs. late-onset split)
- **Ancestry-stratified prevalence** across all 10 gnomAD ancestry groups and population structure of the five reclassified benign alleles
- **Hypomorphic activity gradient** (Spearman correlation between *in vitro* residual $\alpha$-Gal A activity and gnomAD allele frequency; empirical pathogenicity ceiling at AF $3.2 \times 10^{-5}$)
- **Biobank lookups** from the BRaVa consortium (~1.2M individuals) and OpenTargets
- **Clinical exome validation** in 38,117 consecutive patients sequenced at Eurofins Biomnis (Illumina $n = 36,184$ + MGI $n = 1,933$) across `NEPHRO` ($n = 11,161$), `CARDIO` ($n = 688$), and `AUTRE` ($n = 25,935$) strata, including probe-captured `IVS4+919G>A` and prescriber-based index deduplication
- **Publication figures** (Figures 1–7 at 600 dpi)

---

## Repository Structure

```text
gla-hypomorphic-gradient/
├── config.py                            # Genome build, transcript, curation thresholds & evidence tables
├── run_all.py                           # End-to-end execution script
├── requirements.txt                     # Python package dependencies
├── LICENSE                              # MIT License
├── data/
│   ├── public/
│   │   ├── gnomAD_v4.1.2_GLA_ENSG00000102393.csv
│   │   ├── Table_S11_BRaVa_GLA.csv
│   │   └── Table_S12_OpenTargets_GLA.csv
│   └── aggregated/                      # Anonymized variant-level Biomnis counts (no patient/prescriber IDs)
│       ├── Biomnis_ExpertStrategy_carriers.csv
│       └── Biomnis_all_coding_GLA.csv
├── scripts/
│   ├── 01_curation_models.py            # Models A, B (Expert Strategy), C on gnomAD v4.1.2
│   ├── 02_female_carriers.py            # Empirical + theoretical female carrier estimates
│   ├── 03_ancestry_stratification.py    # Ancestry-stratified prevalence & benign allele structure
│   ├── 04_hypomorphic_gradient.py       # Residual enzyme activity vs. AF & pathogenicity ceiling
│   ├── 05_biobank_brava_ot.py           # BRaVa & OpenTargets summaries
│   ├── 06_biomnis_validation.py         # Clinical exome enrichment & gradient (--mode aggregated | raw)
│   ├── 07_tenon_clinical_series.py      # [CONTROLLED ACCESS] Hôpital Tenon EXOME_TNN T1 vs. T3
│   └── 08_figures.py                    # Generate Figures 1-7 (600 dpi)
└── output/
    ├── tables/                          # Reproducible CSV outputs
    └── figures/                         # Generated 600-dpi PNG figures
```

---

## Data Availability & Privacy

- **Publicly executable**: Steps `01` through `05`, `06` (`--mode aggregated`), and `08` run end-to-end using the public gnomAD v4.1.2, BRaVa, OpenTargets, and anonymized variant-level Biomnis counts included in `data/`.
- **Controlled-access data**: Raw Eurofins Biomnis genotype TSV files (which contain internal dossier and prescriber identifiers) and pseudonymized individual-level Hôpital Tenon clinical workbooks (`EXOME_TNN` T1/T3) are **not** included in this repository (`data/raw/` is excluded via `.gitignore`). Step `07` automatically skips when local controlled-access files are absent.

---

## Quick Start

```bash
# 1. Create environment & install dependencies
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# 2. Run the full public + aggregated pipeline and regenerate all figures
python run_all.py
```

### Running controlled-access Biomnis raw parsing (authorized collaborators)

```bash
python scripts/06_biomnis_validation.py --mode raw \
  --illumina-tsv data/raw/GLA_TV2.2026-11-16.GENOTYPES.tsv \
  --mgi-tsv data/raw/GLA_MGI_TV2.2026-11-20.GENOTYPES.tsv
```

---

## Key Curation Parameters (`config.py`)

| Parameter | Value | Description |
|:---|:---|:---|
| Reference build | `GRCh38/hg38` | Genomic coordinate system |
| Gene / Canonical transcript | `GLA` (`ENSG00000102393`), `ENST00000218516.4` / `NM_000169.3` | Ensembl / RefSeq MANE Select |
| Male denominator | `361,371` | gnomAD v4.1.2 exome hemizygous males |
| Female denominator | `425,246` | Derived from `AN_max = 1,211,863` |
| Rare pLoF MAF threshold | `<= 0.0001` (0.01%) | High-confidence LOFTEE (`hc_lof`) |
| *In silico* comparator thresholds | `REVEL >= 0.5` or `CADD >= 20` | Excludes ClinVar Benign/Likely benign |
| Empirical pathogenicity ceiling | `AF = 3.16e-5` | `p.Asn215Ser` maximum among validated P/LP |

---

## Citation

If you use this code or curation framework, please cite:

> Robert T, Ponce F, Vanhoye X, Raymond L, Mesnard L. *The hypomorphic gradient of GLA revealed by population genetics: expert curation resolves the Fabry disease prevalence debate, exposes the frequency ceiling of pathogenicity, and redefines the diagnostic strategy for late-onset disease in chronic kidney disease.* (2026, submitted).

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
