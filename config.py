"""Shared parameters, genomic coordinates, curation rules, and paths for GLA analysis."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA_PUBLIC = ROOT / "data" / "public"
DATA_AGGREGATED = ROOT / "data" / "aggregated"
DATA_RAW = ROOT / "data" / "raw"
TABLES_DIR = ROOT / "output" / "tables"
FIGURES_DIR = ROOT / "output" / "figures"

# Reference genome & transcript
GENOME_BUILD = "GRCh38"
GENE_SYMBOL = "GLA"
ENSEMBL_GENE_ID = "ENSG00000102393"
CANONICAL_TRANSCRIPT_ENSEMBL = "ENST00000218516.4"
CANONICAL_TRANSCRIPT_REFSEQ = "NM_000169.3"
PROTEIN_ID = "NP_000160.1"

# Input files
GNOMAD_CSV = DATA_PUBLIC / "gnomAD_v4.1.2_GLA_ENSG00000102393.csv"
BRAVA_CSV = DATA_PUBLIC / "Table_S11_BRaVa_GLA.csv"
OPENTARGETS_CSV = DATA_PUBLIC / "Table_S12_OpenTargets_GLA.csv"

BIOMNIS_ES_CSV = DATA_AGGREGATED / "Biomnis_ExpertStrategy_carriers.csv"
BIOMNIS_ALL_CSV = DATA_AGGREGATED / "Biomnis_all_coding_GLA.csv"

# gnomAD v4.1.2 denominators (X-chromosome exomes/genomes)
N_MALES = 361371
AN_MAX = 1211863
N_FEMALES = int((AN_MAX - N_MALES) / 2)  # 425,246
GLOBAL_MALE_FRAC = N_MALES / 807186  # 0.44744

ANCESTRIES = [
    "African/African American",
    "Admixed American",
    "Ashkenazi Jewish",
    "East Asian",
    "European (Finnish)",
    "Middle Eastern",
    "European (non-Finnish)",
    "Amish",
    "South Asian",
    "Remaining",
]

# Biomnis combined cohort denominators (Illumina 36,184 + MGI 1,933)
BIOMNIS_DENOMS = {
    "TOTAL": 38117,
    "NEPHRO": 11161,
    "CARDIO": 688,
    "AUTRE": 25935,
    "NA": 333,
}

# Variant curation thresholds
PLP_LABELS = ["Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic"]
PLOF_VEP = [
    "frameshift_variant",
    "stop_gained",
    "splice_acceptor_variant",
    "splice_donor_variant",
]
CODING_VEP = PLOF_VEP + [
    "missense_variant",
    "stop_lost",
    "inframe_deletion",
    "inframe_insertion",
]
BENIGN_LABELS = [
    "Benign",
    "Benign/Likely benign",
    "Benign/Likely benign; other",
    "Likely benign",
]
MAF_MAX = 0.0001  # 0.01%
REVEL_MIN = 0.5
CADD_MIN = 20.0
FREQUENCY_CEILING_AF = 3.15663645433e-5  # p.Asn215Ser

# Literature evidence table for key contested/late-onset GLA variants
EVIDENCE = {
    "X-101398432-C-A": (
        "BENIGN",
        "p.D313Y pseudodeficiency (~75% residual activity, normal lyso-Gb3)",
        "Effraimidis 2020 PMID 32246457; Germain 2022 PMID 34927718",
    ),
    "X-101401752-C-T": (
        "BENIGN",
        "p.A143T reclassified benign (25-72% activity, normal lyso-Gb3, null UKBB)",
        "Germain 2022 PMID 34927718; Giammanco 2025 PMID 41111061",
    ),
    "X-101401803-T-C": (
        "BENIGN",
        "p.S126G reclassified benign (cellular characterization, null UKBB)",
        "Breyer 2024 PMID 38469097; Giammanco 2025 PMID 41111061",
    ),
    "X-101403984-C-G": (
        "BENIGN",
        "p.E66Q benign East Asian polymorphism (0.83% Japan, normal lyso-Gb3)",
        "Sawada 2025 PMID 40823509; Togawa 2012 PMID 22305854",
    ),
    "X-101403828-G-A": (
        "BENIGN",
        "p.R118C reclassified benign (expert panel, 24.5% activity, null UKBB)",
        "Germain 2022 PMID 34927718; Giammanco 2025 PMID 41111061",
    ),
    "X-101398032-C-T": (
        "VUS_CONTESTED",
        "p.R356Q negative biomarkers (Sawada 2025), 89% HEK293 activity",
        "Sawada 2025 PMID 40823509; Gilchrist 2023 PMID 35977816",
    ),
    "X-101399747-C-T": (
        "LATE_ONSET_VALIDATED",
        "IVS4+919G>A (c.640-801G>A) validated late-onset cardiac founder variant (ClinVar Pathogenic 2-star; genome-only in gnomAD)",
        "Hwu 2009 PMID 19621417; Chien 2012 PMID 22437327",
    ),
    "X-101400712-A-G": (
        "LATE_ONSET_VALIDATED",
        "p.I198T validated mild late-onset cardiac variant (UKBB, Lantern, VVCC)",
        "Gilchrist 2023 PMID 35977816; Pillai 2025 PMID 40840051",
    ),
}

# Residual alpha-Gal A activity (% WT baseline, Lukas et al. 2020 + Togawa 2012)
RESIDUAL_ACTIVITY = {
    "X-101398432-C-A": ("p.D313Y", 83.9, 59.0, "Lukas 2020"),
    "X-101397997-C-T": ("p.A368T", 103.7, 54.6, "Lukas 2020"),
    "X-101403933-C-T": ("p.D83N", 62.9, 69.2, "Lukas 2020"),
    "X-101401654-G-C": ("p.D175E", 89.8, 44.3, "Lukas 2020"),
    "X-101403984-C-G": ("p.E66Q", 55.0, 55.0, "Togawa 2012"),
    "X-101401803-T-C": ("p.S126G", 51.3, 83.7, "Lukas 2020"),
    "X-101401752-C-T": ("p.A143T", 31.3, 21.4, "Lukas 2020"),
    "X-101398942-T-C": ("p.N215S", 36.7, 15.6, "Lukas 2020"),
    "X-101400712-A-G": ("p.I198T", 38.7, 64.7, "Lukas 2020"),
    "X-101398011-C-T": ("p.R363H", 31.9, 20.0, "Lukas 2020"),
    "X-101398032-C-T": ("p.R356Q", 33.3, 36.1, "Lukas 2020"),
    "X-101403828-G-A": ("p.R118C", 20.0, 24.0, "Lukas 2020"),
    "X-101398468-G-C": ("p.R301G", 19.3, 19.1, "Lukas 2020"),
    "X-101398033-G-A": ("p.R356W", 16.9, 11.0, "Lukas 2020"),
    "X-101398891-A-G": ("p.I232T", 11.5, 15.0, "Lukas 2020"),
    "X-101398945-G-A": ("p.P214L", 19.4, 33.0, "Lukas 2020"),
    "X-101403845-C-T": ("p.R112H", 1.6, 2.6, "Lukas 2020"),
    "X-101398939-T-C": ("p.Y216C", 2.3, 2.0, "Lukas 2020"),
    "X-101400682-A-T": ("p.M208K", 2.7, None, "Lukas 2020"),
    "X-101398558-C-A": ("p.G271C", 0.2, 0.0, "Lukas 2020"),
    "X-101398074-C-T": ("p.R342Q", 0.0, 0.0, "Lukas 2020"),
    "X-101398068-A-G": ("p.L344P", 0.0, 0.0, "Lukas 2020"),
    "X-101398390-G-T": ("p.Q327K", 0.0, 0.0, "Lukas 2020"),
}


def ensure_dirs() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
