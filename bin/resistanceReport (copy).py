#!/usr/bin/env python3

# ============================================================
# TB Resistance Final Report Generator
#
# Usage:
#   python3 resistanceReport.py <biosample ID>
#
# Input:
#   resistance/<biosample ID>/<biosample ID>_OMStarget.xlsx
#       → merged variant list from all callers (GATK, NORM, LOFREQ, DELLY)
#         including PASS/FAIL, AF, zygosity, annotation, and OMS catalog matches
#
# Output:
#   results/resistance/<biosample ID>.xlsx
#       → Final deduplicated and curated resistance report
#
# Rules:
#   1. Keep both PASS/FAIL variants.
#   2. Remove decomposed SNPs inside the span of a MNP
#         (e.g. CA→TG vs C→T and A→G after vt norm).
#   3. If multiple variants overlap a position:
#         - If any is HET → keep all.
#         - If MNP vs SNP → keep MNP, remove SNP.
#   4. Caller priority for duplicates:
#         GATK > NORM > LOFREQ > DELLY
#   5. Evidence values converted to:
#         R / r / u / s / S
#
# ============================================================

import sys
import pandas as pd
from pathlib import Path

# ============================================================
# PROJECT STRUCTURE (projectDir-aware)
# ============================================================
PROJECT_DIR = Path(__file__).resolve().parent.parent

INPUT_DIR  = PROJECT_DIR / "resistance"
OUTPUT_DIR = PROJECT_DIR / "results" / "resistance"

CALLER_PRIORITY = ["GATK", "NORM", "LOFREQ", "DELLY"]

FINAL_COLUMNS = [
    "Drug","Gene","Tier","Variant","Effect","Evidence","Comment",
    "AF","ALT_READS","Heteroresistance","Caller",
    "Filter_Status","Filter_Method"
]

EVIDENCE_MAP = {
    "Assoc w R": "R",
    "Assoc w R - Interim": "r",
    "Uncertain significance": "u",
    "Not assoc w R - Interim": "s",
    "Not assoc w R": "S"
}

# ============================================================
# HELPERS
# ============================================================
def convert_evidence(raw):
    if not isinstance(raw, str):
        return "S"
    if ") " in raw:
        raw = raw.split(") ", 1)[1].strip()
    return EVIDENCE_MAP.get(raw, "S")

def choose_best_caller(group):
    for caller in CALLER_PRIORITY:
        sel = group[group["caller"] == caller]
        if not sel.empty:
            return sel.iloc[0]
    return group.iloc[0]

def is_mnp(ref, alt):
    return len(str(ref)) > 1 or len(str(alt)) > 1

# ============================================================
# MNP RESOLUTION LOGIC
# ============================================================
def resolve_mnp_span(df):

    to_keep = set()
    to_remove = set()

    mnp_rows = df[(df["ref"].str.len() > 1) | (df["alt"].str.len() > 1)]

    for mnp_idx, mnp in mnp_rows.iterrows():
        mnp_start = mnp["position"]
        mnp_end = mnp["position"] + len(str(mnp["ref"])) - 1
        mnp_is_het = (mnp["zygosity"] == "HET")

        snps = df[
            (df["position"] >= mnp_start) &
            (df["position"] <= mnp_end) &
            ~((df["ref"].str.len() > 1) | (df["alt"].str.len() > 1))
        ]

        if snps.empty:
            to_keep.add(mnp_idx)
            continue

        snps_het = snps[snps["zygosity"] == "HET"]
        snps_hom = snps[snps["zygosity"] == "HOM"]

        if not mnp_is_het and not snps_hom.empty and snps_het.empty:
            to_keep.add(mnp_idx)
            to_remove.update(snps.index)

        elif mnp_is_het and not snps_hom.empty and snps_het.empty:
            to_keep.add(mnp_idx)
            to_remove.update(snps.index)

        elif not mnp_is_het and not snps_het.empty and snps_hom.empty:
            to_keep.add(mnp_idx)
            to_keep.update(snps_het.index)
            to_remove.update(snps_hom.index)

        elif mnp_is_het and not snps_het.empty:
            to_keep.add(mnp_idx)
            to_remove.update(snps.index)

        else:
            to_keep.add(mnp_idx)
            if not mnp_is_het:
                to_keep.update(snps_het.index)
            to_remove.update(snps_hom.index)

    final_idx = set(df.index) - to_remove
    return df.loc[sorted(final_idx)].copy()


def resolve_mnp_vs_snp(df):
    result = []

    for pos, group in df.groupby("position"):

        if any(group["zygosity"] == "HET"):
            result.append(group)
            continue

        mnp_mask = (group["ref"].str.len() > 1) | (group["alt"].str.len() > 1)
        snp_mask = ~mnp_mask

        if mnp_mask.any() and snp_mask.any():
            result.append(group[mnp_mask])
        else:
            result.append(group)

    return pd.concat(result, ignore_index=True)

# ============================================================
# FAST VARIANT RESOLUTION (vectorized)
# ============================================================
def resolve_variant(df):

    aa = df["aa_change"]
    master = df["master_change"]
    nt = df["nt_change"]

    # --------------------------------------------------------
    # Extract genes (vectorized)
    # --------------------------------------------------------
    gene_master = master.str.split("_", n=1).str[0]
    gene_aa = aa.str.split("_", n=1).str[0]
    gene_nt = nt.str.split("_", n=1).str[0]

    gene_master = gene_master.where(master != "NA")
    gene_aa = gene_aa.where(aa != "NA")
    gene_nt = gene_nt.where(nt != "NA")

    # --------------------------------------------------------
    # Detect gene mismatch with catalog
    # --------------------------------------------------------
    mismatch = (
        (gene_aa.notna() & (gene_aa != gene_master)) |
        (gene_nt.notna() & (gene_nt != gene_master))
    )

    # --------------------------------------------------------
    # Apply variant priority
    # aa_change > master_change > nt_change
    # --------------------------------------------------------
    variant = aa.copy()

    mask = variant == "NA"
    variant[mask] = master[mask]

    mask = variant == "NA"
    variant[mask] = nt[mask]

    # --------------------------------------------------------
    # Force catalog variant if gene mismatch
    # --------------------------------------------------------
    variant[mismatch] = master[mismatch]

    return variant

# ============================================================
# CORE PROCESS
# ============================================================
def process_biosample(biosample):

    input_xlsx  = INPUT_DIR / biosample / f"{biosample}_OMStarget.xlsx"
    output_xlsx = OUTPUT_DIR / f"{biosample}.xlsx"

    # ===================== SKIP =====================
    if output_xlsx.exists():
        print(f"[SKIP] Final resistance report already exists → {output_xlsx}")
        return

    # ===================== FAIL FAST =====================
    if not input_xlsx.exists():
        print(f"[ERROR] Required input not found: {input_xlsx}")
        sys.exit(1)

    print(f"[RUN] Generating resistance report for {biosample}")
    df = pd.read_excel(input_xlsx)

    if df.empty:
        print(f"[ERROR] OMStarget is empty: {input_xlsx}")
        sys.exit(1)

    # ===================== CLEANING =====================
    df = resolve_mnp_span(df)
    df = resolve_mnp_vs_snp(df)

    # ===================== COLUMN MAPPING =====================
    df["Drug"]  = df["drug"].astype(str)
    df["Gene"]  = df["gene"].astype(str)
    df["Tier"]  = df["tier"]

    # Normalize NA
    df[["aa_change","master_change","nt_change"]] = df[
        ["aa_change","master_change","nt_change"]
    ].fillna("NA")

    # Resolve change annotation priority
    df["Variant"] = resolve_variant(df)

    df["Effect"]  = df["effect"]
    df["Evidence"] = df["FINAL CONFIDENCE GRADING"].apply(convert_evidence)
    df["AF"]       = df["AF"]
    df["ALT_READS"] = df["ALT_reads"]
    df["Heteroresistance"] = df["zygosity"]
    df["Caller"]   = df["caller"]
    df["Filter_Status"] = df["Filter_Status"]
    df["Filter_Method"] = df["Filter_Method"]

    # ===================== DEDUPLICATION =====================
    final_rows = []
    for _, group in df.groupby(["Drug", "Gene", "Variant"], dropna=False):
        final_rows.append(choose_best_caller(group))

    final = pd.DataFrame(final_rows)[FINAL_COLUMNS]

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    final.to_excel(output_xlsx, index=False)

    print(f"[OK] Final resistance report written → {output_xlsx}")

# ============================================================
def main():
    if len(sys.argv) != 2:
        print("Usage: python3 resistanceReport.py <biosample>")
        sys.exit(1)

    process_biosample(sys.argv[1])

if __name__ == "__main__":
    main()


