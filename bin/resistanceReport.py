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

def deduplicate_variants(df):

    # garantir tipo consistente para agrupamento
    df["alt"] = df["alt"].astype(str)
    df["position"] = df["position"].astype(int)
    df["Drug"] = df["Drug"].astype(str)
    df["Variant"] = df["Variant"].astype(str)

    final_rows = []

    for _, group in df.groupby(["Drug","Variant","position","alt"], dropna=False):

        # Check if any caller passed filtering
        pass_rows = group[group["Filter_Status"] == "PASS"]
        if not pass_rows.empty:
            search_space = pass_rows
        else:
            search_space = group

        # Apply caller priority
        selected = None
        for caller in CALLER_PRIORITY:
            rows = search_space[search_space["caller"] == caller]
            if not rows.empty:
                selected = rows.iloc[0]
                break

        # Safety fallback (should rarely happen)
        if selected is None:
            selected = search_space.head(1).iloc[0]

        final_rows.append(selected)

    return pd.DataFrame(final_rows)


def resolve_complex_variants(df):

    # COMPLEX VARIANT RESOLUTION (MNP / INDEL / SNP)
    to_remove = set()

    complex_rows = df[
        (df["ref"].astype(str).str.len() > 1) |
        (df["alt"].astype(str).str.len() > 1)
    ]

    for idx, row in complex_rows.iterrows():

        pos = int(row["position"])
        ref = str(row["ref"])
        alt = str(row["alt"])

        ref_len = len(ref)
        alt_len = len(alt)

        # CASE 1: MNP (same length substitution)
        # example: ATC -> GGG
        if ref_len == alt_len:

            mnp_zyg = str(row["zygosity"]).upper()

            for i in range(ref_len):

                mnp_pos = pos + i
                mnp_ref = ref[i]
                mnp_alt = alt[i]

                snps = df[
                    (df["position"] == mnp_pos) &
                    (df["ref"].astype(str).str.len() == 1) &
                    (df["alt"].astype(str).str.len() == 1)
                ]

                for snp_idx, snp in snps.iterrows():

                    snp_zyg = str(snp["zygosity"]).upper()

                    # SNP identical to decomposed MNP → remove
                    if snp["ref"] == mnp_ref and snp["alt"] == mnp_alt:
                        to_remove.add(snp_idx)

                    # SNP mismatch
                    else:

                        # remove only if both are HOM
                        if mnp_zyg == "HOM" and snp_zyg == "HOM":
                            to_remove.add(snp_idx)

        # CASE 2: DELETION
        # example: ATC -> A
        # deleted bases occur after anchor
        elif ref_len > alt_len:

            del_zyg = str(row["zygosity"]).upper()

            start = pos + alt_len
            end = pos + ref_len - 1

            snps = df[
                (df["position"] >= start) &
                (df["position"] <= end) &
                (df["ref"].str.len() == 1) &
                (df["alt"].str.len() == 1)
            ]

            for snp_idx, snp in snps.iterrows():

                snp_zyg = str(snp["zygosity"]).upper()

                # remove only if both are HOM
                if del_zyg == "HOM" and snp_zyg == "HOM":
                    to_remove.add(snp_idx)

        # CASE 3: INSERTION
        # example: A -> ATC
        # insertion does not invalidate SNPs
        elif alt_len > ref_len:

            continue

    final_idx = set(df.index) - to_remove

    return df.loc[sorted(final_idx)].copy()


def resolve_variant(df):

    # VARIANT RESOLUTION
    aa = df["aa_change"]
    master = df["master_change"]
    nt = df["nt_change"]

    # Extract genes
    gene_master = master.str.split("_", n=1).str[0]
    gene_aa = aa.str.split("_", n=1).str[0]
    gene_nt = nt.str.split("_", n=1).str[0]

    gene_master = gene_master.where(master != "NA")
    gene_aa = gene_aa.where(aa != "NA")
    gene_nt = gene_nt.where(nt != "NA")

    # Detect gene mismatch with catalog
    mismatch = (
        (gene_aa.notna() & (gene_aa != gene_master)) |
        (gene_nt.notna() & (gene_nt != gene_master))
    )

    # Apply variant priority
    # aa_change > master_change > nt_change
    variant = aa.copy()

    mask = variant == "NA"
    variant[mask] = master[mask]

    mask = variant == "NA"
    variant[mask] = nt[mask]

    # Force catalog variant if gene mismatch
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

    # ===================== RESOLVE MNPs-INDELs-SNPs overlay =====================
    df = resolve_complex_variants(df)

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

    # ordenar para garantir prioridade determinística
    priority_map = {c:i for i,c in enumerate(CALLER_PRIORITY)}

    df = df.sort_values(
        by=["Drug","Variant","position","alt","caller"],
        key=lambda col: col.map(priority_map) if col.name=="caller" else col
    )
    
    final = deduplicate_variants(df)[FINAL_COLUMNS]

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


