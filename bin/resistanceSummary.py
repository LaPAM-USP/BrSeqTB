#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ============================================================
# Summarize TB resistance profiles into a single master table
#
# NEW INPUT:
#   results/resistance/<biosample>.xlsx
#
# Output:
#   results/resistance_summary.xlsx
#
# Rules:
#   - Only PASS variants
#   - Read curated variants (already deduplicated/MNP-corrected)
#   - annotation format:  Variant_flagX
#   - phenotype R/r/u/s/S
#   - drug columns interleaved with annotation columns
# ============================================================

import os
import pandas as pd

# ============================================================
# PROJECT DIR (script always lives in bin/)
# ============================================================

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# ============================================================
# CONFIG
# ============================================================

INPUT_DIR = os.path.join(PROJECT_DIR, "results", "resistance")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "results")
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "resistance_summary.xlsx")

# ============================================================
# PATTERNS — convert Evidence to flags
# ============================================================

patterns = {
    "R": "flagR",
    "r": "flagr",
    "u": "flagu",
    "s": "flagnr",
    "S": "flagnR"
}

def determine_final_resistance(string):
    """Return phenotype based on annotation flags."""
    if "flagR" in string:
        return "R"
    elif "flagr" in string:
        return "r"
    elif "flagu" in string:
        return "u"
    elif "flagnR" in string:
        return "S"
    elif "flagnr" in string:
        return "s"
    return "S"


# ============================================================
# PROCESS ONE BIOSAMPLE
# ============================================================

def process_biosample(biosample):

    target_file = os.path.join(INPUT_DIR, f"{biosample}.xlsx")
    if not os.path.exists(target_file):
        return None

    df = pd.read_excel(target_file)

    # --------------------------------------------------------
    # Keep only PASS variants
    # --------------------------------------------------------
    df = df[df["Filter_Status"] == "PASS"]
    if df.empty:
        return {"biosample": biosample}

    # Drugs appear already as normalized strings in CAPITAL/LOWER policy
    drugs = sorted(df["Drug"].dropna().unique())

    result = {"biosample": biosample}

    # ========================================================
    # Build phenotype and annotation per drug
    # ========================================================
    for drug in drugs:

        subset = df[df["Drug"] == drug]
        if subset.empty:
            continue

        annotations = []

        for _, row in subset.iterrows():

            # row["Evidence"] already contains: R/r/u/s/S
            phenotype_code = str(row["Evidence"]).strip()
            flag_code = patterns.get(phenotype_code, "")

            variant = row["Variant"]
            if pd.isna(variant):
                continue
            variant = str(variant).strip()

            annot = f"{variant}_{flag_code}"
            annotations.append(annot)

        # Remove duplicates
        unique_ann = sorted(set(annotations))
        joined_ann = ",".join(unique_ann)

        # Determine final phenotype for the drug
        phenotype = determine_final_resistance(joined_ann)

        drug_col = drug.lower()
        result[drug_col] = phenotype
        result[drug_col + "_annotation"] = joined_ann

    return result


# ============================================================
# MAIN
# ============================================================

def main():

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    biosamples = [
        os.path.splitext(f)[0]
        for f in os.listdir(INPUT_DIR)
        if f.endswith(".xlsx")
    ]

    all_rows = []

    for biosample in biosamples:
        data = process_biosample(biosample)
        if data:
            all_rows.append(data)

    if not all_rows:
        print("No biosample results found.")
        return

    df = pd.DataFrame(all_rows).fillna("")

    # ============================================================
    # ORDER COLUMNS — biosample + interleaved drug/annotation
    # ============================================================

    fixed_cols = ["biosample"]

    drug_cols = sorted([
        c for c in df.columns
        if c not in fixed_cols and not c.endswith("_annotation")
    ])

    ordered_cols = fixed_cols[:]
    for d in drug_cols:
        ordered_cols.append(d)
        ann = d + "_annotation"
        if ann in df.columns:
            ordered_cols.append(ann)

    df = df.reindex(columns=ordered_cols)
    
    # Fill empty phenotype cells with "S", but do not touch annotation columns
    for col in df.columns:
        if col != "biosample" and not col.endswith("_annotation"):
            df[col] = df[col].replace("", "S")


    df.to_excel(OUTPUT_FILE, index=False)
    print("Summary written to:", OUTPUT_FILE)


if __name__ == "__main__":
    main()

