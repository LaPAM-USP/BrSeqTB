#!/usr/bin/env python3

# ============================================================
# OMS TB-DR Catalogue Processor
# Usage: python3 omsCatalog.py
# Input:  omsCatalog/WHO-UCN-TB-2023.7-eng.xlsx
# Output: 
#   - tbdr.bed  (non-redundant genomic intervals)
#   - tbdr_genomic_coordinates.csv
#   - tbdr_catalogue_master_file.csv
#   - tbdrR.csv (POS (R and r only) + resistance: variant_drug1_drug2; ...)
# Reference: OMS TB Drug Resistance Catalogue (v2023.7)
# ============================================================


import pandas as pd
import numpy as np
import os
import sys
import gc

# ============================================================
# FIXED CONFIGURATION
# ============================================================

CATALOG_DIR = "database/omsCatalog"
EXCEL_NAME = "WHO-UCN-TB-2023.7-eng.xlsx"
EXCEL_PATH = os.path.join(CATALOG_DIR, EXCEL_NAME)

REFERENCE_NAME = "NC_000962.3"

BED_OUT = os.path.join(CATALOG_DIR, "tbdr.bed")
GENOMIC_CSV_OUT = os.path.join(CATALOG_DIR, "tbdr_genomic_coordinates.csv")
MASTER_CSV_OUT = os.path.join(CATALOG_DIR, "tbdr_catalogue_master_file.csv")

R_POS_CSV = os.path.join(CATALOG_DIR, "tbdrR.csv")
INVALID_LOG = os.path.join(CATALOG_DIR, "invalid_rows.log")


# ============================================================
# FUNCTION: Generate non-redundant genomic positions (BED)
# ============================================================

def gen_nr_positions(genomic_coordinates, output_dir):

    invalid_rows = []

    def safe_int(x):
        try:
            return int(str(x).strip())
        except:
            return None

    genomic_coordinates["position_clean"] = genomic_coordinates["position"].apply(safe_int)
    invalid_mask = genomic_coordinates["position_clean"].isna()

    if invalid_mask.any():
        invalid_subset = genomic_coordinates[invalid_mask]
        invalid_rows.extend(invalid_subset.to_dict(orient="records"))
        genomic_coordinates = genomic_coordinates[~invalid_mask]

    genomic_coordinates["position"] = genomic_coordinates["position_clean"]
    genomic_coordinates.drop(columns=["position_clean"], inplace=True)

    if invalid_rows:
        with open(INVALID_LOG, "w") as log:
            log.write("Invalid rows (non-numeric POS):\n")
            for row in invalid_rows:
                log.write(str(row) + "\n")
        print(f"[WARN] Invalid rows saved to {INVALID_LOG}")

    # Identify MNPs
    ref_len = genomic_coordinates["reference_nucleotide"].str.len()
    alt_len = genomic_coordinates["alternative_nucleotide"].str.len()

    # MNP positions
    mnps_mask = (ref_len == alt_len) & (ref_len > 1)
    mnps = genomic_coordinates[mnps_mask]

    if not mnps.empty:
        starts = mnps["position"].to_numpy()
        lengths = mnps["reference_nucleotide"].str.len().to_numpy()
        offsets = np.arange(lengths.max())
        mask = offsets < lengths[:, None]
        expanded = (starts[:, None] + offsets)[mask]
        mnps_positions = pd.Series(expanded, dtype=int)
    else:
        mnps_positions = pd.Series([], dtype=int)

    # SNPs + INDELs
    snp_mask = (ref_len != alt_len) | ((ref_len == alt_len) & (ref_len == 1))
    snp_positions = genomic_coordinates.loc[snp_mask, "position"].astype(int)

    # Merge
    positions = pd.concat([mnps_positions, snp_positions], ignore_index=True)
    positions = positions.drop_duplicates().sort_values()

    if positions.empty:
        print("[WARN] No positions to write BED.")
        return

    # Build BED intervals
    bed_intervals = []
    start = positions.iloc[0]
    end = start

    for pos in positions.iloc[1:]:
        if pos == end + 1:
            end = pos
        else:
            bed_intervals.append([REFERENCE_NAME, start, end + 1])
            start = pos
            end = pos

    bed_intervals.append([REFERENCE_NAME, start, end + 1])

    df_bed = pd.DataFrame(bed_intervals, columns=["chrom", "start", "end"])
    df_bed.to_csv(BED_OUT, sep="\t", header=False, index=False)

    print(f"[OUT] BED saved: {BED_OUT}")


# ============================================================
# FUNCTION: Generate tbdrR.csv (ONE COLUMN: POS)
# ============================================================

def gen_resistance_csv(genomic_coordinates, master_df):

    print("[RUN] Generating tbdrR.csv (POS (R and r) + resistance annotation)")

    # Identify resistance-associated variants
    fcg = master_df["FINAL CONFIDENCE GRADING"].astype(str)
    mask1 = fcg.str.contains("1) Assoc w R", regex=False)
    mask2 = fcg.str.contains("2) Assoc w R - Interim", regex=False)

    resistance_df = master_df.loc[mask1 | mask2].copy()

    if resistance_df.empty:
        print("[WARN] No resistance variants found.")
        return

    # We need variant→drugs mapping
    # one variant may have multiple drugs
    resistance_df["drug"] = resistance_df["drug"].astype(str).str.strip()
    resistance_df["variant"] = resistance_df["variant"].astype(str).str.strip()

    # Build map: variant → list of drugs
    var_to_drugs = (
        resistance_df.groupby("variant")["drug"]
        .apply(lambda x: sorted(set(x)))
        .to_dict()
    )

    # Keep only genomic rows that match resistance variants
    genomic_coordinates["variant"] = genomic_coordinates["variant"].astype(str).str.strip()
    df = genomic_coordinates[genomic_coordinates["variant"].isin(var_to_drugs.keys())].copy()

    if df.empty:
        print("[WARN] No positions mapped to resistance variants.")
        return

    # Clean POS
    df["position"] = df["position"].astype(str).str.strip()
    df = df[df["position"].str.isdigit()]
    df["position"] = df["position"].astype(int)

    # For each POS: collect all variants and their drugs
    pos_map = {}

    grouped = df.groupby("position")

    for pos, subset in grouped:
        entries = []
        for variant in subset["variant"].unique():
            drugs = var_to_drugs.get(variant, [])

            # variant_drug1_drug2
            if drugs:
                entry = variant + "_" + "_".join(drugs)
            else:
                entry = variant  # no drug? (unlikely, but safe)
            entries.append(entry)

        # join variants separated by ";"
        pos_map[pos] = "; ".join(entries)

    # Build final df
    out_df = (
        pd.DataFrame({"POS": list(pos_map.keys()), "resistance": list(pos_map.values())})
        .sort_values("POS")
    )

    out_df.to_csv(R_POS_CSV, index=False)
    print(f"[OUT] Resistance positions with annotations CSV: {R_POS_CSV}")


# ============================================================
# MAIN PIPELINE ENTRYPOINT
# ============================================================

def main():

    print("[RUN] OMS Catalog processing started")

    # --------------------------------------------------------
    # SKIP IF *ALL* OUTPUT FILES ALREADY EXIST
    # --------------------------------------------------------
    required_outputs = [
        BED_OUT,
        GENOMIC_CSV_OUT,
        MASTER_CSV_OUT,
        R_POS_CSV
    ]

    if all(os.path.isfile(f) for f in required_outputs):
        print("[SKIP] All catalog output files already exist. Nothing to do.")
        print("       Outputs:")
        for f in required_outputs:
            print(f"       - {f}")
        return

    # --------------------------------------------------------
    # CHECK INPUTS
    # --------------------------------------------------------
    if not os.path.isdir(CATALOG_DIR):
        print(f"[ERROR] Directory '{CATALOG_DIR}' missing.")
        sys.exit(1)

    if not os.path.isfile(EXCEL_PATH):
        print(f"[ERROR] Excel file missing: {EXCEL_PATH}")
        sys.exit(1)

    print("---------------------------------------------")
    print(f"[IN] Loading Excel: {EXCEL_PATH}")

    # Load genomic coordinates
    genomic_coordinates = pd.read_excel(EXCEL_PATH, sheet_name="Genomic_coordinates")
    genomic_coordinates.to_csv(GENOMIC_CSV_OUT, index=False)
    print(f"[OUT] Saved: {GENOMIC_CSV_OUT}")

    # Load master file
    master_file = pd.read_excel(EXCEL_PATH, sheet_name="Catalogue_master_file", header=2)
    keep_cols = ["drug","variant","gene","tier","effect","FINAL CONFIDENCE GRADING","Comment"]
    master_file = master_file[[c for c in keep_cols if c in master_file.columns]]
    master_file.to_csv(MASTER_CSV_OUT, index=False)
    print(f"[OUT] Saved: {MASTER_CSV_OUT}")

    # Generate BED (unchanged behavior)
    gen_nr_positions(genomic_coordinates.copy(), CATALOG_DIR)

    # Generate tbdrR.csv (NEW)
    gen_resistance_csv(genomic_coordinates.copy(), master_file)

    print("---------------------------------------------")
    print("[OK] Finished OMS Catalog processing.")


if __name__ == "__main__":
    main()

