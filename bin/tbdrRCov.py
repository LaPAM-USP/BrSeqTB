#!/usr/bin/env python3
# ============================================================
# TBDR-R Coverage Checker (interval expansion for SNP/INDEL/MNP)
# Usage: python3 tbdrRCov.py <biosample>
#
# Inputs:
#   database/omsCatalog/tbdrR.csv   (columns: POS, resistance)
#   gatk/<biosample>/<biosample>.g.vcf.gz
#
# Output:
#   tbdrRCov/<biosample>/<biosample>_tbdrRcov_summary.csv
# ============================================================

import os
import sys
import pandas as pd
import pysam
import numpy as np

PROJECT_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..")
)

CATALOG_DIR = os.path.join(PROJECT_DIR, "database", "omsCatalog")
TBDR_CSV = os.path.join(CATALOG_DIR, "tbdrR.csv")

GATK_DIR = os.path.join(PROJECT_DIR, "gatk")
OUT_DIR  = os.path.join(PROJECT_DIR, "tbdrRCov")

REFERENCE_NAME = "NC_000962.3"


# ============================================================
# LOAD tbdrR.csv (POS + resistance)
# ============================================================
def load_resistance_positions(path):
    df = pd.read_csv(path)

    if "POS" not in df.columns or "resistance" not in df.columns:
        raise ValueError("tbdrR.csv must contain POS and resistance columns.")

    df["POS"] = df["POS"].astype(int)
    df = df.sort_values("POS")

    return df


# ============================================================
# EXTRACT INTERVALS FROM GVCF
# ============================================================
def load_gvcf_intervals(gvcf_path):
    print("[LOAD] Extracting intervals from gVCF...")

    vcf = pysam.VariantFile(gvcf_path)

    starts = []
    ends = []
    dps = []

    import re
    pattern_end = re.compile(r"END=(\d+)")

    for rec in vcf.fetch(REFERENCE_NAME):
        sample = rec.samples[0]

        dp = sample.get("DP")
        dp = int(dp) if dp is not None else 0

        real_alts = [a for a in (rec.alts or []) if a != "<NON_REF>"]

        raw = str(rec).strip()

        # --- CASE 1: reference block ---
        if not real_alts:
            end_val = rec.info.get("END", None)

            if end_val is None:
                m = pattern_end.search(raw)
                end_val = int(m.group(1)) if m else rec.pos

            start = int(rec.pos)
            end = int(end_val)

        # --- CASE 2: real variant ---
        else:
            ref = rec.ref
            alt = real_alts[0]
            start = int(rec.pos)

            ref_len = len(ref)
            alt_len = len(alt)

            if ref_len > alt_len:
                end = start + ref_len - 1       # deletion / MNP
            elif ref_len == 1 and alt_len == 1:
                end = start                    # SNP
            else:
                end = start                    # insertion

        starts.append(start)
        ends.append(end)
        dps.append(dp)

    print(f"[INFO] Loaded {len(starts)} intervals.")
    return np.array(starts), np.array(ends), np.array(dps)


# ============================================================
# BINARY SEARCH FOR INTERVAL MATCHING
# ============================================================
def get_dp_for_position(pos, starts, ends, dps):
    idx = np.searchsorted(starts, pos, side="right") - 1
    if idx < 0:
        return 0
    if starts[idx] <= pos <= ends[idx]:
        return int(dps[idx])
    return 0


# ============================================================
# MAIN
# ============================================================
def main():
    if len(sys.argv) != 2:
        print("Usage: python3 tbdrRCov.py <biosample>")
        sys.exit(1)

    biosample = sys.argv[1]

    print(f"[RUN] Checking coverage for biosample: {biosample}")

    gvcf_path = os.path.join(GATK_DIR, biosample, f"{biosample}.g.vcf.gz")

    out_folder = os.path.join(OUT_DIR, biosample)
    os.makedirs(out_folder, exist_ok=True)

    out_file = os.path.join(out_folder, f"{biosample}_tbdrRcov_summary.csv")

    # ===================== SKIP IF ALREADY DONE =====================
    if os.path.exists(out_file):
        print(f"[SKIP] Coverage summary already exists: {out_file}")
        sys.exit(0)

    # ===================== INPUT CHECKS =====================
    if not os.path.exists(TBDR_CSV):
        print(f"[ERROR] tbdrR.csv not found: {TBDR_CSV}")
        sys.exit(1)

    if not os.path.exists(gvcf_path):
        print(f"[ERROR] gVCF not found: {gvcf_path}")
        sys.exit(1)

    # ===================== LOAD DATA =====================
    df_tbdr = load_resistance_positions(TBDR_CSV)
    print(f"[INFO] Resistance positions loaded: {len(df_tbdr)}")

    starts, ends, dps = load_gvcf_intervals(gvcf_path)

    order = np.argsort(starts)
    starts = starts[order]
    ends = ends[order]
    dps = dps[order]

    # ===================== ANALYSIS =====================
    low_cov = []

    for _, row in df_tbdr.iterrows():
        pos = int(row["POS"])
        variant = row["resistance"]

        dp = get_dp_for_position(pos, starts, ends, dps)

        if dp < 10:
            low_cov.append([REFERENCE_NAME, pos, dp, variant])

    # ===================== OUTPUT =====================
    df = pd.DataFrame(low_cov, columns=["CHROM", "POS", "DP", "VARIANT"])
    df.to_csv(out_file, index=False)

    print(f"[OUT] Saved: {out_file}")
    print("[DONE]")


if __name__ == "__main__":
    main()

