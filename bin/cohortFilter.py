#!/usr/bin/env python3
# ============================================================
# Cohort Variant Summary – Excel Generator (Pipeline-safe)
#
# Behavior:
#   - Detects project root automatically
#   - Always reads from: <PROJECT_DIR>/cohort/
#
# Input:
#   cohort/cohort_snps_filtered.vcf.gz
#   cohort/cohort_indels_filtered.vcf.gz
#
# Output:
#   cohort/filter/filter.xlsx
# ============================================================

import gzip
import sys
import pandas as pd
from pathlib import Path
import os

# ============================================================
# PROJECT DIR (same pattern as tbdrRCov)
# ============================================================
PROJECT_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..")
)

COHORT_DIR = Path(PROJECT_DIR) / "cohort"

SNPS_VCF   = COHORT_DIR / "cohort_snps_filtered.vcf.gz"
INDELS_VCF = COHORT_DIR / "cohort_indels_filtered.vcf.gz"

FILTER_DIR = COHORT_DIR / "filter"
OUTPUT_XLSX = FILTER_DIR / "filter.xlsx"


# ============================================================
# VCF READER
# ============================================================
def read_vcf(vcf_path: Path):
    records = []
    with gzip.open(vcf_path, "rt") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 8:
                records.append(parts[:8])
    return records


# ============================================================
# MAIN
# ============================================================
def main():

    if not COHORT_DIR.exists():
        print(f"[ERROR] cohort directory not found: {COHORT_DIR}", file=sys.stderr)
        sys.exit(1)

    if not SNPS_VCF.exists():
        print(f"[ERROR] Missing SNP VCF: {SNPS_VCF}", file=sys.stderr)
        sys.exit(1)

    if not INDELS_VCF.exists():
        print(f"[ERROR] Missing INDEL VCF: {INDELS_VCF}", file=sys.stderr)
        sys.exit(1)

    FILTER_DIR.mkdir(exist_ok=True)

    rows = []
    rows.extend(read_vcf(SNPS_VCF))
    rows.extend(read_vcf(INDELS_VCF))

    if not rows:
        print("[WARN] No variants found in cohort VCFs.")
        return

    df = pd.DataFrame(
        rows,
        columns=["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    )

    df.to_excel(OUTPUT_XLSX, index=False)

    print(f"[OK] Excel created: {OUTPUT_XLSX}")


if __name__ == "__main__":
    main()

