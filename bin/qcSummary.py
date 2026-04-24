#!/usr/bin/env python3

# ============================================================
# Master QC Summary Generator for TB Pipeline
#
# Usage: python3 qcSummary.py
#
# [IN]
#   fastqc/<biosample>/<biosample>_fastqc_summary.csv
#   trimmomatic/<biosample>/<biosample>_trimmomatic_summary.csv
#   bwa/<biosample>/<biosample>_bwa_summary.csv
#   ntmFilter/<biosample>/<biosample>_ntm_summary.csv
#   mixInfection/<biosample>/<biosample>_mixinfection_summary.csv
#   kaiju/<biosample>/<biosample>_kaiju_summary.csv
#   lineage/<biosample>/<biosample>_lineage_summary.csv
#   tbdrRCov/<biosample>/<biosample>_tbdrRcov_summary.csv   ← NOVO
#
# [OUT]
#   results/qc_summary.xlsx
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

OUTPUT_DIR = os.path.join(PROJECT_DIR, "results")
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "qc_summary.xlsx")

INPUT_DIRS = {
    "fastqc":       os.path.join(PROJECT_DIR, "fastqc"),
    "trimmomatic":  os.path.join(PROJECT_DIR, "trimmomatic"),
    "bwa":          os.path.join(PROJECT_DIR, "bwa"),
    "ntmFilter":    os.path.join(PROJECT_DIR, "ntmFilter"),
    "mixInfection": os.path.join(PROJECT_DIR, "mixInfection"),
    "kaiju":        os.path.join(PROJECT_DIR, "kaiju"),
    "lineage":      os.path.join(PROJECT_DIR, "lineage"),
    "tbdrRCov":     os.path.join(PROJECT_DIR, "tbdrRCov"),
}

SUMMARY_SUFFIX = {
    "fastqc":      "_fastqc_summary.csv",
    "trimmomatic": "_trimmomatic_summary.csv",
    "bwa":         "_bwa_summary.csv",
    "ntmFilter":   "_ntm_summary.csv",
    "mixInfection":"_mixinfection_summary.csv",
    "kaiju":       "_kaiju_summary.csv",
    "lineage":     "_lineage_summary.csv",
    "tbdrRCov":    "_tbdrRcov_summary.csv",
}

KAIJU_IGNORE = {
    "cannot be assigned to a (non-viral) species",
    "unclassified"
}

# ============================================================
# HELPERS
# ============================================================

def read_summary_file(base_dir, biosample, suffix):
    path = os.path.join(base_dir, biosample, f"{biosample}{suffix}")
    if not os.path.exists(path):
        return None
    try:
        return pd.read_csv(path)
    except Exception:
        return None


def process_kaiju(df):
    if df is None or df.empty:
        return None

    df_valid = df.copy()
    df_valid["reads"] = pd.to_numeric(df_valid["reads"], errors="coerce").fillna(0)

    df_main = df_valid[~df_valid["taxon_name"].isin(KAIJU_IGNORE)]

    if not df_main.empty:
        return df_main.loc[df_main["reads"].idxmax()]

    return df_valid.loc[df_valid["reads"].idxmax()]


def process_lineage(df, biosample):
    if df is None or df.empty:
        return {"biosample": biosample, "LINEAGE": "", "LINEAGE_DETAILS": "", "COMMENT": ""}

    lineage = ";".join(sorted(df["LINEAGE"].dropna().astype(str).unique()))
    details = ";".join(sorted(df["LINEAGE_DETAILS"].dropna().astype(str).unique()))
    comment = ";".join(sorted(df["COMMENT"].dropna().astype(str).unique()))

    return {
        "biosample": biosample,
        "LINEAGE": lineage,
        "LINEAGE_DETAILS": details,
        "COMMENT": comment
    }


# ============================================================
# process tbdrRCov
# ============================================================

def process_tbdr_rcov(df, biosample):
    """Return biosample + all VARIANT entries concatenated for dp<10."""
    if df is None or df.empty:
        return {"biosample": biosample, "Cov < 10": ""}

    variants = []
    for v in df["VARIANT"].dropna().astype(str):
        for item in v.split(";"):
            item = item.strip()
            if item not in variants:
                variants.append(item)

    return {
        "biosample": biosample,
        "Cov < 10": ";".join(variants)
    }


# ============================================================
# MAIN
# ============================================================

def main():

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # Determine biosamples from fastqc directory
    fastqc_dir = INPUT_DIRS["fastqc"]
    biosamples = [
        d for d in os.listdir(fastqc_dir)
        if os.path.isdir(os.path.join(fastqc_dir, d))
    ]

    tables = {
        "fastqc": [],
        "trimmomatic": [],
        "bwa": [],
        "ntmFilter": [],
        "mixInfection": [],
        "kaiju": [],
        "lineage": [],
        "tbdrRCov": [],  # ← NEW
    }

    # --------------------------------------------------------
    # Per-biosample iteration
    # --------------------------------------------------------
    for biosample in biosamples:

        # direct concatenation
        for analysis in ["fastqc","trimmomatic","bwa","ntmFilter","mixInfection"]:
            df = read_summary_file(INPUT_DIRS[analysis], biosample, SUMMARY_SUFFIX[analysis])
            if df is not None and not df.empty:
                tables[analysis].append(df)

        # Kaiju
        df_k = read_summary_file(INPUT_DIRS["kaiju"], biosample, SUMMARY_SUFFIX["kaiju"])
        top = process_kaiju(df_k)
        if top is not None:
            tables["kaiju"].append(top.to_frame().T)

        # Lineage
        df_l = read_summary_file(INPUT_DIRS["lineage"], biosample, SUMMARY_SUFFIX["lineage"])
        lin = process_lineage(df_l, biosample)
        tables["lineage"].append(pd.DataFrame([lin]))

        # tbdrRCov  ← NEW
        tbdr_path = os.path.join(INPUT_DIRS["tbdrRCov"], biosample, f"{biosample}_tbdrRcov_summary.csv")
        if os.path.exists(tbdr_path):
            df_tb = pd.read_csv(tbdr_path)
        else:
            df_tb = None

        tbdr_info = process_tbdr_rcov(df_tb, biosample)
        tables["tbdrRCov"].append(pd.DataFrame([tbdr_info]))

    # --------------------------------------------------------
    # Write Excel
    # --------------------------------------------------------
    with pd.ExcelWriter(OUTPUT_FILE, engine="xlsxwriter") as writer:

        # direct sheets
        for analysis in ["fastqc","trimmomatic","bwa","ntmFilter","mixInfection"]:
            if tables[analysis]:
                pd.concat(tables[analysis], ignore_index=True).to_excel(
                    writer, sheet_name=analysis, index=False
                )
            else:
                pd.DataFrame().to_excel(writer, sheet_name=analysis, index=False)

        # Kaiju
        if tables["kaiju"]:
            pd.concat(tables["kaiju"], ignore_index=True).to_excel(
                writer, sheet_name="kaiju", index=False
            )
        else:
            pd.DataFrame().to_excel(writer, sheet_name="kaiju", index=False)

        # Lineage
        if tables["lineage"]:
            pd.concat(tables["lineage"], ignore_index=True).to_excel(
                writer, sheet_name="lineage", index=False
            )
        else:
            pd.DataFrame().to_excel(writer, sheet_name="lineage", index=False)

        # tbdrRCov (NEW)
        if tables["tbdrRCov"]:
            pd.concat(tables["tbdrRCov"], ignore_index=True).to_excel(
                writer, sheet_name="tbdrRCov", index=False
            )
        else:
            pd.DataFrame().to_excel(writer, sheet_name="tbdrRCov", index=False)

    print("QC summary written to:", OUTPUT_FILE)


if __name__ == "__main__":
    main()

