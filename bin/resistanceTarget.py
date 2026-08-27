#!/usr/bin/env python3

# ============================================================
# TB Resistance Profiling from VCFs (GATK / norm / LoFreq / Delly)
# Usage: python3 resistanceTarget.py <biosample ID>
# Input:  snpeff/<biosample ID>/<biosample ID>_{gatk,norm,lofreq,delly}.vcf.gz
# Output:
#   resistance/<biosample ID>/*_{caller}_ANN.xlsx
#   resistance/<biosample ID>/<biosample ID>_OMStarget.xlsx
# Reference:
#   database/omsCatalog/tbdr_genomic_coordinates.csv
#   database/omsCatalog/tbdr_catalogue_master_file.csv
# ============================================================


import pandas as pd
import pysam
import os
import sys

# ============================================================
# PROJECT DIR (script always lives in bin/)
# ============================================================

PROJECT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# ============================================================
# CONSTANTS / PATHS
# ============================================================

REFERENCE_NAME = "NC_000962.3"

CATALOG_DIR = os.path.join(PROJECT_DIR, "database", "omsCatalog")
CATALOG_GENOMIC_COORDS = os.path.join(CATALOG_DIR, "tbdr_genomic_coordinates.csv")
CATALOG_MASTER = os.path.join(CATALOG_DIR, "tbdr_catalogue_master_file.csv")

COHORT_FILTER = os.path.join(PROJECT_DIR, "cohort", "filter", "filter.xlsx")

SNPEFF_BASE = os.path.join(PROJECT_DIR, "snpeff")
RESULTS_BASE = os.path.join(PROJECT_DIR, "resistance")

VCF_CALLERS = {
    "gatk": "_gatk.vcf.gz",
    "norm": "_norm.vcf.gz",
    "lofreq": "_lofreq.vcf.gz",
    "delly": "_delly.vcf.gz"
}

# ============================================================
# FILTER HELPERS
# ============================================================

def is_missing(value):
    return value is None or pd.isna(value)

def get_first_af(value):
    if is_missing(value):
        return None

    if isinstance(value, (list, tuple)):
        if len(value) == 0:
            return None
        return value[0]

    return value

def is_pass_filter(filter_vcf):
    return filter_vcf == "PASS"


def get_record_filter(record):
    filters = list(record.filter.keys())
    if not filters:
        return "."
    return ";".join(filters)


def get_record_qual(record):
    if record.qual is None:
        return "."
    return record.qual


def get_numeric_qual(qual):
    if qual == "." or is_missing(qual):
        return 0
    return qual


def get_cohort_filter(filter_df, pos):
    row = filter_df[filter_df["POS"] == pos]
    if row.empty:
        return "NO_COHORT_SITE"
    return row.iloc[0]["FILTER"]


def gatk_norm_failures(filter_vcf, af, alt_reads, total_depth):

    failures = []

    # VCF hard-filtering level
    if not is_pass_filter(filter_vcf):
        failures.append(filter_vcf)

    # Internal per-sample / allele-level filters
    dp = total_depth if not is_missing(total_depth) else 0

    if is_missing(alt_reads) or alt_reads < 3:
        failures.append("ALT_READS_BELOW_3")

    if is_missing(af) or af < 0.00:
        failures.append("AF_FAIL")

    if dp < 10:
        failures.append("DP_FAIL")

    return failures


def lofreq_failures(filter_vcf, alt_reads, total_depth):

    failures = []

    if not is_pass_filter(filter_vcf):
        failures.append(filter_vcf)

    if is_missing(alt_reads) or alt_reads < 3:
        failures.append("ALT_READS_BELOW_3")

    dp = total_depth if not is_missing(total_depth) else 0
    if dp < 10:
        failures.append("DP_FAIL")

    return failures


def delly_failures(filter_vcf):

    failures = []

    if not is_pass_filter(filter_vcf):
        failures.append(filter_vcf)

    return failures


def define_filter_status(
    caller,
    position,
    af,
    alt_reads,
    total_depth,
    filter_vcf,
    filter_df=None
):

    # ========================================================
    # GATK / NORM
    # First checks cohort filter, then per-sample filters
    # ========================================================
    if caller in ("gatk", "norm"):

        cohort_val = get_cohort_filter(filter_df, position)

        if cohort_val == "PASS":
            return "PASS", "cohort"

        per_sample_fails = gatk_norm_failures(
            filter_vcf=filter_vcf,
            af=af,
            alt_reads=alt_reads,
            total_depth=total_depth,
        )

        if per_sample_fails:
            combined = [cohort_val] + per_sample_fails
            return ";".join(combined), "cohort+per_sample"

        return "PASS", "per_sample"

    # ========================================================
    # LOFREQ
    # Needs PASS in VCF FILTER and at least 3 ALT reads
    # ========================================================
    if caller == "lofreq":

        per_sample_fails = lofreq_failures(
            filter_vcf=filter_vcf,
            alt_reads=alt_reads,
            total_depth=total_depth
        )

        if per_sample_fails:
            return ";".join(per_sample_fails), "per_sample"

        return "PASS", "per_sample"

    # ========================================================
    # DELLY
    # Needs PASS in VCF FILTER
    # ========================================================
    if caller == "delly":

        per_sample_fails = delly_failures(
            filter_vcf=filter_vcf
        )

        if per_sample_fails:
            return ";".join(per_sample_fails), "per_sample"

        return "PASS", "per_sample"

    return "UNKNOWN_CALLER", "filter_error"


# ============================================================
# ANNOTATION
# ============================================================

def recover_annotation(vcf_file, caller):

    vcf = pysam.VariantFile(vcf_file)
    annotations = []

    # To recover variant metrics
    is_delly = caller == "delly"
    is_lofreq = caller == "lofreq"
    is_gatk = caller in ("gatk", "norm")

    filter_df = None

    # If gatk-norm use cohort filter
    if is_gatk:
        filter_df = pd.read_excel(COHORT_FILTER)

    for record in vcf:
        position = record.pos
        ref = record.ref
        alt_alleles = record.alts or []
        qual = get_record_qual(record)
        filter_vcf = get_record_filter(record)

        AF_list = []
        ALT_reads_list = []
        Total_depth = None
        zygosity_list = []

        # --------- DELLY ---------
        if is_delly:
            AF_list = [None] * len(alt_alleles)
            ALT_reads_list = [None] * len(alt_alleles)
            Total_depth = None
            zygosity_list = ["NA"] * len(alt_alleles)

        # --------- GATK / NORM -----
        elif is_gatk:
            sample_name = list(vcf.header.samples)[0]
            fmt = record.samples[sample_name]
            AD = fmt.get("AD", None)
            DP = fmt.get("DP", None)

            if AD and len(AD) >= 1:
                ALT_reads_list = AD[1:]
                Total_depth = sum(AD)

                AF_list = [
                    (ar / Total_depth if Total_depth > 0 else None)
                    for ar in ALT_reads_list
                ]

                zygosity_list = [
                    ("HOM" if af is not None and af >= 0.90 else "HET")
                    if af is not None else "NA"
                    for af in AF_list
                ]
            else:
                Total_depth = DP
                ALT_reads_list = [None] * len(alt_alleles)
                AF_list = [None] * len(alt_alleles)
                zygosity_list = ["NA"] * len(alt_alleles)

        # --------- LOFREQ ----------
        elif is_lofreq:
            Total_depth = record.info.get("DP", None)
            af_value = get_first_af(record.info.get("AF", None))
            AF_list = [af_value]

            dp4 = record.info.get("DP4", None)
            ALT_reads_list = [dp4[2] + dp4[3] if dp4 and len(dp4) == 4 else None]

            zygosity_list = [
                "HOM" if AF_list[0] is not None and AF_list[0] >= 0.90 else "HET"
            ]

        # ============================================================
        # ANN parsing - snpEff
        # ============================================================

        ann_list_raw = record.info.get("ANN", [])
        ann_by_alt = {alt: [] for alt in alt_alleles}

        for ann in ann_list_raw:
            fields = ann.split("|")
            ann_alt = fields[0]
            if ann_alt in ann_by_alt:
                ann_by_alt[ann_alt].append(fields)

        parsed_ann = []

        for alt in alt_alleles:
            entries = ann_by_alt.get(alt, [])

            if not entries:
                parsed_ann.append(("NA", "NA"))
                continue

            selected = None

            for f in entries:
                if f[3] or f[9] or f[10]:
                    selected = f
                    break

            if not selected:
                parsed_ann.append(("NA", "NA"))
                continue

            gene = selected[3] if selected[3] else "NA"
            genic = selected[9] if selected[9] else "NA"
            prot = selected[10] if selected[10] else "NA"

            nt = f"{gene}_{genic}" if genic != "NA" else "NA"
            aa = f"{gene}_{prot}" if prot != "NA" else "NA"

            parsed_ann.append((nt, aa))

        # ============================================================
        # OUTPUT ROWS
        # ============================================================

        for i, alt in enumerate(alt_alleles):
            nt, aa = parsed_ann[i]
            af = AF_list[i] if i < len(AF_list) else None
            ar = ALT_reads_list[i] if i < len(ALT_reads_list) else None
            zy = zygosity_list[i] if i < len(zygosity_list) else "NA"

            filter_status, filter_method = define_filter_status(
                caller=caller,
                position=position,
                af=af,
                alt_reads=ar,
                total_depth=Total_depth,
                filter_vcf=filter_vcf,
                filter_df=filter_df
            )

            annotations.append([
                position,
                ref,
                alt,
                nt,
                aa,
                zy,
                af,
                ar,
                Total_depth,
                qual,
                filter_status,
                filter_method,
                caller
            ])

    return pd.DataFrame(annotations, columns=[
        "position",
        "ref",
        "alt",
        "nt_change",
        "aa_change",
        "zygosity",
        "af",
        "alt_reads",
        "total_depth",
        "qual",
        "filter_status",
        "filter_method",
        "caller",
    ])


# ============================================================
# NORMALIZATION
# ============================================================

def normalize_ann(df):
    df["ref"] = df["ref"].str.upper().str.strip()
    df["alt"] = df["alt"].str.upper().str.strip()
    df["position"] = df["position"].astype(int)
    return df


def normalize_catalog(df):
    df["position"] = df["position"].astype(float).astype(int)
    df["reference_nucleotide"] = df["reference_nucleotide"].str.upper().str.strip()
    df["alternative_nucleotide"] = df["alternative_nucleotide"].str.upper().str.strip()
    return df


# ============================================================
# MATCHING
# ============================================================

def matching(df_ann, coord_file_path, master_file_path):
    coord_file = normalize_catalog(pd.read_csv(coord_file_path))
    master = pd.read_csv(master_file_path)
    df_ann = normalize_ann(df_ann)

    original_columns = list(df_ann.columns)

    # ------------------------------------------------------------
    # Matching method 1 - by coordinate (pos-ref-alt)
    # Tests ALL variants
    # ------------------------------------------------------------
    coord_match = df_ann.merge(
        coord_file[[
            "position",
            "reference_nucleotide",
            "alternative_nucleotide",
            "variant"
        ]],
        left_on=["position", "ref", "alt"],
        right_on=[
            "position",
            "reference_nucleotide",
            "alternative_nucleotide"
        ],
        how="left"
    ).rename(columns={"variant": "master_change"})

    coord_match = coord_match.merge(
        master,
        left_on="master_change",
        right_on="variant",
        how="left"
    )

    coord_match["match_method"] = None
    coord_match.loc[
        coord_match["drug"].notna(),
        "match_method"
    ] = "coord"

    coord_hits = coord_match[
        coord_match["match_method"] == "coord"
    ].copy()

    # ------------------------------------------------------------
    # Matching method 2 - by nt_change
    # Tests ALL variants, regardless of coordinate match
    # ------------------------------------------------------------
    nt_match = df_ann.merge(
        master,
        left_on="nt_change",
        right_on="variant",
        how="left"
    )

    nt_match["master_change"] = None
    nt_match.loc[
        nt_match["drug"].notna(),
        "master_change"
    ] = nt_match.loc[
        nt_match["drug"].notna(),
        "variant"
    ]

    nt_match["match_method"] = None
    nt_match.loc[
        nt_match["drug"].notna(),
        "match_method"
    ] = "nt_change"

    nt_hits = nt_match[
        nt_match["match_method"] == "nt_change"
    ].copy()

    # ------------------------------------------------------------
    # Matching method 3 - by aa_change
    # Tests ALL variants, regardless of coordinate or nt_change match
    # ------------------------------------------------------------
    aa_match = df_ann.merge(
        master,
        left_on="aa_change",
        right_on="variant",
        how="left"
    )

    aa_match["master_change"] = None
    aa_match.loc[
        aa_match["drug"].notna(),
        "master_change"
    ] = aa_match.loc[
        aa_match["drug"].notna(),
        "variant"
    ]

    aa_match["match_method"] = None
    aa_match.loc[
        aa_match["drug"].notna(),
        "match_method"
    ] = "aa_change"

    aa_hits = aa_match[
        aa_match["match_method"] == "aa_change"
    ].copy()

    # ------------------------------------------------------------
    # Final merging WITHOUT hierarchy
    # Keeps all matches from coord, nt_change and aa_change
    # Duplicates are intentionally preserved
    # ------------------------------------------------------------
    final = pd.concat([
        coord_hits,
        nt_hits,
        aa_hits
    ], ignore_index=True)

    final["master_change"] = final["master_change"].fillna("NA")

    # optional - delete this line code if master file is changed
    final = final.dropna(subset=[
        "drug",
        "variant",
        "tier",
        "effect",
        "FINAL CONFIDENCE GRADING"
    ])

    # Remove helper columns from coordinate merge
    final = final.drop(
        columns=[
            "reference_nucleotide",
            "alternative_nucleotide"
        ],
        errors="ignore"
    )

    # Keep current annotation/filter columns first,
    # then append OMS/matching columns at the end
    matching_columns = [
        "drug",
        "variant",
        "gene",
        "tier",
        "effect",
        "FINAL CONFIDENCE GRADING",
        "comment",
        "master_change",
        "match_method"
    ]

    final_columns = original_columns + [
        col for col in matching_columns
        if col in final.columns and col not in original_columns
    ]

    return final[final_columns]

# ============================================================
# VCF TO DF for each caller
# Saves ANN output and returns matched OMS target rows
# ============================================================

def vcf_to_df(biosample, vcf_path):

    caller = os.path.basename(vcf_path).replace(f"{biosample}_", "").replace(".vcf.gz", "")
    results_dir = os.path.join(RESULTS_BASE, biosample)
    os.makedirs(results_dir, exist_ok=True)

    ann_out = os.path.join(results_dir, f"{biosample}_{caller}_ANN.xlsx")

    print(f"[RUN] Profiling caller: {caller}")
    print(f"   VCF   : {vcf_path}")
    print(f"   ANN   : {ann_out}")
    print("---------------------------------------------")

    df_ann = recover_annotation(vcf_path, caller)
    df_ann = normalize_ann(df_ann)

    # Save full ANN table for this caller
    df_ann.to_excel(ann_out, index=False)

    # Match against OMS catalogue
    df_target = matching(
        df_ann=df_ann,
        coord_file_path=CATALOG_GENOMIC_COORDS,
        master_file_path=CATALOG_MASTER
    )

    print(f"[OK] Completed caller: {caller}")
    print(f"   OMS target rows: {len(df_target)}\n")

    return df_target


# ============================================================
# MAIN — GENERATE THE SINGLE FINAL TARGET FILE
# ============================================================

def main():

    if len(sys.argv) != 2:
        print("Usage: python3 resistanceTarget.py <biosample>")
        sys.exit(1)

    biosample = sys.argv[1]
    sample_dir = os.path.join(SNPEFF_BASE, biosample)

    # ============================================================
    # STRICT VALIDATION — ALL CALLERS REQUIRED
    # ============================================================

    missing = []

    for caller, suffix in VCF_CALLERS.items():
        vcf_path = os.path.join(sample_dir, f"{biosample}{suffix}")
        if not os.path.exists(vcf_path):
            missing.append(vcf_path)

    if missing:
        print("[ERROR] Missing required VCF files:")
        for m in missing:
            print(f"  - {m}")
        print("[ABORT] Resistance profiling requires ALL callers (gatk, norm, lofreq, delly).")
        sys.exit(1)

    # ============================================================
    # CONTINUE NORMAL EXECUTION
    # ============================================================

    print(f"[INFO] Scanning VCFs in {sample_dir}")

    merged_targets = []

    for caller, suffix in VCF_CALLERS.items():
        vcf_path = os.path.join(sample_dir, f"{biosample}{suffix}")
        if os.path.exists(vcf_path):
            df_target = vcf_to_df(biosample, vcf_path)
            merged_targets.append(df_target)
        else:
            print(f"[SKIP] {caller}: {vcf_path} not found")

    # ============================================================
    # GENERATE THE SINGLE COMBINED TARGET
    # ============================================================

    if merged_targets:

        final_df = pd.concat(merged_targets, ignore_index=True)

        # Put selected filter/QC columns at the end of the final OMS target file only
        last_columns = [
            "qual",
            "filter_status",
            "filter_method",
            "caller",
        ]

        final_df = final_df[
            [col for col in final_df.columns if col not in last_columns] +
            [col for col in last_columns if col in final_df.columns]
        ]
        
        results_dir = os.path.join(RESULTS_BASE, biosample)
        os.makedirs(results_dir, exist_ok=True)

        final_target_path = os.path.join(results_dir, f"{biosample}_OMStarget.xlsx")
        final_df.to_excel(final_target_path, index=False)

        print(f"[OK] Combined OMStarget saved → {final_target_path}")

    print("[DONE] All callers processed.")


if __name__ == "__main__":
    main()
