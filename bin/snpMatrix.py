#!/usr/bin/env python3
# ============================================================
# SNP Matrix Generator + MixInfection Table Extractor
#
# This script processes SnpEff-annotated and normalized VCF files
# to produce:
#
#   1. snpMatrix/snpmatrix.fasta
#        → multi-sample SNP matrix FASTA
#
#      A position is removed from the entire SNP matrix if any
#      sample has, at that position:
#
#        • a heterozygous genotype
#        • an INDEL or another non-SNP allele
#        • a missing or partially missing genotype
#        • a variant in a forbidden gene
#        • DP lower than MIN_DP (7)
#        • missing DP
#
#      If a sample has no VCF record at a retained position,
#      the reference base is used, assuming absence from the
#      variant-only VCF means homozygous reference.
#
#   2. mixInfection/<sample>/<sample>.tsv
#        → full SNP table including heterozygotes, suitable
#          for MixInfect or mixed-infection analysis.
#
#      The mixInfection logic is unchanged from the original
#      pipeline version.
#
# Variants in forbidden genes are removed. Forbidden genes are
# defined in:
#
#     database/mtbRef/forbidden_genes.txt
#
# The script removes variants if:
#   • Gene_Name matches exactly any forbidden gene
#   • Gene_ID (RvXXXX) matches exactly any forbidden gene
#
# Matching is case-insensitive.
#
# ============================================================

import os
import sys
from pathlib import Path
import pysam
from Bio import SeqIO


# ============================================================
# PROJECT STRUCTURE
# ============================================================
PROJECT_DIR = Path(
    os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
)

MANIFEST        = PROJECT_DIR / "manifest.tsv"
SNPEFF_DIR      = PROJECT_DIR / "snpeff"
REF_FILE        = PROJECT_DIR / "database" / "mtbRef" / "NC0009623.fasta"
FORBIDDEN_FILE  = PROJECT_DIR / "database" / "mtbRef" / "forbidden_genes.txt"

OUT_FASTA_DIR = PROJECT_DIR / "snpMatrix"
OUT_MIX_DIR   = PROJECT_DIR / "mixInfection"

MIN_DP = 7


# ============================================================
# LOAD MANIFEST (SOURCE OF TRUTH)
# ============================================================
if not MANIFEST.is_file():
    sys.exit(f"[ERROR] Manifest not found: {MANIFEST}")

biosamples = []

with open(MANIFEST) as f:

    header = f.readline()

    for line in f:

        s = line.strip()

        if s:
            biosamples.append(s)

if not biosamples:
    sys.exit("[ERROR] Manifest is empty or invalid")


# ============================================================
# LOAD FORBIDDEN GENES
# ============================================================
def load_forbidden_genes(path: Path):

    if not path.is_file():
        sys.exit(f"[ERROR] Forbidden genes file not found: {path}")

    forbidden = set()

    with open(path) as f:

        for line in f:

            line = line.strip()

            if not line or line.startswith("#"):
                continue

            forbidden.add(line.lower())

    return forbidden


FORBIDDEN_SET = load_forbidden_genes(FORBIDDEN_FILE)


# ============================================================
# LOAD REFERENCE
# ============================================================
def load_reference(path: Path):

    if not path.is_file():
        sys.exit(f"[ERROR] Reference FASTA not found: {path}")

    ref = SeqIO.to_dict(
        SeqIO.parse(str(path), "fasta")
    )

    if len(ref) != 1:
        sys.exit(
            "[ERROR] Reference FASTA must contain exactly one sequence"
        )

    key = list(ref.keys())[0]

    return key, list(str(ref[key].seq))


# ============================================================
# ANN / FORBIDDEN LOGIC
# ============================================================
def extract_genes_from_ann(ann_tuple):

    genes = []

    for ann in ann_tuple:

        fields = ann.split("|")

        if len(fields) >= 5:

            gene_name = fields[3].strip().lower()
            gene_id   = fields[4].strip().lower()

            genes.append(
                (gene_name, gene_id)
            )

    return genes


def is_forbidden_ann(ann_tuple):

    genes = extract_genes_from_ann(ann_tuple)

    for gene_name, gene_id in genes:

        if gene_name in FORBIDDEN_SET:
            return True

        if gene_id in FORBIDDEN_SET:
            return True

    return False


# ============================================================
# EXTRACT DP FROM A VCF RECORD
# ============================================================
def get_record_dp(rec, sample_data):
    """
    Obtain depth for SNP-matrix filtering.

    Priority:
        1. FORMAT/DP for the sample
        2. INFO/DP for the record

    Returns:
        Integer/float DP, or None when unavailable or invalid.
    """

    dp = sample_data.get("DP")

    if dp is None:
        dp = rec.info.get("DP")

    if dp is None:
        return None

    # Some VCF fields may be represented as a tuple.
    if isinstance(dp, (tuple, list)):

        if len(dp) == 0:
            return None

        dp = dp[0]

    try:
        return float(dp)

    except (TypeError, ValueError):
        return None


# ============================================================
# SNP EXTRACTION FOR SNP-MATRIX CONSTRUCTION
# ============================================================
def parse_vcf_for_matrix(vcf_path: Path):
    """
    Parse one VCF specifically for SNP-matrix construction.

    Returns:

        valid_snps:
            Dictionary of valid homozygous alternative SNPs:

                {position_0_based: alternative_base}

        candidate_positions:
            Positions containing a valid homozygous alternative
            SNP in this sample.

        invalid_positions:
            Positions that must be removed globally because this
            sample has, at that position:

                - an INDEL or another non-SNP allele
                - a forbidden-gene annotation
                - a missing or partially missing genotype
                - a heterozygous genotype
                - an invalid allele index
                - DP < MIN_DP
                - missing DP
    """

    valid_snps = {}
    candidate_positions = set()
    invalid_positions = set()

    vcf = pysam.VariantFile(str(vcf_path))

    for rec in vcf.fetch():

        pos0 = rec.pos - 1

        # ----------------------------------------------------
        # Remove INDELs and other non-SNP records globally
        # ----------------------------------------------------
        if rec.ref is None or len(rec.ref) != 1:

            invalid_positions.add(pos0)
            continue

        if not rec.alts:

            invalid_positions.add(pos0)
            continue

        if any(
            alt is None or len(alt) != 1
            for alt in rec.alts
        ):

            invalid_positions.add(pos0)
            continue

        # ----------------------------------------------------
        # Remove positions in forbidden genes globally
        # ----------------------------------------------------
        ann = rec.info.get("ANN", ())

        if ann and is_forbidden_ann(ann):

            invalid_positions.add(pos0)
            continue

        # ----------------------------------------------------
        # Require one sample genotype
        # ----------------------------------------------------
        if len(rec.samples) == 0:

            invalid_positions.add(pos0)
            continue

        sample_data = next(
            iter(rec.samples.values())
        )

        gt = sample_data.get("GT")

        # ----------------------------------------------------
        # Missing genotype
        # ----------------------------------------------------
        if gt is None:

            invalid_positions.add(pos0)
            continue

        if len(gt) == 0:

            invalid_positions.add(pos0)
            continue

        if any(
            allele is None
            for allele in gt
        ):

            invalid_positions.add(pos0)
            continue

        # ----------------------------------------------------
        # Heterozygous genotype
        #
        # Examples:
        #   0/1
        #   1/0
        #   1/2
        #
        # Works for haploid/diploid representations:
        #   GT=(1,) is not heterozygous.
        # ----------------------------------------------------
        if len(set(gt)) > 1:

            invalid_positions.add(pos0)
            continue

        # ----------------------------------------------------
        # Depth filter
        # ----------------------------------------------------
        dp = get_record_dp(
            rec,
            sample_data
        )

        if dp is None:

            invalid_positions.add(pos0)
            continue

        if dp < MIN_DP:

            invalid_positions.add(pos0)
            continue

        # ----------------------------------------------------
        # Homozygous/reference genotype
        #
        # Examples:
        #   0
        #   0/0
        #
        # Valid call, but not a candidate SNP.
        # ----------------------------------------------------
        if all(
            allele == 0
            for allele in gt
        ):
            continue

        # ----------------------------------------------------
        # Homozygous alternative genotype
        #
        # Examples:
        #   1
        #   1/1
        # ----------------------------------------------------
        alt_idx = gt[0]

        if (
            alt_idx < 1
            or alt_idx >= len(rec.alleles)
        ):

            invalid_positions.add(pos0)
            continue

        alt_base = rec.alleles[alt_idx]

        if (
            alt_base is None
            or len(alt_base) != 1
        ):

            invalid_positions.add(pos0)
            continue

        valid_snps[pos0] = alt_base
        candidate_positions.add(pos0)

    vcf.close()

    return (
        valid_snps,
        candidate_positions,
        invalid_positions
    )


# ============================================================
# FULL SNP TABLE (MIXED INFECTION)
# ORIGINAL PIPELINE LOGIC UNCHANGED
# ============================================================
def extract_snp_table(vcf_path: Path):

    rows = []

    vcf = pysam.VariantFile(
        str(vcf_path)
    )

    for rec in vcf.fetch():

        if any(
            len(a) != 1
            for a in rec.alleles
        ):
            continue

        ann = rec.info.get("ANN", ())

        if ann and is_forbidden_ann(ann):
            continue

        sample = list(
            rec.samples.values()
        )[0]

        gt = sample.get("GT")

        if gt is None:
            continue

        if all(
            a == 0
            for a in gt
        ):
            continue

        gt_str = "/".join(
            "." if a is None else str(a)
            for a in gt
        )

        alt_idx = (
            gt[0]
            if gt[0] != 0
            else gt[1]
        )

        alt_base = rec.alleles[alt_idx]

        rows.append({
            "chrom": rec.chrom,
            "pos": rec.pos,
            "ref": rec.ref,
            "alt": alt_base,
            "gt": gt_str
        })

    return rows


# ============================================================
# BUILD SNP-MATRIX SEQUENCE
# ============================================================
def build_matrix_sequence(
    ref_seq,
    positions,
    sample_snps
):
    """
    Build one SNP-matrix sequence using only globally retained
    positions.

    For each retained position:

        - ALT is used when this sample has a valid homozygous SNP;

        - reference is used when this sample has no VCF record
          at that position.

    The positions supplied to this function have already passed
    the global cohort filtering.
    """

    sequence = []

    for pos0 in positions:

        if pos0 in sample_snps:

            sequence.append(
                sample_snps[pos0]
            )

        else:

            sequence.append(
                ref_seq[pos0]
            )

    return "".join(sequence)


# ============================================================
# MAIN (COHORT)
# ============================================================
def main():

    OUT_FASTA_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    OUT_MIX_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    if not SNPEFF_DIR.is_dir():

        sys.exit(
            f"[ERROR] snpeff directory not found: {SNPEFF_DIR}"
        )

    # --------------------------------------------------------
    # Load reference genome
    # --------------------------------------------------------
    print(
        "[RUN] Loading reference genome..."
    )

    ref_name, ref_seq = load_reference(
        REF_FILE
    )

    print(
        f"[INFO] Reference: {ref_name}"
    )

    print(
        f"[INFO] Reference length: {len(ref_seq)} bp"
    )

    print(
        f"[INFO] Minimum DP for SNP matrix: {MIN_DP}"
    )

    print(
        f"[INFO] Biosamples from manifest: "
        f"{', '.join(biosamples)}"
    )

    # --------------------------------------------------------
    # Valid homozygous ALT SNPs for each sample
    # --------------------------------------------------------
    all_valid_snps = {}

    # --------------------------------------------------------
    # Positions containing a valid homozygous ALT SNP in at
    # least one sample
    # --------------------------------------------------------
    candidate_positions = set()

    # --------------------------------------------------------
    # Positions invalidated by at least one sample
    # --------------------------------------------------------
    invalid_positions = set()

    print(
        "[RUN] Parsing VCFs (cohort mode)..."
    )

    for sample in biosamples:

        sample_dir = (
            SNPEFF_DIR / sample
        )

        if not sample_dir.is_dir():

            sys.exit(
                f"[ERROR] Missing snpeff directory "
                f"for biosample: {sample}"
            )

        vcf_files = sorted(
            sample_dir.glob(
                "*_norm.vcf.gz"
            )
        )

        if not vcf_files:

            sys.exit(
                f"[ERROR] Missing *_norm.vcf.gz "
                f"for biosample: {sample}"
            )

        vcf_path = vcf_files[0]

        # ====================================================
        # SNP-MATRIX-SPECIFIC PARSING
        # ====================================================
        (
            valid_snps,
            sample_candidate_positions,
            sample_invalid_positions
        ) = parse_vcf_for_matrix(
            vcf_path
        )

        all_valid_snps[sample] = (
            valid_snps
        )

        candidate_positions.update(
            sample_candidate_positions
        )

        invalid_positions.update(
            sample_invalid_positions
        )

        print(
            f"[INFO] {sample}: "
            f"{len(valid_snps)} valid homozygous SNPs; "
            f"{len(sample_invalid_positions)} invalid positions"
        )

        # ====================================================
        # MIXED-INFECTION TABLE
        # Original pipeline logic unchanged
        # ====================================================
        rows = extract_snp_table(
            vcf_path
        )

        outdir = (
            OUT_MIX_DIR / sample
        )

        outdir.mkdir(
            parents=True,
            exist_ok=True
        )

        outfile = (
            outdir / f"{sample}.tsv"
        )

        with open(outfile, "w") as out:

            out.write(
                "chrom\tpos\tref\talt\tgt\n"
            )

            for r in rows:

                out.write(
                    f"{r['chrom']}\t"
                    f"{r['pos']}\t"
                    f"{r['ref']}\t"
                    f"{r['alt']}\t"
                    f"{r['gt']}\n"
                )

        print(
            f"[OK] mixInfection table written: "
            f"{outfile}"
        )

    # ========================================================
    # Check if any valid homozygous SNP was found
    # ========================================================
    if not candidate_positions:

        sys.exit(
            "[ERROR] No valid homozygous SNPs found in cohort"
        )

    # ========================================================
    # GLOBAL SNP-MATRIX FILTERING
    #
    # Remove every candidate position that was invalid in at
    # least one sample/VCF record.
    #
    # Therefore:
    #
    #     retained positions
    #         =
    #     candidate positions
    #         -
    #     globally invalid positions
    #
    # A sample without a VCF record at a retained position
    # receives the reference nucleotide.
    # ========================================================
    excluded_candidate_positions = (
        candidate_positions
        & invalid_positions
    )

    valid_positions = (
        candidate_positions
        - invalid_positions
    )

    positions_sorted = sorted(
        valid_positions
    )

    print("")

    print(
        "[INFO] SNP-matrix global filtering:"
    )

    print(
        f"[INFO] Candidate SNP positions: "
        f"{len(candidate_positions)}"
    )

    print(
        f"[INFO] Candidate positions excluded globally: "
        f"{len(excluded_candidate_positions)}"
    )

    print(
        f"[INFO] Final retained SNP positions: "
        f"{len(positions_sorted)}"
    )

    if not positions_sorted:

        sys.exit(
            "[ERROR] No SNP positions remained "
            "after global filtering"
        )

    # ========================================================
    # WRITE SNP MATRIX FASTA
    # ========================================================
    fasta_path = (
        OUT_FASTA_DIR
        / "snpmatrix.fasta"
    )

    print(
        "[RUN] Writing SNP matrix FASTA..."
    )

    with open(fasta_path, "w") as o:

        for sample in biosamples:

            sample_snps = (
                all_valid_snps.get(
                    sample,
                    {}
                )
            )

            matrix_sequence = (
                build_matrix_sequence(
                    ref_seq=ref_seq,
                    positions=positions_sorted,
                    sample_snps=sample_snps
                )
            )

            o.write(
                f">{sample}\n"
                f"{matrix_sequence}\n"
            )

    print(
        f"[DONE] SNP matrix written: "
        f"{fasta_path}"
    )

    print(
        f"[DONE] SNP-matrix alignment length: "
        f"{len(positions_sorted)} positions"
    )

    print(
        "[DONE] mixInfection tables generated."
    )


# ============================================================
if __name__ == "__main__":
    main()
