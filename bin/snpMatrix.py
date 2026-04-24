#!/usr/bin/env python3
# ============================================================
# SNP Matrix Generator + MixInfection Table Extractor
#
# This script processes SnpEff-annotated and normalized VCF files
# to produce:
#
#   1. snpMatrix/snpmatrix.fasta
#        → multi-sample SNP matrix FASTA (no heterozygotes)
#
#   2. mixInfection/<sample>/<sample>.tsv
#        → full SNP table (including heterozygotes), suitable
#          for MixInfect or mixed infection analysis.
#
# Variants in forbidden genes are removed. Forbidden genes are
# defined in:
#
#     dataset/mtbRef/forbidden_genes.txt
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

    ref = SeqIO.to_dict(SeqIO.parse(str(path), "fasta"))
    if len(ref) != 1:
        sys.exit("[ERROR] Reference FASTA must contain exactly one sequence")

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
            genes.append((gene_name, gene_id))
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
# SNP EXTRACTION
# ============================================================
def extract_snps(vcf_path: Path, allow_het: bool):

    snps = []
    vcf = pysam.VariantFile(str(vcf_path))

    for rec in vcf.fetch():

        if any(len(a) != 1 for a in rec.alleles):
            continue

        ann = rec.info.get("ANN", ())
        if ann and is_forbidden_ann(ann):
            continue

        sample = list(rec.samples.values())[0]
        gt = sample.get("GT")
        if gt is None:
            continue

        is_het = (len(gt) == 2 and gt[0] != gt[1])
        if is_het and not allow_het:
            continue

        if all(a == 0 for a in gt):
            continue

        alt_idx  = gt[0] if gt[0] != 0 else gt[1]
        alt_base = rec.alleles[alt_idx]

        pos0 = rec.pos - 1
        snps.append((pos0, alt_base))

    return snps

# ============================================================
# FULL SNP TABLE (MIXED INFECTION)
# ============================================================
def extract_snp_table(vcf_path: Path):

    rows = []
    vcf = pysam.VariantFile(str(vcf_path))

    for rec in vcf.fetch():

        if any(len(a) != 1 for a in rec.alleles):
            continue

        ann = rec.info.get("ANN", ())
        if ann and is_forbidden_ann(ann):
            continue

        sample = list(rec.samples.values())[0]
        gt = sample.get("GT")
        if gt is None:
            continue

        if all(a == 0 for a in gt):
            continue

        gt_str = "/".join("." if a is None else str(a) for a in gt)
        alt_idx  = gt[0] if gt[0] != 0 else gt[1]
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
# BUILD GENOME
# ============================================================
def build_modified_genome(ref_seq, snps):

    seq = ref_seq[:]
    for pos, base in snps:
        if 0 <= pos < len(seq):
            seq[pos] = base
    return seq

# ============================================================
# MAIN (COHORT)
# ============================================================
def main():

    OUT_FASTA_DIR.mkdir(parents=True, exist_ok=True)
    OUT_MIX_DIR.mkdir(parents=True, exist_ok=True)

    if not SNPEFF_DIR.is_dir():
        sys.exit(f"[ERROR] snpeff directory not found: {SNPEFF_DIR}")

    print("[RUN] Loading reference genome...")
    ref_name, ref_seq = load_reference(REF_FILE)

    print(f"[INFO] Biosamples from manifest: {', '.join(biosamples)}")

    all_snps_nohet = {}
    all_positions = set()

    print("[RUN] Parsing VCFs (cohort mode)...")

    for sample in biosamples:

        sample_dir = SNPEFF_DIR / sample
        if not sample_dir.is_dir():
            sys.exit(f"[ERROR] Missing snpeff directory for biosample: {sample}")

        vcf_files = sorted(sample_dir.glob("*_norm.vcf.gz"))
        if not vcf_files:
            sys.exit(f"[ERROR] Missing *_norm.vcf.gz for biosample: {sample}")

        vcf_path = vcf_files[0]

        snps_nohet = extract_snps(vcf_path, allow_het=False)
        all_snps_nohet[sample] = snps_nohet

        for pos, _ in extract_snps(vcf_path, allow_het=True):
            all_positions.add(pos)

        rows = extract_snp_table(vcf_path)
        outdir = OUT_MIX_DIR / sample
        outdir.mkdir(parents=True, exist_ok=True)
        outfile = outdir / f"{sample}.tsv"

        with open(outfile, "w") as out:
            out.write("chrom\tpos\tref\talt\tgt\n")
            for r in rows:
                out.write(
                    f"{r['chrom']}\t{r['pos']}\t{r['ref']}\t{r['alt']}\t{r['gt']}\n"
                )

        print(f"[OK] mixInfection table written: {outfile}")

    if not all_positions:
        sys.exit("[ERROR] No SNPs found in cohort")

    positions_sorted = sorted(all_positions)

    fasta_path = OUT_FASTA_DIR / "snpmatrix.fasta"
    print("[RUN] Writing SNP matrix FASTA...")

    with open(fasta_path, "w") as o:
        for sample in biosamples:
            snps = all_snps_nohet.get(sample, [])
            seq = build_modified_genome(ref_seq, snps)
            trimmed = ''.join(seq[p] for p in positions_sorted)
            o.write(f">{sample}\n{trimmed}\n")

    print(f"[DONE] SNP matrix written: {fasta_path}")
    print("[DONE] mixInfection tables generated.")

# ============================================================
if __name__ == "__main__":
    main()

