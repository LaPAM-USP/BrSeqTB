#!/usr/bin/env bash
# ============================================================
# Variant calling with LoFreq (single-thread, no indelqual)
# Usage: ./lofreq.sh <biosample>
# Input:  bwa/<biosample>/<biosample>.bam
# Output: lofreq/<biosample>/<biosample>_lofreq.vcf.gz
# Reference: database/mtbRef/NC0009623.fasta
#
# Notes:
# - LoFreq is inherently single-threaded
# - Parallelization is handled by Nextflow per biosample
# ============================================================

export LC_ALL=C
set -euo pipefail

START_TIME=$SECONDS

BIOSAMPLE="${1:-}"

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

BWA_DIR="${PROJECT_DIR}/bwa/${BIOSAMPLE}"
REF="${PROJECT_DIR}/database/mtbRef/NC0009623.fasta"
OUTPUT_DIR="${PROJECT_DIR}/lofreq/${BIOSAMPLE}"

# ================== CHECK INPUT ==================
if [[ -z "$BIOSAMPLE" ]]; then
    echo "Usage: ./lofreq.sh <biosample>"
    exit 1
fi

if [[ ! -d "$BWA_DIR" ]]; then
    echo "[ERROR] BAM directory not found: ${BWA_DIR}"
    exit 1
fi

if [[ ! -f "$REF" ]]; then
    echo "[ERROR] Reference genome not found: ${REF}"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# ================== DEPENDENCY CHECKS ==================
for cmd in lofreq bgzip bcftools; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: $cmd"
        exit 1
    fi
done

echo "[RUN] Running LoFreq for biosample: ${BIOSAMPLE}"
echo "[REF] Reference genome: ${REF}"
echo "[OUT] Output directory: ${OUTPUT_DIR}"
echo "---------------------------------------------"

# ================== FIND BAM FILE ==================
BAM_FILE="${BWA_DIR}/${BIOSAMPLE}.bam"

if [[ ! -f "$BAM_FILE" ]]; then
    echo "[ERROR] BAM file not found: ${BAM_FILE}"
    exit 1
fi

# ================== DEFINE OUTPUT ==================
VCF_OUTPUT="${OUTPUT_DIR}/${BIOSAMPLE}_lofreq.vcf"
VCF_GZ="${VCF_OUTPUT}.gz"

# ================== SKIP IF RESULTS ALREADY EXIST ==================
if [[ -f "$VCF_GZ" ]] && [[ -f "${VCF_GZ}.csi" ]]; then
    echo "[SKIP] LoFreq output already exists:"
    echo "       $(basename "$VCF_GZ")"
    echo "       $(basename "${VCF_GZ}.csi")"
    exit 0
fi

# ================== RUN LoFreq ==================
echo "[RUN] Calling variants (single-thread)..."

lofreq call \
    --call-indels \
    --no-default-filter \
    -f "$REF" \
    -o "$VCF_OUTPUT" \
    "$BAM_FILE"

echo "[OK] Variant calling complete."

# ================== COMPRESS & INDEX ==================
bgzip -f "$VCF_OUTPUT"
bcftools index -f "$VCF_GZ"

echo "[OK] Compressed and indexed VCF."

# ================== COMPLETION ==================
ELAPSED=$(( SECONDS - START_TIME ))
printf "[DONE] LoFreq completed for %s (%02d min %02d sec)\n" \
    "$BIOSAMPLE" $((ELAPSED/60)) $((ELAPSED%60))

echo "[OUT] ${VCF_GZ}"

