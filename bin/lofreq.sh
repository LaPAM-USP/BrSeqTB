#!/usr/bin/env bash
# ============================================================
# Variant calling with LoFreq (parallel, with indelqual)
# Usage: ./lofreq.sh <biosample> [threads] [min_af]
# Input:  bwa/<biosample>/<biosample>.bam
# Output: lofreq/<biosample>/<biosample>_lofreq.vcf.gz
# Reference: database/mtbRef/NC0009623.fasta
#
# Notes:
# - LoFreq call-parallel is used for within-sample parallelization
# - Indel qualities are added before variant calling for INDEL support
# - Final variants are filtered by AF >= min_af
# - Default min_af is 0.05
# ============================================================

export LC_ALL=C
set -euo pipefail

START_TIME=$SECONDS

BIOSAMPLE="${1:-}"
THREADS="${2:-1}"
MIN_AF="${3:-0.00}"

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

BWA_DIR="${PROJECT_DIR}/bwa/${BIOSAMPLE}"
REF="${PROJECT_DIR}/database/mtbRef/NC0009623.fasta"
OUTPUT_DIR="${PROJECT_DIR}/lofreq/${BIOSAMPLE}"

# ================== VALIDATE MIN_AF ==================
if ! [[ "$MIN_AF" =~ ^([0-9]+([.][0-9]+)?|[.][0-9]+)$ ]]; then
    echo "[ERROR] MIN_AF must be a decimal number."
    echo "        Examples: 0, 0.05, 0.10, 0.5, 1.0"
    echo "        Received: ${MIN_AF}"
    exit 1
fi

if ! awk -v af="$MIN_AF" 'BEGIN { exit !(af >= 0 && af <= 1) }'; then
    echo "[ERROR] MIN_AF must be between 0 and 1."
    echo "        Received: ${MIN_AF}"
    exit 1
fi

# ================== CHECK INPUT ==================
if [[ -z "$BIOSAMPLE" ]]; then
    echo "Usage: ./lofreq.sh <biosample> [threads] [min_af]"
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
for cmd in lofreq samtools bgzip bcftools; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: $cmd"
        exit 1
    fi
done

echo "[RUN] Running LoFreq for biosample: ${BIOSAMPLE}"
echo "[REF] Reference genome: ${REF}"
echo "[OUT] Output directory: ${OUTPUT_DIR}"
echo "[THREADS] ${THREADS}"
echo "[MIN_AF] ${MIN_AF}"
echo "---------------------------------------------"

# ================== FIND BAM FILE ==================
BAM_FILE="${BWA_DIR}/${BIOSAMPLE}.bam"

if [[ ! -f "$BAM_FILE" ]]; then
    echo "[ERROR] BAM file not found: ${BAM_FILE}"
    exit 1
fi

if [[ ! -f "${BAM_FILE}.bai" ]] && [[ ! -f "${BAM_FILE%.bam}.bai" ]]; then
    echo "[ERROR] BAM index not found for: ${BAM_FILE}"
    echo "        Expected: ${BAM_FILE}.bai or ${BAM_FILE%.bam}.bai"
    exit 1
fi

# ================== DEFINE OUTPUT ==================
INDELQUAL_BAM="${OUTPUT_DIR}/${BIOSAMPLE}.indelqual.bam"
RAW_VCF_OUTPUT="${OUTPUT_DIR}/${BIOSAMPLE}_lofreq.raw.vcf"
VCF_OUTPUT="${OUTPUT_DIR}/${BIOSAMPLE}_lofreq.vcf"
VCF_GZ="${VCF_OUTPUT}.gz"

# ================== SKIP IF RESULTS ALREADY EXIST ==================
if [[ -f "$VCF_GZ" ]] && [[ -f "${VCF_GZ}.csi" ]]; then
    echo "[SKIP] LoFreq output already exists:"
    echo "       $(basename "$VCF_GZ")"
    echo "       $(basename "${VCF_GZ}.csi")"
    exit 0
fi

# ================== PREPARE REFERENCE INDEX ==================
if [[ ! -f "${REF}.fai" ]]; then
    echo "[RUN] FASTA index not found. Creating: ${REF}.fai"
    samtools faidx "$REF"
fi

# ================== ADD INDEL QUALITIES ==================
echo "[RUN] Adding indel qualities with LoFreq indelqual..."

lofreq indelqual \
    --dindel \
    -f "$REF" \
    -o "$INDELQUAL_BAM" \
    "$BAM_FILE"

rm -f "${INDELQUAL_BAM}.bai" "${INDELQUAL_BAM%.bam}.bai"
samtools index -@ "$THREADS" "$INDELQUAL_BAM"

echo "[OK] Indel qualities added."

# ================== RUN LoFreq ==================
echo "[RUN] Calling variants with LoFreq call-parallel..."

lofreq call-parallel \
    --pp-threads "$THREADS" \
    --call-indels \
    -m 60 \
    -Q 30 \
    -f "$REF" \
    -o "$RAW_VCF_OUTPUT" \
    "$INDELQUAL_BAM"

echo "[OK] Raw variant calling complete."

# ================== FILTER BY ALLELE FREQUENCY ==================
echo "[RUN] Filtering variants by allele frequency >= ${MIN_AF}..."

lofreq filter \
    -a "$MIN_AF" \
    -i "$RAW_VCF_OUTPUT" \
    -o "$VCF_OUTPUT"

echo "[OK] Allele frequency filtering complete."

# ================== COMPRESS & INDEX ==================
bgzip -f "$VCF_OUTPUT"
bcftools index -f "$VCF_GZ"

echo "[OK] Compressed and indexed VCF."

# ================== COMPLETION ==================
ELAPSED=$(( SECONDS - START_TIME ))
printf "[DONE] LoFreq completed for %s (%02d min %02d sec)\n" \
    "$BIOSAMPLE" $((ELAPSED/60)) $((ELAPSED%60))

echo "[OUT] ${VCF_GZ}"
