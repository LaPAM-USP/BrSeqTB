#!/usr/bin/env bash
# ============================================================
# Structural Variant Calling with DELLY
# ============================================================

export LC_ALL=C
set -euo pipefail

BIOSAMPLE="${1:-}"

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

BWA_DIR="${PROJECT_DIR}/bwa/${BIOSAMPLE}"
REF="${PROJECT_DIR}/database/mtbRef/NC0009623.fasta"
OUTPUT_DIR="${PROJECT_DIR}/delly/${BIOSAMPLE}"
LOG_FILE="${OUTPUT_DIR}/${BIOSAMPLE}_delly.log"

if [[ -z "$BIOSAMPLE" ]]; then
    echo "Usage: ./delly.sh <biosample>"
    exit 1
fi

[[ -d "$BWA_DIR" ]] || { echo "[ERROR] BAM directory not found: $BWA_DIR"; exit 1; }
[[ -f "$REF" ]] || { echo "[ERROR] Reference genome not found: $REF"; exit 1; }

mkdir -p "$OUTPUT_DIR"

# ==================== DEPENDENCIES ====================
for cmd in delly bcftools bgzip tabix; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "[ERROR] Required command not found: $cmd"
        exit 1
    }
done

BAM_FILE="${BWA_DIR}/${BIOSAMPLE}.bam"
[[ -f "$BAM_FILE" ]] || { echo "[ERROR] BAM file not found: $BAM_FILE"; exit 1; }

# ==================== SKIP IF ALREADY DONE ====================
VCF_FILE="${OUTPUT_DIR}/${BIOSAMPLE}_delly.vcf.gz"

if [[ -f "$VCF_FILE" && -f "${VCF_FILE}.tbi" ]]; then
    echo "[SKIP] DELLY output already exists:"
    echo "       $(basename "$VCF_FILE")"
    echo "       $(basename "${VCF_FILE}.tbi")"
    exit 0
fi

echo "[RUN] DELLY for ${BIOSAMPLE}"
echo "---------------------------------------------"

RAW_BCF="${OUTPUT_DIR}/${BIOSAMPLE}_delly.bcf"

# ==================== RUN DELLY ====================
delly call \
    -g "$REF" \
    -o "$RAW_BCF" \
    "$BAM_FILE" \
    2>&1 | tee "$LOG_FILE"

[[ -s "$RAW_BCF" ]] || { echo "[ERROR] DELLY did not produce a valid BCF."; exit 1; }

# ==================== CONVERT TO VCF ====================
bcftools view "$RAW_BCF" -Ov | bgzip -c > "$VCF_FILE"
tabix -f -p vcf "$VCF_FILE"

# ==================== CLEANUP ====================
rm -f "$RAW_BCF" "${RAW_BCF}.csi" 2>/dev/null || true

echo "[DONE] DELLY completed for ${BIOSAMPLE}"
echo "[OUT]  $VCF_FILE"
echo "[LOG]  $LOG_FILE"

