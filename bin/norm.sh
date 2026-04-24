#!/usr/bin/env bash
# ============================================================
# Normalize and decompose variants using bcftools + vt
# Usage: ./norm.sh <biosample ID>
# Input:  gatk/<biosample ID>/*_gatk.vcf.gz
# Output: norm/<biosample ID>/<biosample ID>_norm.vcf.gz
# Reference: database/mtbRef/NC0009623.fasta
#
# Notes:
# - This script is designed to run on one biosample at a time
# - Parallelization is handled externally by Nextflow
# - Uses tools from Conda environment (activated by Nextflow)
# ============================================================

export LC_ALL=C
set -euo pipefail

# ================= PROJECT DIR =================
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ================= INPUTS =================
BIOSAMPLE="${1:-}"

INPUT_DIR="${PROJECT_DIR}/gatk/${BIOSAMPLE}"
REF="${PROJECT_DIR}/database/mtbRef/NC0009623.fasta"
OUTPUT_DIR="${PROJECT_DIR}/norm/${BIOSAMPLE}"

# ================= CHECK INPUT =================
if [[ -z "$BIOSAMPLE" ]]; then
    echo "Usage: ./norm.sh <biosample>"
    exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "[ERROR] Input directory not found: ${INPUT_DIR}"
    exit 1
fi

if [[ ! -f "$REF" ]]; then
    echo "[ERROR] Reference genome not found: ${REF}"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# ================= OUTPUT (SKIP IF EXISTS) =================
FINAL_VCF="${OUTPUT_DIR}/${BIOSAMPLE}_norm.vcf.gz"

if [[ -f "$FINAL_VCF" ]]; then
    echo "[SKIP] Normalized VCF already exists: ${FINAL_VCF}"
    exit 0
fi

# ================= DEPENDENCY CHECKS =================
for cmd in bcftools vt bgzip; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found in Conda environment: $cmd"
        exit 1
    fi
done

echo "[RUN] Normalization and decomposition for biosample: ${BIOSAMPLE}"
echo "[IN]  Input directory: ${INPUT_DIR}"
echo "[REF] Reference genome: ${REF}"
echo "[OUT] Output directory: ${OUTPUT_DIR}"
echo "---------------------------------------------"

# ================= LOCATE INPUT VCF =================
VCF_FILE=$(find "$INPUT_DIR" -type f -name "*_gatk.vcf.gz" | head -n 1)

if [[ -z "$VCF_FILE" ]]; then
    echo "[ERROR] No input VCF (.vcf.gz) found in ${INPUT_DIR}"
    exit 1
fi

echo "[INFO] Using input VCF: $(basename "$VCF_FILE")"

# ================= RUN NORMALIZATION =================
OUTPUT_VCF="${OUTPUT_DIR}/${BIOSAMPLE}_norm.vcf"

echo "[RUN] Normalizing and decomposing variants..."
bcftools norm --fasta-ref "$REF" -m-any "$VCF_FILE" | \
vt decompose_blocksub - -o "$OUTPUT_VCF"

echo "[OK] Normalization complete: ${OUTPUT_VCF}"

# ================= COMPRESS AND INDEX =================
echo "[RUN] Compressing and indexing normalized VCF..."
bgzip -f "$OUTPUT_VCF"
bcftools index -f "${OUTPUT_VCF}.gz"

echo "[OK] Compressed and indexed: ${OUTPUT_VCF}.gz"

# ================= COMPLETION =================
echo "---------------------------------------------"
echo "[DONE] Normalization and decomposition completed for biosample: ${BIOSAMPLE}"
echo "[OUT] Final file: ${OUTPUT_VCF}.gz"
echo "[OUT] Index: ${OUTPUT_VCF}.gz.csi"
echo "[OUT] Folder: ${OUTPUT_DIR}/"
