#!/usr/bin/env bash
# ============================================================
# Variant annotation with SnpEff per biosample
# Usage: ./snpeff.sh <biosample>
#
# Uses:
#   - Custom DB: NC_0009623
#   - snpEff from Conda/Mamba environment (Nextflow)
# ============================================================

export LC_ALL=C
set -euo pipefail

BIOSAMPLE="${1:-}"

if [[ -z "$BIOSAMPLE" ]]; then
    echo "Usage: ./snpeff.sh <biosample>"
    exit 1
fi

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

GENOME_DB="NC_0009623"
CALLERS=("gatk" "lofreq" "norm" "delly")

# -------------------- Locate SnpEff (from PATH) --------------------
if ! command -v snpEff >/dev/null 2>&1; then
    echo "[ERROR] snpEff not found in PATH (Conda environment not active?)"
    exit 1
fi

SNPEFF_BIN="$(which snpEff)"

SNPEFF_PROJECT_DIR="${PROJECT_DIR}/database/snpeff"
SNPEFF_CONFIG="${SNPEFF_PROJECT_DIR}/snpEff.config"
SNPEFF_DATA_ROOT="${SNPEFF_PROJECT_DIR}/data"

# -------------------- DEPENDENCY CHECKS --------------------
for cmd in bgzip tabix; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: $cmd"
        exit 1
    fi
done

[[ -f "$SNPEFF_CONFIG" ]] || {
    echo "[ERROR] snpEff.config not found:"
    echo "        ${SNPEFF_CONFIG}"
    exit 1
}

OUTPUT_DIR="${PROJECT_DIR}/snpeff/${BIOSAMPLE}"
mkdir -p "$OUTPUT_DIR"

echo "[INFO] snpEff binary : $SNPEFF_BIN"
echo "[INFO] Config file   : $SNPEFF_CONFIG"
echo "[INFO] Database      : ${GENOME_DB}"
echo "---------------------------------------------"

# ============================================================
# STRICT INPUT CHECK — all 4 callers MUST exist
# ============================================================
declare -A CALLER_VCF

for CALLER in "${CALLERS[@]}"; do
    INPUT_DIR="${PROJECT_DIR}/${CALLER}/${BIOSAMPLE}"

    [[ -d "$INPUT_DIR" ]] || {
        echo "[ERROR] Missing caller directory: ${INPUT_DIR}"
        exit 1
    }

    mapfile -t VCF_FILES < <(
        find "$INPUT_DIR" -type f -name "*.vcf.gz" ! -name "*.g.vcf.gz" | sort
    )

    [[ ${#VCF_FILES[@]} -gt 0 ]] || {
        echo "[ERROR] No VCF found for caller '${CALLER}' in ${INPUT_DIR}"
        exit 1
    }

    CALLER_VCF["$CALLER"]="${VCF_FILES[0]}"
done

# ============================================================
# GLOBAL SKIP — only if ALL 4 are already annotated
# ============================================================
ALL_DONE=true

for CALLER in "${CALLERS[@]}"; do
    VCF_FILE="${CALLER_VCF[$CALLER]}"
    VCF_NAME=$(basename "$VCF_FILE")
    ANNOTATED="${OUTPUT_DIR}/${VCF_NAME}"

    if [[ ! -f "$ANNOTATED" || ! -f "${ANNOTATED}.tbi" ]]; then
        ALL_DONE=false
        break
    fi
done

if [[ "$ALL_DONE" == true ]]; then
    echo "[SKIP] All 4 callers already annotated for ${BIOSAMPLE}"
    exit 0
fi

echo "[RUN] Starting SnpEff annotation for biosample: ${BIOSAMPLE}"
echo "---------------------------------------------"

# ============================================================
# PROCESS — annotate only what is missing
# ============================================================
for CALLER in "${CALLERS[@]}"; do
    VCF_FILE="${CALLER_VCF[$CALLER]}"
    VCF_NAME=$(basename "$VCF_FILE")

    OUTPUT_FILE="${OUTPUT_DIR}/${VCF_NAME}"
    LOG_FILE="${OUTPUT_DIR}/${VCF_NAME%.vcf.gz}_snpeff.log"

    if [[ -f "$OUTPUT_FILE" && -f "${OUTPUT_FILE}.tbi" ]]; then
        echo "[SKIP] Already annotated (${CALLER}): ${VCF_NAME}"
        continue
    fi

    echo "[RUN] Annotating (${CALLER}): ${VCF_NAME}"

    TMP_VCF=$(mktemp --suffix=".vcf")
    gunzip -c "$VCF_FILE" > "$TMP_VCF"

    snpEff eff \
        -c "$SNPEFF_CONFIG" \
        -dataDir "$SNPEFF_DATA_ROOT" \
        -v \
        -ud 100 \
        "$GENOME_DB" \
        "$TMP_VCF" \
        2> "$LOG_FILE" | bgzip -c > "$OUTPUT_FILE"

    tabix -f -p vcf "$OUTPUT_FILE"
    rm -f "$TMP_VCF"

    echo "[OK] Annotated: ${OUTPUT_FILE}"
    echo "---------------------------------------------"
done

echo "[DONE] SnpEff annotation completed for biosample: ${BIOSAMPLE}"
echo "[OUT]  Results in: ${OUTPUT_DIR}/"
