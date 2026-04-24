#!/usr/bin/env bash
# ============================================================
# Variant calling with GATK HaplotypeCaller (GVCF)
# Usage: ./gatkGvcf.sh <biosample>
# Input:  bwa/<biosample>/<biosample>.bam
# Output: gatk/<biosample>/*.g.vcf.gz
# Reference: database/mtbRef/NC0009623.fasta
# ============================================================

export LC_ALL=C
set -euo pipefail

START_TIME=$SECONDS

BIOSAMPLE="${1:-}"
THREADS="${NXF_TASK_CPUS:-1}"

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

BWA_DIR="${PROJECT_DIR}/bwa/${BIOSAMPLE}"
REF="${PROJECT_DIR}/database/mtbRef/NC0009623.fasta"
OUTPUT_DIR="${PROJECT_DIR}/gatk/${BIOSAMPLE}"


# ================== CHECK INPUT ==================
if [[ -z "$BIOSAMPLE" ]]; then
    echo "Usage: ./gatkGvcf.sh <biosample>"
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

# ================== DEPENDENCIES ==================
for cmd in gatk samtools; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: $cmd"
        exit 1
    fi
done

# ================== DEFINE OUTPUTS ==================
GVCF="${OUTPUT_DIR}/${BIOSAMPLE}.g.vcf.gz"

# ================== SKIP IF RESULTS ALREADY EXIST ==================
if [[ -f "$GVCF" && -f "${GVCF}.tbi" ]]; then
    echo "[SKIP] GATK GVCF output already exists for ${BIOSAMPLE}:"
    echo "       - $(basename "$GVCF")"
    exit 0
fi

# ================== FIND BAM ==================
BAM_FILE="${BWA_DIR}/${BIOSAMPLE}.bam"

if [[ ! -f "$BAM_FILE" ]]; then
    echo "[ERROR] BAM file not found: ${BAM_FILE}"
    exit 1
fi

echo "[RUN] GATK HaplotypeCaller (GVCF) for ${BIOSAMPLE}"
echo "[INFO] Threads: ${THREADS}"
echo "---------------------------------------------"

# ============================================================
# GVCF MODE
# ============================================================
LOG_GVCF="${OUTPUT_DIR}/${BIOSAMPLE}_gvcf.log"

gatk HaplotypeCaller \
    -R "$REF" \
    -I "$BAM_FILE" \
    -O "$GVCF" \
    -ERC GVCF \
    --native-pair-hmm-threads "$THREADS" \
    -A FisherStrand \
    -A StrandOddsRatio \
    -A QualByDepth \
    -A MappingQualityRankSumTest \
    -A ReadPosRankSumTest \
    -A DepthPerAlleleBySample \
    -A Coverage \
    2>&1 | tee "$LOG_GVCF"

# ================== COMPLETION ==================
ELAPSED=$(( SECONDS - START_TIME ))
printf "[DONE] GATK GVCF completed for %s (%02d min %02d sec)\n" \
    "$BIOSAMPLE" $((ELAPSED/60)) $((ELAPSED%60))

echo "[OUT] GVCF: $GVCF"

