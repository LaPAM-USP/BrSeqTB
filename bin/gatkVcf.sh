#!/usr/bin/env bash
# ============================================================
# Variant calling with GATK HaplotypeCaller (VCF)
# Usage: ./gatkVcf.sh <biosample>
# Input:  bwa/<biosample>/<biosample>.bam
# Output: gatk/<biosample>/*_gatk.vcf.gz
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
    echo "Usage: ./gatkVcf.sh <biosample>"
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
RAW_VCF="${OUTPUT_DIR}/${BIOSAMPLE}_gatk.raw.vcf.gz"
VCF="${OUTPUT_DIR}/${BIOSAMPLE}_gatk.vcf.gz"

# ================== SKIP IF RESULTS ALREADY EXIST ==================
if [[ -f "$VCF" && -f "${VCF}.tbi" ]]; then
    echo "[SKIP] GATK VCF output already exists for ${BIOSAMPLE}:"
    echo "       - $(basename "$VCF")"
    exit 0
fi

# ================== FIND BAM ==================
BAM_FILE="${BWA_DIR}/${BIOSAMPLE}.bam"

if [[ ! -f "$BAM_FILE" ]]; then
    echo "[ERROR] BAM file not found: ${BAM_FILE}"
    exit 1
fi

echo "[RUN] GATK HaplotypeCaller (VCF) for ${BIOSAMPLE}"
echo "[INFO] Threads: ${THREADS}"
echo "---------------------------------------------"

# ============================================================
# VCF MODE (MNPs)
# ============================================================
LOG_VCF="${OUTPUT_DIR}/${BIOSAMPLE}_gatk_vcf.log"

gatk HaplotypeCaller \
    -R "$REF" \
    -I "$BAM_FILE" \
    -O "$RAW_VCF" \
    --max-mnp-distance 1 \
    --native-pair-hmm-threads "$THREADS" \
    -A FisherStrand \
    -A StrandOddsRatio \
    -A QualByDepth \
    -A MappingQualityRankSumTest \
    -A ReadPosRankSumTest \
    -A DepthPerAlleleBySample \
    -A Coverage \
    2>&1 | tee "$LOG_VCF"

gatk VariantFiltration \
    -R "$REF" \
    -V "$RAW_VCF" \
    -O "$VCF" \
    --filter-expression "vc.isSNP() && QD < 2.0" \
    --filter-name "QD2" \
    --filter-expression "vc.isSNP() && QUAL < 30.0" \
    --filter-name "QUAL30" \
    --filter-expression "vc.isSNP() && SOR > 3.0" \
    --filter-name "SOR3" \
    --filter-expression "vc.isSNP() && FS > 60.0" \
    --filter-name "FS60" \
    --filter-expression "vc.isSNP() && MQ < 40.0" \
    --filter-name "MQ40" \
    --filter-expression "vc.isSNP() && MQRankSum < -12.5" \
    --filter-name "MQRankSum" \
    --filter-expression "vc.isSNP() && ReadPosRankSum < -8.0" \
    --filter-name "ReadPosRankSum" \
    --filter-expression "vc.isIndel() && QUAL < 30.0" \
    --filter-name "QUAL30" \
    --filter-expression "vc.isIndel() && FS > 200.0" \
    --filter-name "FS200" \
    --filter-expression "vc.isIndel() && ReadPosRankSum < -20.0" \
    --filter-name "ReadPosRankSum" \
    2>&1 | tee -a "$LOG_VCF"
    
# ================== CLEAN RAW VCF ==================
if [[ -f "$VCF" && -f "${VCF}.tbi" ]]; then
    rm -f "$RAW_VCF" "${RAW_VCF}.tbi"
    echo "[CLEAN] Removed raw VCF: $(basename "$RAW_VCF")"
fi

# ================== COMPLETION ==================
ELAPSED=$(( SECONDS - START_TIME ))
printf "[DONE] GATK VCF completed for %s (%02d min %02d sec)\n" \
    "$BIOSAMPLE" $((ELAPSED/60)) $((ELAPSED%60))

echo "[OUT] VCF: $VCF"
