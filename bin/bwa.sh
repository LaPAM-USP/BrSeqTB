#!/usr/bin/env bash
# ============================================================
# Mapping paired-end reads using BWA-MEM and generating BAMs
# Compatible with Nextflow Conda/Mamba environment
# ============================================================

export LC_ALL=C
set -euo pipefail

BIOSAMPLE="${1:-}"
THREADS="${NXF_TASK_CPUS:-1}"

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

TRIM_DIR="${PROJECT_DIR}/trimmomatic/${BIOSAMPLE}"
REF="${PROJECT_DIR}/database/mtbRef/NC0009623.fasta"
OUTPUT_DIR="${PROJECT_DIR}/bwa/${BIOSAMPLE}"
SUMMARY_CSV="${OUTPUT_DIR}/${BIOSAMPLE}_bwa_summary.csv"

MIN_MAPPED=95
MIN_COVERAGE=90

# ===================== CHECK INPUT =====================
if [[ -z "$BIOSAMPLE" ]]; then
    echo "Usage: ./bwa.sh <biosample>"
    exit 1
fi

if [[ ! -d "$TRIM_DIR" ]]; then
    echo "[ERROR] Trimmed reads directory not found: $TRIM_DIR"
    exit 1
fi

if [[ ! -f "$REF" ]]; then
    echo "[ERROR] Reference genome not found: $REF"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# ===================== SKIP IF ALREADY DONE =====================
if [[ -s "$SUMMARY_CSV" ]]; then
    echo "[SKIP] BWA already completed for ${BIOSAMPLE}"
    echo "[OUT] ${SUMMARY_CSV}"
    exit 0
fi

# ===================== DEPENDENCY CHECKS =====================
for cmd in bwa samtools bc java picard; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: $cmd"
        exit 1
    fi
done

echo "[INFO] Using bwa: $(which bwa)"
echo "[INFO] Using samtools: $(which samtools)"
echo "[INFO] Using picard: $(which picard)"
echo "---------------------------------------------"

# ===================== HEADER OF CSV =====================
echo "biosample,filename,total_reads,mapped_reads,mapped_pct,duplicates_pct,properly_paired_pct,coverage_pct,status" > "$SUMMARY_CSV"

# ===================== LOCATE READS =====================
R1_FILE=$(find "$TRIM_DIR" -type f -name "*_R1_001.fastq.gz" | head -n 1)
R2_FILE="${R1_FILE/_R1_/_R2_}"

if [[ -z "$R1_FILE" || ! -f "$R2_FILE" ]]; then
    echo "[ERROR] Could not find paired R1/R2 FASTQs."
    exit 1
fi

TMP_DIR="${OUTPUT_DIR}/tmp"
mkdir -p "$TMP_DIR"

RAW_BAM="${TMP_DIR}/${BIOSAMPLE}.raw.bam"
SORTED_BAM="${TMP_DIR}/${BIOSAMPLE}.sorted.bam"
NAMED_BAM="${TMP_DIR}/${BIOSAMPLE}.named.bam"
FINAL_BAM="${OUTPUT_DIR}/${BIOSAMPLE}.bam"
METRICS_FILE="${TMP_DIR}/${BIOSAMPLE}_dupMetrics.txt"
FLAGSTAT_FILE="${TMP_DIR}/${BIOSAMPLE}_flagstat.txt"

echo "[RUN] Mapping ${BIOSAMPLE} with ${THREADS} threads..."

# ===================== MAPPING =====================
bwa mem -t "$THREADS" "$REF" "$R1_FILE" "$R2_FILE" \
| samtools view -bS - > "$RAW_BAM"

# ===================== SORT =====================
samtools sort -@ "$THREADS" "$RAW_BAM" -o "$SORTED_BAM"
rm -f "$RAW_BAM"

samtools index "$SORTED_BAM"

# ===================== ADD READ GROUPS =====================
picard AddOrReplaceReadGroups \
    I="$SORTED_BAM" \
    O="$NAMED_BAM" \
    RGID=1 \
    RGLB="$BIOSAMPLE" \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM="$BIOSAMPLE" \
    VALIDATION_STRINGENCY=LENIENT

# ===================== MARK DUPLICATES =====================
picard MarkDuplicates \
    I="$NAMED_BAM" \
    O="$FINAL_BAM" \
    M="$METRICS_FILE" \
    VALIDATION_STRINGENCY=LENIENT

samtools index "$FINAL_BAM"

# ===================== FLAGSTAT =====================
samtools flagstat "$FINAL_BAM" > "$FLAGSTAT_FILE"

# ===================== METRICS =====================
total_reads=$(grep "in total" "$FLAGSTAT_FILE" | awk '{print $1}')
mapped_reads=$(grep " mapped (" "$FLAGSTAT_FILE" | awk '{print $1}')
mapped_pct=$(grep " mapped (" "$FLAGSTAT_FILE" | sed 's/.*(\(.*%\).*/\1/')
paired_pct=$(grep "properly paired (" "$FLAGSTAT_FILE" | sed 's/.*(\(.*%\).*/\1/')
dup_pct=$(grep " duplicates" "$FLAGSTAT_FILE" | sed 's/.*(\(.*%\).*/\1/')

coverage_pct=$(samtools coverage "$FINAL_BAM" 2>/dev/null | awk 'NR==2 {print $6}')
[[ -z "$coverage_pct" ]] && coverage_pct="0.00"

mapped_value=$(echo "$mapped_pct" | tr -d '%')
coverage_value=$(echo "$coverage_pct" | tr -d '%')

status="PASS"
if (( $(echo "$mapped_value < $MIN_MAPPED" | bc -l) )) || \
   (( $(echo "$coverage_value < $MIN_COVERAGE" | bc -l) )); then
    status="FAIL"
fi

echo "${BIOSAMPLE},${BIOSAMPLE}.bam,${total_reads},${mapped_reads},${mapped_pct},${dup_pct},${paired_pct},${coverage_pct},${status}" >> "$SUMMARY_CSV"

rm -rf "$TMP_DIR"

echo "[DONE] BWA finished for ${BIOSAMPLE} → ${status}"
echo "[OUT] ${FINAL_BAM}"
echo "[OUT] ${SUMMARY_CSV}"
