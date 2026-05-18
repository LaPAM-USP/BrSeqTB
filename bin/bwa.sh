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

TMP_DIR="${OUTPUT_DIR}/tmp"
RAW_BAM="${TMP_DIR}/${BIOSAMPLE}.raw.bam"
SORTED_BAM="${TMP_DIR}/${BIOSAMPLE}.sorted.bam"
FINAL_BAM="${OUTPUT_DIR}/${BIOSAMPLE}.bam"
METRICS_FILE="${TMP_DIR}/${BIOSAMPLE}_dupMetrics.txt"
FLAGSTAT_FILE="${TMP_DIR}/${BIOSAMPLE}_flagstat.txt"

# ===================== SKIP IF ALREADY DONE =====================
if [[ -s "$SUMMARY_CSV" && -f "$FINAL_BAM" ]]; then
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

# ===================== LOCATE READS =====================
mapfile -t R1_FILES < <(find "$TRIM_DIR" -type f -name "*_R1_001.fastq.gz" | sort)

if [[ "${#R1_FILES[@]}" -eq 0 ]]; then
    echo "[ERROR] Could not find any R1 FASTQs in: $TRIM_DIR"
    exit 1
fi

mkdir -p "$TMP_DIR"

echo "[INFO] Found ${#R1_FILES[@]} R1 FASTQ file(s) for ${BIOSAMPLE}"
echo "[RUN] Mapping ${BIOSAMPLE} with ${THREADS} threads..."

RG_BAMS=()

# ===================== MAPPING and ADD READ GROUPS =====================
for R1_FILE in "${R1_FILES[@]}"; do
    R2_FILE="${R1_FILE/_R1_/_R2_}"

    if [[ ! -f "$R2_FILE" ]]; then
        echo "[ERROR] Could not find matching R2 for: $R1_FILE"
        exit 1
    fi

    R1_BASENAME="$(basename "$R1_FILE")"

    if [[ "$R1_BASENAME" =~ ^(.+)_S([0-9]+)_L([0-9]{3})_R1_001\.fastq\.gz$ ]]; then
        SAMPLE_NUMBER="${BASH_REMATCH[2]}"
        LANE="${BASH_REMATCH[3]}"
    else
        echo "[ERROR] Could not extract sample number and lane from filename: $R1_BASENAME"
        exit 1
    fi

    READ_GROUP_ID="${BIOSAMPLE}_S${SAMPLE_NUMBER}_L${LANE}"
    READ_GROUP_PU="S${SAMPLE_NUMBER}_L${LANE}"

    RG_BAM="${TMP_DIR}/${READ_GROUP_ID}.raw.bam"
    RG_BAMS+=("$RG_BAM")

    echo "[RUN] Mapping ${R1_BASENAME}"
    echo "[INFO] Read group ID: ${READ_GROUP_ID}"
    echo "[INFO] Sample name: ${BIOSAMPLE}"
    echo "[INFO] Platform unit: ${READ_GROUP_PU}"

    bwa mem -t "$THREADS" \
      -R "@RG\tID:${READ_GROUP_ID}\tSM:${BIOSAMPLE}\tPL:ILLUMINA\tLB:${BIOSAMPLE}\tPU:${READ_GROUP_PU}" \
      "$REF" "$R1_FILE" "$R2_FILE" \
    | samtools view -b - > "$RG_BAM"
done

# ===================== MERGE RAW BAMS FROM ALL READ PAIRS / LANES =====================
if [[ "${#RG_BAMS[@]}" -eq 1 ]]; then
    mv "${RG_BAMS[0]}" "$RAW_BAM"
else
    echo "[RUN] Merging ${#RG_BAMS[@]} BAM files into one raw BAM for ${BIOSAMPLE}"

    samtools merge \
        -@ "$THREADS" \
        -f \
        "$RAW_BAM" \
        "${RG_BAMS[@]}"

    rm -f "${RG_BAMS[@]}"
fi

# ===================== SORT =====================
samtools sort -@ "$THREADS" "$RAW_BAM" -o "$SORTED_BAM"
rm -f "$RAW_BAM"

samtools index "$SORTED_BAM"

# ===================== MARK DUPLICATES =====================
picard MarkDuplicates \
    I="$SORTED_BAM" \
    O="$FINAL_BAM" \
    M="$METRICS_FILE" \
    VALIDATION_STRINGENCY=LENIENT

samtools index "$FINAL_BAM"

# ===================== FLAGSTAT =====================
samtools flagstat "$FINAL_BAM" > "$FLAGSTAT_FILE"

# ===================== HEADER OF CSV =====================
echo "biosample,filename,total_reads,mapped_reads,mapped_pct,duplicates_pct,properly_paired_pct,coverage_pct,status" > "$SUMMARY_CSV"

# ===================== METRICS =====================
total_reads=$(awk '$0 ~ / in total / {print $1; exit}' "$FLAGSTAT_FILE")

mapped_reads=$(awk '$0 ~ / primary mapped / {print $1; exit}' "$FLAGSTAT_FILE")
mapped_pct=$(awk -F'[()%]' '$0 ~ / primary mapped / {print $2; exit}' "$FLAGSTAT_FILE")

paired_pct=$(awk -F'[()%]' '$0 ~ / properly paired / {print $2; exit}' "$FLAGSTAT_FILE")

primary_reads=$(awk '$0 ~ / primary$/ {print $1; exit}' "$FLAGSTAT_FILE")
primary_duplicates=$(awk '$0 ~ / primary duplicates$/ {print $1; exit}' "$FLAGSTAT_FILE")

dup_total=$(awk -v dup="$primary_duplicates" -v total="$primary_reads" 'BEGIN {
    if (total > 0) printf "%.2f%%", (dup / total) * 100;
    else print "0.00%"
}')

coverage_pct=$(samtools coverage "$FINAL_BAM" 2>/dev/null | awk 'NR==2 {print $6; exit}')
[[ -z "$coverage_pct" ]] && coverage_pct="0.00"

mapped_value=$(echo "$mapped_pct" | tr -d '%')
coverage_value=$(echo "$coverage_pct" | tr -d '%')

status="PASS"
if (( $(echo "$mapped_value < $MIN_MAPPED" | bc -l) )) || \
   (( $(echo "$coverage_value < $MIN_COVERAGE" | bc -l) )); then
    status="FAIL"
fi

echo "${BIOSAMPLE},${BIOSAMPLE}.bam,${total_reads},${mapped_reads},${mapped_pct}%,${dup_total},${paired_pct}%,${coverage_pct},${status}" >> "$SUMMARY_CSV"
rm -rf "$TMP_DIR"

echo "[DONE] BWA finished for ${BIOSAMPLE} → ${status}"
echo "[OUT] ${FINAL_BAM}"
echo "[OUT] ${SUMMARY_CSV}"
