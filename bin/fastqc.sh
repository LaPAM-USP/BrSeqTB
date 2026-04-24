#!/usr/bin/env bash
# ============================================================
# Run FastQC for all files belonging to a given biosample
# Usage: ./fastqc.sh <biosample ID>
# Input: reads/<biosample ID>_SX_L00X_RX_00X.fastq.gz
# Input files must use the Illumina read naming convention:
# 1827-22_S1_L001_R1_001.fastq.gz
# ││││││││││││││││││││││││││
# └─── biosample ID (1827-22)
#        └── biosample number (S1)
#           └── lane (L002)
#                └── read direction (R1/R2)
#                   └── file index (_001)
# Output: fastqc/<biosample ID>/
# ============================================================

export LC_ALL=C
set -euo pipefail

BIOSAMPLE="${1:-}"

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

READS_DIR="${PROJECT_DIR}/reads"
OUTPUT_DIR="${PROJECT_DIR}/fastqc/${BIOSAMPLE}"
REPORT_CSV="${OUTPUT_DIR}/${BIOSAMPLE}_fastqc_summary.csv"

MIN_QUALITY=20      # minimum mean per-base quality for PASS
GC_MIN=64.5         # lower GC threshold for M. tuberculosis
GC_MAX=66.5         # upper GC threshold for M. tuberculosis

# ================== CHECK INPUT ==================
if [[ -z "$BIOSAMPLE" ]]; then
    >&2 echo "Usage: ./fastqc.sh <biosample ID>"
    exit 1
fi

if [[ ! -d "$READS_DIR" ]]; then
    >&2 echo "[ERROR] Reads directory not found: $READS_DIR"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# ===================== SKIP IF ALREADY DONE =====================
if [[ -s "$REPORT_CSV" ]]; then
    >&2 echo "[SKIP] FastQC already completed for ${BIOSAMPLE}"
    >&2 echo "[OUT] ${REPORT_CSV}"
    exit 0
fi
# ===============================================================


>&2 echo "[RUN] Running FastQC for biosample: $BIOSAMPLE"
>&2 echo "[DIR] Searching files in: $READS_DIR"
>&2 echo "[OUT] Output directory: $OUTPUT_DIR"
>&2 echo "---------------------------------------------"

# ======== FIND FILES (anchored prefix) ========
# Strictly match Illumina-style names for this biosample (supports multiple lanes)
FILES=$(find "$READS_DIR" -maxdepth 1 -type f \
    \( -name "${BIOSAMPLE}_S*_L*_R*_001.fastq.gz" -o -name "${BIOSAMPLE}_S*_L*_R*_001.fastq" \) \
    | sort)

if [[ -z "${FILES}" ]]; then
    >&2 echo "[ERROR] No FASTQ files found for biosample '${BIOSAMPLE}' in '${READS_DIR}'."
    >&2 echo "[HINT] Expected e.g.: ${BIOSAMPLE}_S1_L001_R1_001.fastq.gz and ${BIOSAMPLE}_S1_L001_R2_001.fastq.gz"
    exit 1
fi

FILE_COUNT=$(echo "$FILES" | wc -l | tr -d ' ')

# ======== CHECK FILE COUNT IS EVEN ========
if (( FILE_COUNT % 2 != 0 )); then
    >&2 echo "[ERROR] Found an odd number of files (${FILE_COUNT}) for biosample '${BIOSAMPLE}'."
    >&2 echo "[HINT] Each biosample must have paired-end reads (R1 and R2 per lane)."
    >&2 echo "[INFO] Files found:"
    echo "$FILES" | sed 's/^/   • /' >&2
    exit 1
fi

>&2 echo "[INFO] Files found (${FILE_COUNT}):"
echo "$FILES" | sed 's/^/   • /' >&2
>&2 echo "---------------------------------------------"

# ======== RUN FASTQC ========
for FILE in $FILES; do
    >&2 echo "[RUN] FastQC on: $FILE"
    fastqc -o "$OUTPUT_DIR" "$FILE"
done

>&2 echo "[OK] FastQC completed for all files."
>&2 echo "[INFO] Unzipping FastQC reports..."
for ZIP in "${OUTPUT_DIR}"/*_fastqc.zip; do
    unzip -qo "$ZIP" -d "$OUTPUT_DIR"
done

>&2 echo "[INFO] Generating summary report: ${REPORT_CSV}"
>&2 echo "---------------------------------------------"

echo "biosample,filename,total_sequences,poor_quality_sequences,sequence_length,percent_GC,mean_per_base_quality,status" > "$REPORT_CSV"

for QC_DIR in "${OUTPUT_DIR}"/*_fastqc; do
    DATA_FILE="${QC_DIR}/fastqc_data.txt"
    [[ ! -f "$DATA_FILE" ]] && continue

    filename=$(grep "^Filename" "$DATA_FILE" | cut -f2)
    total=$(grep "^Total Sequences" "$DATA_FILE" | cut -f2)
    poor=$(grep "^Sequences flagged as poor quality" "$DATA_FILE" | cut -f2)
    length=$(grep "^Sequence length" "$DATA_FILE" | cut -f2)
    gc=$(grep "^%GC" "$DATA_FILE" | cut -f2)

    # Mean "Per base sequence quality"
    mean_quality=$(awk '
        /^#Base/{getline; while ($1 ~ /^[0-9]/) {sum+=$2; n++; getline} }
        END{if (n>0) printf "%.2f", sum/n; else print "NA"}' "$DATA_FILE")

    # Determine PASS/FAIL
    status="PASS"

    if [[ "$mean_quality" == "NA" || -z "$gc" ]]; then
        status="FAIL"
    else
        # Check mean quality
        if (( $(echo "$mean_quality < $MIN_QUALITY" | bc -l) )); then
            status="FAIL"
        fi
        # Check GC range
        if (( $(echo "$gc < $GC_MIN" | bc -l) )) || (( $(echo "$gc > $GC_MAX" | bc -l) )); then
            status="FAIL"
        fi
    fi

    echo "${BIOSAMPLE},${filename},${total},${poor},${length},${gc},${mean_quality},${status}" >> "$REPORT_CSV"
done

>&2 echo "[OK] Summary CSV generated: $REPORT_CSV"
>&2 echo "---------------------------------------------"
>&2 echo "[DONE] Analysis completed for biosample: $BIOSAMPLE"
>&2 echo "[OUT] Reports available in: $OUTPUT_DIR"

