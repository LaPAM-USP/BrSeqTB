#!/usr/bin/env bash
# ============================================================
# Kaiju classification of paired-end reads
# Usage: ./kaiju.sh <biosample ID>
# Requires: database/kaiju/db already prepared
# ============================================================

export LC_ALL=C
set -euo pipefail

BIOSAMPLE="${1:-}"
THREADS="${NXF_TASK_CPUS:-1}"

# IMPORTANT:
# When executed by Nextflow, cwd = work/<hash>
# So all paths must be resolved from project root
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

TRIM_DIR="${PROJECT_DIR}/trimmomatic/${BIOSAMPLE}"
DB_DIR="${PROJECT_DIR}/database/kaiju/db"
OUTPUT_DIR="${PROJECT_DIR}/kaiju/${BIOSAMPLE}"

# ===================== CHECK INPUT =====================
[[ -z "$BIOSAMPLE" ]] && { echo "Usage: ./kaiju.sh <biosample ID>"; exit 1; }
[[ ! -d "$TRIM_DIR" ]] && { echo "[ERR] ${TRIM_DIR} not found."; exit 1; }

mkdir -p "$OUTPUT_DIR"

# ===================== SKIP IF ALREADY DONE =====================
SUMMARY_CSV="${OUTPUT_DIR}/${BIOSAMPLE}_kaiju_summary.csv"

if [[ -s "$SUMMARY_CSV" ]]; then
    echo "[SKIP] Kaiju already completed for ${BIOSAMPLE}"
    echo "[OUT] ${SUMMARY_CSV}"
    exit 0
fi

# ===================== CHECK DATABASE =====================
DB_FMI=$(find "$DB_DIR" -maxdepth 1 -name "*.fmi" | head -n 1)

[[ ! -s "${DB_DIR}/nodes.dmp" ]] && { echo "[ERR] nodes.dmp missing"; exit 1; }
[[ ! -s "${DB_DIR}/names.dmp" ]] && { echo "[ERR] names.dmp missing"; exit 1; }
[[ -z "$DB_FMI" ]] && { echo "[ERR] FMI index missing"; exit 1; }

# ===================== LOCATE READS =====================
R1=$(find "$TRIM_DIR" -type f -name "*_R1_001.fastq.gz" | head -n 1)
R2=$(find "$TRIM_DIR" -type f -name "*_R2_001.fastq.gz" | head -n 1)

[[ -z "$R1" || -z "$R2" ]] && { echo "[ERR] FASTQs not found."; exit 1; }

# ===================== RUN KAIJU =====================
OUTFILE="${OUTPUT_DIR}/${BIOSAMPLE}_kaiju.out"
LOGFILE="${OUTPUT_DIR}/${BIOSAMPLE}_kaiju.log"

echo "[INFO] Running Kaiju for ${BIOSAMPLE} with ${THREADS} threads"

kaiju \
    -t "${DB_DIR}/nodes.dmp" \
    -f "$DB_FMI" \
    -i "$R1" -j "$R2" \
    -o "$OUTFILE" \
    -z "$THREADS" -v \
    | tee "$LOGFILE"

# ===================== HTML REPORT =====================
KRONA_TXT="${OUTPUT_DIR}/${BIOSAMPLE}_kaiju_krona.txt"
HTML_REPORT="${OUTPUT_DIR}/${BIOSAMPLE}_kaiju_report.html"

if command -v kaiju2krona >/dev/null && command -v ktImportText >/dev/null; then
    kaiju2krona \
        -t "${DB_DIR}/nodes.dmp" \
        -n "${DB_DIR}/names.dmp" \
        -i "$OUTFILE" \
        -o "$KRONA_TXT"

    ktImportText \
        -o "$HTML_REPORT" \
        "$KRONA_TXT"
else
    echo "[WARN] kaiju2krona or ktImportText not found; HTML report not generated."
fi

# ===================== SUMMARY =====================
SUMMARY_CSV="${OUTPUT_DIR}/${BIOSAMPLE}_kaiju_summary.csv"

if command -v kaiju2table >/dev/null; then
    TMP_TSV=$(mktemp)

    kaiju2table \
        -t "${DB_DIR}/nodes.dmp" \
        -n "${DB_DIR}/names.dmp" \
        -r species \
        -o "$TMP_TSV" \
        "$OUTFILE"

    awk -v sample="$BIOSAMPLE" '
        BEGIN { FS="\t"; OFS="," }
        NR==1 { print "biosample","file","percent","reads","taxon_id","taxon_name"; next }
        { print sample,$1,$2,$3,$4,$5 }
    ' "$TMP_TSV" > "$SUMMARY_CSV"

    rm -f "$TMP_TSV"
fi

echo "[DONE] Kaiju finished for ${BIOSAMPLE}"
echo "[OUT] ${OUTFILE}"
echo "[OUT] ${HTML_REPORT}"
echo "[OUT] ${SUMMARY_CSV}"


