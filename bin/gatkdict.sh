#!/usr/bin/env bash
# ============================================================
# Prepare reference for GATK (FAI + sequence dictionary)
# Usage: ./gatkdict.sh
# Requires: database/mtbRef/NC0009623.fasta
# ============================================================

export LC_ALL=C
set -euo pipefail

REF="database/mtbRef/NC0009623.fasta"
DICT="database/mtbRef/NC0009623.dict"

# CHECK REFERENCE
if [[ ! -f "$REF" ]]; then
    echo "[ERROR] Reference genome not found: ${REF}"
    exit 1
fi

# DEPENDENCIES
for cmd in gatk samtools; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: $cmd"
        exit 1
    fi
done

echo "[INFO] gatk version: $(gatk --version 2>&1 | head -n1)"
echo "[INFO] samtools version: $(samtools --version | head -n1)"

# FASTA INDEX
if [[ -f "${REF}.fai" && -s "${REF}.fai" ]]; then
    echo "[OK] FASTA index already exists."
else
    echo "[INFO] Creating FASTA index (.fai)..."
    rm -f "${REF}.fai"
    samtools faidx "$REF"
fi

# SEQUENCE DICTIONARY
if [[ -f "$DICT" && -s "$DICT" ]]; then
    # minimal sanity check
    grep -q "^@SQ" "$DICT" || {
        echo "[ERROR] Invalid sequence dictionary format."
        exit 1
    }
    echo "[OK] Sequence dictionary already exists."
else
    echo "[INFO] Creating sequence dictionary (.dict)..."
    rm -f "$DICT"
    gatk CreateSequenceDictionary \
        -R "$REF" \
        -O "$DICT"
fi

echo "[DONE] GATK reference preparation completed."

