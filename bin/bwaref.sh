#!/usr/bin/env bash
# ============================================================
# Prepare BWA reference (check ref + check index)
# Usage: ./bwaref.sh
# Requires: database/mtbRef/NC0009623.fasta
#
# Compatible with Nextflow Conda/Micromamba environment
# ============================================================

export LC_ALL=C
set -euo pipefail

REF_DIR="database/mtbRef"
REF="${REF_DIR}/NC0009623.fasta"


# CHECK REFERENCE

if [[ ! -f "$REF" ]]; then
    echo "[ERROR] Reference genome not found: ${REF}"
    exit 1
fi


# LOCATE BWA (FROM CONDA ENV)

if ! command -v bwa >/dev/null 2>&1; then
    echo "[ERROR] bwa not found in PATH (Conda environment not active?)"
    exit 1
fi

BWA_BIN="$(command -v bwa || true)"
echo "[INFO] Using bwa binary: ${BWA_BIN:-not found}"

bwa_version="$(bwa 2>&1 | head -n1 || true)"
echo "[INFO] bwa version: ${bwa_version:-unknown}"


# CHECK INDEX

INDEX_FILES=(
    "${REF}.bwt"
    "${REF}.pac"
    "${REF}.ann"
    "${REF}.amb"
    "${REF}.sa"
)

existing=0

for f in "${INDEX_FILES[@]}"; do
    if [[ -f "$f" ]]; then
        existing=$((existing + 1))
    fi
done

if [[ "$existing" -eq 5 ]]; then
    echo "[OK] Complete BWA index already present for ${REF}"

elif [[ "$existing" -gt 0 ]]; then
    echo "[ERROR] Partial BWA index detected for ${REF}"
    echo "[ERROR] Remove the following files and re-run:"
    for f in "${INDEX_FILES[@]}"; do
        [[ -f "$f" ]] && echo "  - $f"
    done
    exit 1

else
    echo "[INFO] Building BWA index..."
    bwa index "$REF"
    echo "[OK] BWA index created."
fi

echo "[DONE] Reference ready."
