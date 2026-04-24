#!/usr/bin/env bash
# ============================================================
# Build custom SnpEff database for M. tuberculosis (H37Rv)
# Database ID: NC_0009623
#
# Uses ONLY:
#   - genes.gff
#   - sequences.fa
#
# Explicitly disables CDS / protein validation
# Safe to run multiple times (SKIP if already built)
#
# Compatible with Nextflow Conda/Mamba environment
# (uses snpEff from PATH)
#
# IMPORTANT:
#   - Does NOT modify Conda environment
#   - Stores custom DB inside project directory
# ============================================================

export LC_ALL=C
set -euo pipefail

START_TIME=$SECONDS
GENOME_DB="NC_0009623"

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

# -------------------- Locate SnpEff (from Conda env) --------------------
if ! command -v snpEff >/dev/null 2>&1; then
    echo "[ERROR] snpEff not found in PATH (Conda environment not active?)"
    exit 1
fi

SNPEFF_BIN="$(which snpEff)"

# -------------------- Custom SnpEff project location --------------------
# We NEVER write inside Conda environment.
# Custom database is stored inside project.
SNPEFF_PROJECT_DIR="${PROJECT_DIR}/database/snpeff"
SNPEFF_CONFIG="${SNPEFF_PROJECT_DIR}/snpEff.config"
SNPEFF_DATA_ROOT="${SNPEFF_PROJECT_DIR}/data"
SNPEFF_DATA_DIR="${SNPEFF_DATA_ROOT}/${GENOME_DB}"

# -------------------- Reference inputs --------------------
REF_DIR="${PROJECT_DIR}/database/mtbRef"
GFF="${REF_DIR}/genes.gff"
FASTA="${REF_DIR}/sequences.fa"

# -------------------- Checks --------------------
[[ -f "$GFF" ]]   || { echo "[ERROR] genes.gff not found: $GFF"; exit 1; }
[[ -f "$FASTA" ]] || { echo "[ERROR] sequences.fa not found: $FASTA"; exit 1; }

echo "[INFO] snpEff binary : $SNPEFF_BIN"
echo "[INFO] Project DB dir: $SNPEFF_PROJECT_DIR"
echo "[INFO] Data directory: $SNPEFF_DATA_DIR"
echo "---------------------------------------------"

# -------------------- Prepare project directories --------------------
mkdir -p "$SNPEFF_DATA_DIR"

# -------------------- Initialize config if missing --------------------
# Copy default config from Conda package only once
if [[ ! -f "$SNPEFF_CONFIG" ]]; then
    echo "[INFO] Creating local snpEff.config"

    DEFAULT_CONFIG=$(find "$(dirname "$SNPEFF_BIN")/.." -name snpEff.config 2>/dev/null | head -n 1)

    if [[ -z "$DEFAULT_CONFIG" ]]; then
        echo "[ERROR] Could not locate default snpEff.config inside Conda environment"
        exit 1
    fi

    cp "$DEFAULT_CONFIG" "$SNPEFF_CONFIG"
fi

# -------------------- SKIP IF ALREADY BUILT --------------------
BIN_FILE="${SNPEFF_DATA_DIR}/snpEffectPredictor.bin"

if [[ -f "$BIN_FILE" ]]; then
    echo "[SKIP] SnpEff database already exists for ${GENOME_DB}"
    echo "[OK]   ${BIN_FILE}"
    exit 0
fi

echo "[RUN] Building SnpEff database: ${GENOME_DB}"
echo "---------------------------------------------"

# -------------------- Copy reference files --------------------
cp -f "$GFF"   "${SNPEFF_DATA_DIR}/genes.gff"
cp -f "$FASTA" "${SNPEFF_DATA_DIR}/sequences.fa"

echo "[INFO] Copied genes.gff and sequences.fa"

# -------------------- Register genome in LOCAL config --------------------
if ! grep -q "^${GENOME_DB}\.genome" "$SNPEFF_CONFIG"; then
    echo "[INFO] Registering genome in snpEff.config"
    echo "${GENOME_DB}.genome : Mycobacterium_tuberculosis_H37Rv" >> "$SNPEFF_CONFIG"
else
    echo "[INFO] Genome already registered in snpEff.config"
fi

# -------------------- Build database --------------------
echo "[RUN] Executing snpEff build..."
echo "---------------------------------------------"

snpEff build \
    -c "$SNPEFF_CONFIG" \
    -dataDir "$SNPEFF_DATA_ROOT" \
    -gff3 \
    -noCheckCds \
    -noCheckProtein \
    -v "$GENOME_DB"

# -------------------- Validate --------------------
[[ -f "$BIN_FILE" ]] || { echo "[ERROR] Build failed: ${BIN_FILE} not found"; exit 1; }

ELAPSED=$(( SECONDS - START_TIME ))
printf "[DONE] SnpEff DB built for %s (%02d min %02d sec)\n" \
    "$GENOME_DB" $((ELAPSED/60)) $((ELAPSED%60))

echo "[OK]  ${BIN_FILE}"
