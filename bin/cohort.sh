#!/usr/bin/env bash
# ============================================================
# GATK Cohort Variant Calling and Hard Filtering (Pipeline-safe)
#
# Usage:
#   ./cohort.sh <manifest.tsv> [--demo]
#
# Behavior:
#   - DEFAULT: usa apenas biosamples do manifest.tsv
#   - --demo  : adiciona gVCFs de assets/demo/gatk/
#
# Notes:
#   - Reprodutível (usa Conda environment ativado pelo Nextflow)
#   - Seguro para execução via Nextflow (work/)
# ============================================================

export LC_ALL=C
set -euo pipefail

# ============================================================
# PROJECT
# ============================================================
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ============================================================
# ARGUMENT PARSING
# ============================================================
if [[ $# -lt 1 ]]; then
    echo "Usage: ./cohort.sh <manifest.tsv> [--demo]"
    exit 1
fi

MANIFEST="$1"
shift

USE_DEMO=false
for arg in "$@"; do
    case "$arg" in
        --demo)
            USE_DEMO=true
            echo "[OPT] Running WITH demo gVCFs."
            ;;
    esac
done

# ============================================================
# PATHS
# ============================================================
OUTPUT_DIR="${PROJECT_DIR}/cohort"
DB_PATH="${OUTPUT_DIR}/dbgatk"
BED="${OUTPUT_DIR}/variant_intervals.bed"

SNPS_FILTERED="${OUTPUT_DIR}/cohort_snps_filtered.vcf.gz"
INDELS_FILTERED="${OUTPUT_DIR}/cohort_indels_filtered.vcf.gz"

REF="${PROJECT_DIR}/database/mtbRef/NC0009623.fasta"
GATK_DIR="${PROJECT_DIR}/gatk"
DEMO_DIR="${PROJECT_DIR}/assets/demo/gatk"

THREADS="${NXF_TASK_CPUS:-1}"

mkdir -p "${OUTPUT_DIR}"

# ============================================================
# DEPENDENCY CHECKS
# ============================================================
for cmd in gatk; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found in Conda environment: $cmd"
        exit 1
    fi
done

# ============================================================
# INPUT CHECKS
# ============================================================
if [[ ! -f "$MANIFEST" ]]; then
    echo "[ERROR] Manifest not found: $MANIFEST"
    exit 1
fi

if [[ ! -f "$REF" ]]; then
    echo "[ERROR] Reference genome not found: $REF"
    exit 1
fi

# ============================================================
# BED (WHOLE GENOME — TB)
# ============================================================
if [[ ! -f "$BED" ]]; then
    echo -e "NC_000962.3\t1\t4411532" > "$BED"
fi

# ============================================================
# COLLECT BIOSAMPLES FROM MANIFEST
# ============================================================
echo "[INFO] Reading biosamples from manifest..."
mapfile -t BIOSAMPLES < <(
    awk -F'\t' 'NR>1 {print $1}' "$MANIFEST" | sort -u
)

if [[ ${#BIOSAMPLES[@]} -eq 0 ]]; then
    echo "[ERROR] No biosamples found in manifest."
    exit 1
fi

# ============================================================
# COLLECT GVCFs (STRICT FROM MANIFEST)
# ============================================================
GVCF_FILES=()

for b in "${BIOSAMPLES[@]}"; do
    gvcf="${GATK_DIR}/${b}/${b}.g.vcf.gz"
    if [[ ! -f "$gvcf" ]]; then
        echo "[ERROR] Missing gVCF for biosample: $b"
        exit 1
    fi
    GVCF_FILES+=("$gvcf")
done

# ============================================================
# OPTIONAL DEMO GVCFs
# ============================================================
if [[ "$USE_DEMO" == true && -d "$DEMO_DIR" ]]; then
    echo "[INFO] Adding DEMO gVCFs from assets/demo/"
    mapfile -t DEMO_GVCF < <(find "$DEMO_DIR" -type f -name "*.g.vcf.gz" | sort)
    GVCF_FILES+=("${DEMO_GVCF[@]}")
fi

echo "[INFO] Total gVCFs in cohort: ${#GVCF_FILES[@]}"

# ============================================================
# GENOMICSDB
# ============================================================
rm -rf "$DB_PATH"

echo "[RUN] GenomicsDBImport..."
gatk GenomicsDBImport \
    --genomicsdb-workspace-path "$DB_PATH" \
    --intervals "$BED" \
    $(for v in "${GVCF_FILES[@]}"; do echo "-V $v"; done) \
    --reader-threads "$THREADS"

# ============================================================
# GENOTYPE
# ============================================================
COHORT_RAW="${OUTPUT_DIR}/cohort_raw.vcf.gz"

gatk GenotypeGVCFs \
    -R "$REF" \
    -V "gendb://$DB_PATH" \
    -O "$COHORT_RAW"

# ============================================================
# SPLIT VARIANTS
# ============================================================
COHORT_SNPS="${OUTPUT_DIR}/cohort_snps.vcf.gz"
COHORT_INDELS="${OUTPUT_DIR}/cohort_indels.vcf.gz"

gatk SelectVariants -V "$COHORT_RAW" --select-type-to-include SNP   -O "$COHORT_SNPS"
gatk SelectVariants -V "$COHORT_RAW" --select-type-to-include INDEL -O "$COHORT_INDELS"

# ============================================================
# HARD FILTERING
# ============================================================
gatk VariantFiltration \
    -V "$COHORT_SNPS" \
    -filter "QD < 2.0"   --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 60.0"  --filter-name "FS60" \
    -O "$SNPS_FILTERED"

gatk VariantFiltration \
    -V "$COHORT_INDELS" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0"  --filter-name "FS200" \
    -O "$INDELS_FILTERED"

# ============================================================
# CLEANUP
# ============================================================
rm -f "$COHORT_RAW" "$COHORT_SNPS" "$COHORT_INDELS"
find "$OUTPUT_DIR" -type f -name "*.tbi" -delete

echo "[DONE] Cohort completed successfully."
echo "[OUT] ${OUTPUT_DIR}"
