#!/usr/bin/env bash
# ============================================================
# Prepare Kaiju database (download, extract, verify)
#
# Usage:
#   ./kaijudb.sh               -> automatic mode (download if needed)
#   ./kaijudb.sh true          -> manual mode (use existing archive)
#
# Expected archive location (manual or automatic):
#   database/kaiju/kaiju_mycobacterium.tar.gz
#
# Output:
#   database/kaiju/db/
# ============================================================

export LC_ALL=C
set -euo pipefail

ADD_MANUAL="${1:-false}"
DB_DIR="database/kaiju/db"

DB_TAR="database/kaiju/db.tar.gz"
ZENODO_URL="https://zenodo.org/records/18064127/files/db.tar.gz"
EXPECTED_SHA256="74b05e77a5b43a4d0e6c81cc1dfe826596458889fec3b874ecd2a27a0a36eab6"


mkdir -p "$(dirname "$DB_TAR")"

# VERIFY DATABASE
verify_database() {
    echo "[DB] Verifying Kaiju database integrity..."

    local required=("nodes.dmp" "names.dmp")
    local fmi_file

    # Require exactly one .fmi file
    mapfile -t fmis < <(find "$DB_DIR" -maxdepth 1 -name "*.fmi" 2>/dev/null || true)

    if [[ "${#fmis[@]}" -ne 1 ]]; then
        echo "[DB] Expected exactly 1 .fmi file, found ${#fmis[@]}"
        return 1
    fi

    fmi_file="${fmis[0]}"

    for f in "${required[@]}"; do
        [[ -s "${DB_DIR}/${f}" ]] || {
            echo "[DB] Missing ${f}"
            return 1
        }
    done

    [[ -s "$fmi_file" ]] || {
        echo "[DB] Missing FMI index"
        return 1
    }

    # Basic corruption check
    head -c 8 "$fmi_file" >/dev/null 2>&1 || {
        echo "[DB] FMI appears corrupted"
        return 1
    }

    # Size sanity check (detect truncated files)
    if [[ "$(stat -c%s "$fmi_file")" -lt 1000000 ]]; then
        echo "[DB] FMI file size too small — possibly corrupted"
        return 1
    fi

    echo "[DB] Database OK."
    return 0
}

# VERIFY ARCHIVE 
verify_archive() {
    [[ -f "$DB_TAR" && -s "$DB_TAR" ]] || {
        echo "[DB] Archive missing: $DB_TAR"
        return 1
    }

    echo "[DB] Verifying archive checksum..."

    local sha
    sha=$(sha256sum "$DB_TAR" | awk '{print $1}')

    [[ "$sha" == "$EXPECTED_SHA256" ]] || {
        echo "[DB] Checksum mismatch"
        echo "     Expected: $EXPECTED_SHA256"
        echo "     Found:    $sha"
        return 1
    }

    echo "[DB] Archive checksum OK."
    return 0
}

# DOWNLOAD DATABASE
download_database() {
    echo "[DB] Downloading Kaiju DB from Zenodo..."

    for i in {1..3}; do
        echo "[DB] Attempt $i/3"
        if command -v wget >/dev/null; then
            wget -O "$DB_TAR" "$ZENODO_URL" && return 0
        else
            curl -L -o "$DB_TAR" "$ZENODO_URL" && return 0
        fi
        sleep 10
    done

    echo "[DB] Download failed after multiple attempts."
    return 1
}

# EXTRACT DATABASE
extract_database() {
    echo "[DB] Extracting Kaiju DB..."

    mkdir -p "$DB_DIR"
    rm -rf "${DB_DIR:?}"/*

    tar --strip-components=1 -xzf "$DB_TAR" -C "$DB_DIR"

    echo "[DB] Extraction complete."
}


# MAIN
echo "[DB] Checking Kaiju database..."

if verify_database; then
    echo "[DB] Database already prepared. Skipping."
    exit 0
fi

if [[ "$ADD_MANUAL" == "true" ]]; then
    echo "[DB] Manual mode enabled."

    if ! verify_archive; then
        echo "[DB] Manual archive not found or invalid."
        echo "     Expected at: $DB_TAR"
        exit 1
    fi
else
    if ! verify_archive; then
        rm -f "$DB_TAR"
        download_database
        verify_archive || {
            echo "[DB] Archive verification failed after download."
            exit 1
        }
    fi
fi

extract_database
verify_database || {
    echo "[DB] Database preparation failed."
    exit 1
}

echo "[DB] Kaiju database ready in: $DB_DIR"


