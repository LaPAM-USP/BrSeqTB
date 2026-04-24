#!/usr/bin/env bash

set -e

echo "========================================="
echo "        BrSeqTB Installer"
echo "========================================="

# Detecta diretório real da pipeline
PIPELINE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BIN_DIR="$PIPELINE_DIR/bin"

echo ""
echo "Pipeline directory detected:"
echo "  $PIPELINE_DIR"
echo ""

# --------------------------------------------------
# Check Nextflow
# --------------------------------------------------
if ! command -v nextflow &> /dev/null; then
    echo "[ERROR] Nextflow not found in PATH."
    echo "Please install Nextflow before using BrSeqTB."
    exit 1
fi

echo "[OK] Nextflow detected"

# --------------------------------------------------
# Check Java
# --------------------------------------------------
if ! command -v java &> /dev/null; then
    echo "[ERROR] Java not found in PATH."
    echo "Please install Java (>=11) before using BrSeqTB."
    exit 1
fi

echo "[OK] Java detected"

# --------------------------------------------------
# Make wrapper executable
# --------------------------------------------------
chmod +x "$BIN_DIR/brseqtb"

# --------------------------------------------------
# Add to PATH safely (dynamic path)
# --------------------------------------------------
if grep -Fq "$BIN_DIR" ~/.bashrc; then
    echo "[INFO] BrSeqTB already present in PATH."
else
    echo ""
    echo "Adding BrSeqTB to your PATH..."
    echo "" >> ~/.bashrc
    echo "# BrSeqTB" >> ~/.bashrc
    echo "export PATH=\"$BIN_DIR:\$PATH\"" >> ~/.bashrc
    echo "[OK] Added to ~/.bashrc"
fi

echo ""
echo "========================================="
echo " Installation complete!"
echo "========================================="
echo ""
echo "Please open a new terminal or run:"
echo "  source ~/.bashrc"
echo ""
echo "Then you can run:"
echo "  brseqtb"
echo ""
