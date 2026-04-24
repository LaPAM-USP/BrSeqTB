#!/usr/bin/env python3
# =============================================================
# SNP cluster detection from SNP FASTA (output of snpMatrix.py)
#
# Method:
#   Distance = number of SNP mismatches
#   Cluster = connected components of MST filtered by cutoff
# =============================================================

import pandas as pd
import networkx as nx
from Bio import SeqIO
from pathlib import Path
import os

# =============================================================
# PROJECT STRUCTURE
# =============================================================
PROJECT_DIR = Path(
    os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
)

SNP_MATRIX_DIR   = PROJECT_DIR / "snpMatrix"
TRANSMISSION_DIR = PROJECT_DIR / "transmission"
RESULTS_DIR      = PROJECT_DIR / "results"

OUT_FASTA_FILE       = SNP_MATRIX_DIR / "snpmatrix.fasta"
OUT_DIST_MATRIX_TSV  = SNP_MATRIX_DIR / "snp_matrix.tsv"
OUT_CLUSTER_CSV      = TRANSMISSION_DIR / "transmission_clusters.csv"
OUT_CLUSTER_XLSX     = TRANSMISSION_DIR / "snp_distance_matrix.xlsx"
OUT_RESULTS_CSV      = RESULTS_DIR / "transmission_clusters.csv"

# Ensure output directories exist
TRANSMISSION_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================
# PARAMETERS
# =============================================================
CUTOFF = 12

# =============================================================
# 1. LOAD FASTA
# =============================================================
if not OUT_FASTA_FILE.exists():
    raise FileNotFoundError(f"Missing SNP FASTA: {OUT_FASTA_FILE}")

seqs = list(SeqIO.parse(OUT_FASTA_FILE, "fasta"))
samples = [s.id for s in seqs]
seq_dict = {s.id: str(s.seq) for s in seqs}

# Check consistency
lengths = set(len(s) for s in seq_dict.values())
if len(lengths) != 1:
    raise ValueError("All FASTA sequences must be same length.")

# =============================================================
# 2. COMPUTE PAIRWISE DISTANCES
# =============================================================
def snp_distance(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)

n = len(samples)
dist_mat = pd.DataFrame(0, index=samples, columns=samples, dtype=int)

for i in range(n):
    for j in range(i + 1, n):
        d = snp_distance(seq_dict[samples[i]], seq_dict[samples[j]])
        dist_mat.loc[samples[i], samples[j]] = d
        dist_mat.loc[samples[j], samples[i]] = d

# Save distance matrix (TSV)
dist_mat.to_csv(OUT_DIST_MATRIX_TSV, sep="\t")
print(f"[OK] Distance matrix written to {OUT_DIST_MATRIX_TSV}")

# Save distance matrix (Excel)
dist_mat.to_excel(OUT_CLUSTER_XLSX)
print(f"[OK] Excel SNP distance matrix saved to {OUT_CLUSTER_XLSX}")

# =============================================================
# 3. BUILD MST AND CLUSTERS
# =============================================================
G = nx.Graph()
for i in samples:
    for j in samples:
        if i != j:
            G.add_edge(i, j, weight=int(dist_mat.loc[i, j]))

mst = nx.minimum_spanning_tree(G, weight="weight")

# Keep edges <= cutoff
edges_short = [
    (u, v, d)
    for u, v, d in mst.edges(data=True)
    if int(d["weight"]) <= CUTOFF
]

G_cut = nx.Graph()
G_cut.add_edges_from(edges_short)

clusters = list(nx.connected_components(G_cut))

# Map sample -> cluster_id
cluster_id = {}
for cid, comp in enumerate(clusters, 1):
    for s in comp:
        cluster_id[s] = cid

# Isolated samples
for s in samples:
    if s not in cluster_id:
        cluster_id[s] = 0

# =============================================================
# 4. OUTPUT TRANSMISSION CLUSTERS
# =============================================================
df_out = pd.DataFrame({
    "sample": samples,
    "cluster_id": [cluster_id[s] for s in samples]
})

df_out.to_csv(OUT_CLUSTER_CSV, index=False)
print(f"[OK] Transmission clusters saved to {OUT_CLUSTER_CSV}")

df_out.to_csv(OUT_RESULTS_CSV, index=False)
print(f"[OK] Transmission clusters also saved to {OUT_RESULTS_CSV}")

print("[DONE] Transmission clustering complete.")

