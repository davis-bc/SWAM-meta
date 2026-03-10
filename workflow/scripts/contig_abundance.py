"""
contig_abundance.py

Snakemake script (called via script: directive).

Maps clean paired-end reads to assembled contigs with minimap2, computes
per-contig mean depth with samtools depth, then normalises to copies-per-genome
(cpg) using the SCG DIAMOND alignment produced by the short_reads rule.

cpg formula (same as short_reads_processing.R):
    n_genomes = sum(alignment_length / gene_length) / 40
    cpg       = mean_depth / n_genomes
"""

import os
import subprocess
import pandas as pd

# Snakemake-injected variables
r1      = snakemake.input.r1
r2      = snakemake.input.r2
contigs = snakemake.input.contigs
scgs    = str(snakemake.input.scgs) if snakemake.input.get("scgs") else ""
bam     = snakemake.output.bam
bai     = snakemake.output.bai
tsv     = snakemake.output.tsv
threads = snakemake.resources.threads

os.makedirs(os.path.dirname(bam), exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Map reads to all contigs
# ---------------------------------------------------------------------------
map_cmd = (
    f"minimap2 -t {threads} -ax sr {contigs} {r1} {r2} "
    f"| samtools sort -@ 4 -o {bam} - "
)
subprocess.run(map_cmd, shell=True, check=True, executable="/bin/bash")
subprocess.run(f"samtools index {bam}", shell=True, check=True)

# ---------------------------------------------------------------------------
# 2. Compute per-contig mean depth
# ---------------------------------------------------------------------------
depth_proc = subprocess.run(
    f"samtools depth -a {bam}",
    shell=True, capture_output=True, text=True, check=True
)

# samtools depth output: contig  pos  depth  (one line per base, no header)
depth_data: dict[str, list] = {}
for line in depth_proc.stdout.splitlines():
    parts = line.split("\t")
    if len(parts) < 3:
        continue
    contig_id = parts[0]
    d = float(parts[2])
    if contig_id not in depth_data:
        depth_data[contig_id] = []
    depth_data[contig_id].append(d)

mean_depth = {cid: sum(vals) / len(vals) for cid, vals in depth_data.items()}

# ---------------------------------------------------------------------------
# 3. Estimate n_genomes from SCG alignments
#    Falls back to n_genomes=1 when scgs is unavailable (run_short_reads=False)
# ---------------------------------------------------------------------------
if scgs and os.path.exists(scgs) and os.path.getsize(scgs) > 0:
    scg_df = pd.read_csv(
        scgs, sep="\t", header=None,
        names=["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "slen"]
    )
    scg_df = scg_df.drop_duplicates("qseqid")
    if len(scg_df) > 0:
        n_genomes = (scg_df["length"] / scg_df["slen"]).sum() / 40
    else:
        n_genomes = 1.0
else:
    n_genomes = 1.0  # fallback: no SCG data available

n_genomes = max(n_genomes, 1e-6)

# ---------------------------------------------------------------------------
# 4. Write output TSV
# ---------------------------------------------------------------------------
rows = [
    {"contig_id": cid, "mean_depth": d, "abundance_cpg": d / n_genomes}
    for cid, d in mean_depth.items()
]
out_df = pd.DataFrame(rows, columns=["contig_id", "mean_depth", "abundance_cpg"])
out_df.to_csv(tsv, sep="\t", index=False)
