"""
contig_abundance.py

Snakemake script (called via script: directive).

Maps clean paired-end reads to assembled contigs with minimap2, computes
per-contig mean depth with samtools depth, then normalises to copies-per-genome
(cpg) using the SCG DIAMOND alignment produced by the short_reads rule.

cpg formula (matches short_reads_processing.R for directly comparable values):
    SCG_depth_i  = n_reads_mapped_to_SCG_i / (slen_i × 3)   [reads/base]
    median_SCG_depth = median(SCG_depth_i) across detected SCGs
    cpg = samtools_mean_depth / median_SCG_depth

Both numerator and denominator are in reads/base so the ratio is dimensionless
and directly comparable to the short-reads CPG from the KMA alignment.
Returns NaN when no SCG data is available.
"""

import math
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
# 3. Estimate median SCG depth (reads / base) for CPG normalisation.
#    Uses the same DIAMOND SCG alignment as short_reads_processing.R so that
#    contig CPG values are on the same scale as short-reads CPG values.
#
#    SCG_depth_i = n_reads_mapped_to_SCG_i / (slen_i × 3)   [reads/base]
#    median_SCG_depth = median across all detected SCGs (non-zero by construction)
#
#    CPG = samtools_mean_depth / median_SCG_depth
#    Both are in reads/base → ratio is directly comparable to short-reads CPG.
#
#    Falls back to NaN (propagated as empty in TSV) when no SCG hits exist.
# ---------------------------------------------------------------------------
median_scg_depth = float("nan")  # default: no SCG data

if scgs and os.path.exists(scgs) and os.path.getsize(scgs) > 0:
    scg_df = pd.read_csv(
        scgs, sep="\t", header=None,
        names=["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "slen"]
    )
    scg_df = scg_df.drop_duplicates("qseqid")  # one alignment per read (best hit)
    if len(scg_df) > 0:
        scg_per_gene = (
            scg_df.groupby("sseqid")
            .agg(scg_reads=("qseqid", "count"), slen=("slen", "first"))
            .reset_index()
        )
        scg_per_gene["scg_depth"] = scg_per_gene["scg_reads"] / (scg_per_gene["slen"] * 3)
        median_scg_depth = float(scg_per_gene["scg_depth"].median())

if median_scg_depth <= 0 or median_scg_depth != median_scg_depth:  # 0 or NaN
    median_scg_depth = float("nan")

# ---------------------------------------------------------------------------
# 4. Write output TSV
# ---------------------------------------------------------------------------
rows = [
    {
        "contig_id":     cid,
        "mean_depth":    d,
        "abundance_cpg": d / median_scg_depth if not math.isnan(median_scg_depth) else float("nan"),
    }
    for cid, d in mean_depth.items()
]
out_df = pd.DataFrame(rows, columns=["contig_id", "mean_depth", "abundance_cpg"])
out_df.to_csv(tsv, sep="\t", index=False)
