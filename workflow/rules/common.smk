import os
import glob
import pandas as pd
import re
import csv

# ---------------------------------------------------------------------------
#   Test-mode path overrides
#   Activated with:  snakemake --config test=True
# ---------------------------------------------------------------------------

_TEST = config.get("test", False)
_REPO = os.path.dirname(os.path.dirname(os.path.dirname(workflow.snakefile)))

if _TEST:
    input_dir  = os.path.join(_REPO, "test", "data")
    output_dir = os.path.join(_REPO, "test", "output")
    # Override all database paths to mini test versions
    config["scg_db"]      = os.path.join(_REPO, "test", "dbs", "scg", "SCGs_40.fasta")
    config["uniref50_db"]  = os.path.join(_REPO, "test", "dbs", "uniref50", "uniref50_mmseqs")
    config["markers_db"]   = os.path.join(_REPO, "test", "dbs", "markers")
    # GTDB-tk and METABOLIC are skipped in test mode (see Snakefile rule all)
    config["skip_gtdbtk"]    = True
    config["skip_metabolic"] = True
    config["skip_checkm2"]   = True
    # Increase geNomad splits to reduce peak memory on low-RAM machines
    config.setdefault("genomad_splits", 4)
else:
    # Load config variables
    input_dir  = config.get("in_dir", "")
    output_dir = config.get("out_dir", "")

# ---------------------------------------------------------------------------
#   Resource helper: returns test value when running in test mode
# ---------------------------------------------------------------------------

def res(production, test=4000):
    """Return test-mode resource value when config test=True, else production."""
    return test if _TEST else production

# ---------------------------------------------------------------------------
#   Module flags — control which pipeline stages are active
# ---------------------------------------------------------------------------

_RUN_SR   = config.get("run_short_reads", True)
_RUN_CTG  = config.get("run_contigs",     True)
_RUN_MAGS = config.get("run_mags",        True)

# ---------------------------------------------------------------------------
#   Config validation for external-input scenarios
# ---------------------------------------------------------------------------

from snakemake.exceptions import WorkflowError

if not _RUN_SR and (_RUN_CTG or _RUN_MAGS) and not config.get("clean_reads_dir"):
    raise WorkflowError(
        "run_short_reads: False requires 'clean_reads_dir' in config.yaml "
        "(directory with pre-filtered {sample}*R1*.fastq* read pairs)."
    )

if not _RUN_CTG and _RUN_MAGS and not config.get("contigs_dir"):
    raise WorkflowError(
        "run_contigs: False with run_mags: True requires 'contigs_dir' in config.yaml "
        "(directory with pre-assembled {sample}.contigs.fa files)."
    )

# ---------------------------------------------------------------------------
#   Input helper functions for dynamic stage routing
# ---------------------------------------------------------------------------

def get_clean_r1(wildcards):
    """Return path to host-filtered R1 reads."""
    if _RUN_SR:
        return os.path.join(output_dir, "data", "clean_reads", f"{wildcards.sample}_R1.clean.fastq.gz")
    for pattern in [
        os.path.join(config["clean_reads_dir"], f"{wildcards.sample}*R1*.fastq*"),
        os.path.join(config["clean_reads_dir"], f"{wildcards.sample}*_1.fastq*"),
    ]:
        hits = glob.glob(pattern)
        if hits:
            return sorted(hits)[0]
    raise WorkflowError(
        f"No R1 reads found for sample '{wildcards.sample}' in clean_reads_dir: {config['clean_reads_dir']}"
    )

def get_clean_r2(wildcards):
    """Return path to host-filtered R2 reads."""
    if _RUN_SR:
        return os.path.join(output_dir, "data", "clean_reads", f"{wildcards.sample}_R2.clean.fastq.gz")
    for pattern in [
        os.path.join(config["clean_reads_dir"], f"{wildcards.sample}*R2*.fastq*"),
        os.path.join(config["clean_reads_dir"], f"{wildcards.sample}*_2.fastq*"),
    ]:
        hits = glob.glob(pattern)
        if hits:
            return sorted(hits)[0]
    raise WorkflowError(
        f"No R2 reads found for sample '{wildcards.sample}' in clean_reads_dir: {config['clean_reads_dir']}"
    )

def get_contigs(wildcards):
    """Return path to assembled contigs (MEGAHIT output or user-provided)."""
    if _RUN_CTG:
        return os.path.join(output_dir, "data", "megahit", f"{wildcards.sample}.contigs.fa")
    return os.path.join(config["contigs_dir"], f"{wildcards.sample}.contigs.fa")

# ----------------------------------------------
#   Build a mapping from sample R1 and R2 files
# ----------------------------------------------

# Helper: extract sample name from fastq filename
def extract_sample_name(filename):
    base = os.path.basename(filename)
    if base.endswith('.gz'):
        base = base[:-3]
    m = re.match(r"(.+?)(?:_R?1(?:_001)?|_1)(?:\.fastq)?$", base)
    if m:
        return m.group(1)
    else:
        return re.sub(r"(_R?1(?:_001)?|_1)(\.fastq)?(\.gz)?$", "", base)

# Find R1/R2 files
r1_files = glob.glob(os.path.join(input_dir, "*R1*.fastq*")) + glob.glob(os.path.join(input_dir, "*_1.fastq*"))
r2_files = glob.glob(os.path.join(input_dir, "*R2*.fastq*")) + glob.glob(os.path.join(input_dir, "*_2.fastq*"))

# Map sample names to paired files
sample_to_files = {}
for r1 in r1_files:
    sample = extract_sample_name(r1)
    possible_r2s = [r1.replace('R1', 'R2'), r1.replace('_1', '_2')]
    found_r2 = None
    for r2_candidate in possible_r2s:
        if os.path.exists(r2_candidate):
            found_r2 = r2_candidate
            break
    if not found_r2:
        for r2 in r2_files:
            if extract_sample_name(r2) == sample:
                found_r2 = r2
                break
    if found_r2:
        sample_to_files[sample] = {'r1': r1, 'r2': found_r2}

samples = sorted(sample_to_files.keys())

def get_r1(wildcards):
    return sample_to_files[wildcards.sample]['r1']

def get_r2(wildcards):
    return sample_to_files[wildcards.sample]['r2']