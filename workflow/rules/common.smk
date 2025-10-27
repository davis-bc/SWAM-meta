import os
import glob
import pandas as pd
import re
import csv


# ----------------------------------------------
#   Build a mapping from sample R1 and R2 files
# ----------------------------------------------

# Load config variables
input_dir = config["in_dir"]
output_dir = config["out_dir"]

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