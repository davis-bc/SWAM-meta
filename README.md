# SWAM-meta

Snakemake workflow for end-to-end metagenomic analysis of antibiotic resistance in environmental bacterial communities.

## Overview

SWAM-meta takes paired-end FASTQ files and produces:

- **Short-read summaries** — AMR gene abundance (copies per genome, cpg) and metagenomic coverage per sample
- **Assembled contigs** — de novo assembly, contig classification (chromosome / plasmid / phage), plasmid system inference
- **Contig-level annotation table** — taxonomy, AMR genes, mobile genetic elements (MGEs), and abundance in one TSV
- **Metagenome-Assembled Genomes (MAGs)** — binning, AMR/MGE annotation, taxonomy, metabolic potential, and relative abundance per MAG

## Pipeline stages

```
raw reads
  └─ fastp QC + minimap2 human-read removal
       ├─ KMA → AMR gene alignment (AMRFinderPlus DB)
       ├─ DIAMOND → 40 single-copy genes (genome count estimation)
       ├─ Nonpareil → metagenomic coverage
       └─ summary CSVs (fastp_summary.csv, short_reads_output.csv)

clean reads
  └─ MEGAHIT assembly (meta-large, min contig 1 000 bp)
       ├─ geNomad → chromosome / plasmid / phage classification
       ├─ MobMess → plasmid circularity + system inference
       ├─ Prodigal → ORF prediction (shared by AMR + summary)
       ├─ AMRFinderPlus → contig AMR annotation
       ├─ MobileElementFinder → contig MGE annotation
       ├─ MMseqs2 easy-taxonomy (UniRef50) → contig taxonomy
       ├─ minimap2 + samtools depth → contig cpg abundance
       └─ contig_summary.tsv (all per-contig annotations joined)

contigs + BAM
  └─ MetaBAT2 → MAG binning
       ├─ Prodigal → per-bin ORFs
       ├─ AMRFinderPlus → MAG AMR
       ├─ MobileElementFinder → MAG MGE
       ├─ GTDB-tk classify_wf → MAG taxonomy (optional)
       ├─ METABOLIC-G → MAG metabolic potential (optional)
       └─ CoverM genome → MAG abundance (trimmed mean)
```

## Requirements

- [Conda](https://docs.conda.io/) / [Mamba](https://mamba.readthedocs.io/) (recommended)
- [Snakemake](https://snakemake.readthedocs.io/) ≥ 7

All other software is installed automatically into per-rule conda environments on first run.

## Installation

```bash
git clone https://github.com/<org>/SWAM-meta.git
cd SWAM-meta
```

No additional installation step is required; `--use-conda` handles all dependencies.

## Configuration

Edit `config/config.yaml` before running:

| Key | Description |
|-----|-------------|
| `in_dir` | Directory containing paired-end FASTQ files (`*R1*.fastq*` / `*_1.fastq*`) |
| `out_dir` | Output directory (created if absent) |
| `scg_db` | Path to the 40 single-copy genes FASTA (`SCGs_40_All.fasta`) |
| `uniref50_db` | Prefix of a pre-built MMseqs2 UniRef50 taxonomy database (see below) |
| `gtdbtk_db` | Path to GTDB-tk reference data directory |
| `metabolic_dir` | Path to cloned METABOLIC repository |
| `genomad_splits` | Number of geNomad DB splits (increase to reduce memory; default `1`) |
| `skip_gtdbtk` | Set `True` to skip GTDB-tk (default `False`) |
| `skip_metabolic` | Set `True` to skip METABOLIC (default `False`) |

### Building the UniRef50 MMseqs2 database

```bash
mmseqs createdb uniref50.fasta PREFIX
mmseqs createtaxdb PREFIX tmp --ncbi-tax-dump TAXDIR
mmseqs createindex PREFIX tmp --search-type 2
```

Set `uniref50_db` to `PREFIX`.

### Downloading the GTDB-tk database

```bash
download-db.sh /path/to/gtdbtk_db
```

Set `gtdbtk_db` to the directory used above.

## Running the workflow

```bash
# Dry run (validate DAG without executing)
snakemake -n --use-conda

# Local run
snakemake --use-conda --cores <N>

# SLURM cluster
snakemake --use-conda --slurm --jobs <N>

# Force re-run a specific rule
snakemake --use-conda --forcerun <rule_name>
```

The first run downloads and indexes the AMRFinderPlus database (~300 MB) and the human reference genome (~900 MB) automatically via the `initiate_dbs` rule.

## Test mode

A self-contained test dataset (two mock samples, ≤ 16 GB RAM, ≤ 4 cores) is included under `test/`:

```bash
# Dry run
snakemake -n --use-conda --cores 4 --config test=True

# Full test run
snakemake --use-conda --cores 4 --config test=True
```

Test mode overrides all database paths to mini versions in `test/dbs/` and sets reduced resource limits. GTDB-tk and METABOLIC are skipped automatically.

## Outputs

| Path | Description |
|------|-------------|
| `{out_dir}/fastp_summary.csv` | Per-sample fastp QC metrics |
| `{out_dir}/short_reads_output.csv` | AMR gene abundances (cpg) per sample |
| `{out_dir}/data/megahit/{sample}.contigs.fa` | Assembled contigs (headers prefixed with sample name) |
| `{out_dir}/data/genomad/{sample}/` | geNomad classification output |
| `{out_dir}/data/mobmess/{sample}-mobmess_contigs.txt` | MobMess plasmid systems |
| `{out_dir}/contig_summary.tsv` | Joined per-contig table: abundance, molecule type, taxonomy, AMR, MGE |
| `{out_dir}/data/bins/{sample}/bins/` | MAG FASTA files |
| `{out_dir}/data/bins/{sample}/mag_amr.tsv` | AMR genes per MAG |
| `{out_dir}/data/bins/{sample}/mag_mge.tsv` | MGEs per MAG |
| `{out_dir}/data/bins/{sample}/gtdbtk.summary.tsv` | GTDB-tk taxonomy per MAG |
| `{out_dir}/data/bins/{sample}/mag_abundance.tsv` | Trimmed-mean relative abundance per MAG |

Clean reads (`.clean.fastq.gz`), Nonpareil `.npo` files, and assembled contigs are marked `protected()` to prevent accidental deletion.

## Abundance normalisation

AMR gene and contig abundances are expressed in **copies per genome (cpg)**:

```
n_genomes = Σ(alignment_length / gene_length) / 40    # across 40 SCGs
cpg       = depth / n_genomes
```

The same `n_genomes` estimate from the SCG DIAMOND alignment is reused for contig-level abundance.

## Repository structure

```
config/
  config.yaml          # user configuration
workflow/
  Snakefile            # entry point; defines rule all
  rules/
    common.smk         # sample discovery, test-mode overrides, resource helper
    short_reads.smk    # QC, host filtering, AMR + SCG alignment, coverage
    contigs.smk        # assembly, geNomad, MobMess, Prodigal, taxonomy, MGE, abundance
    mags.smk           # binning, per-MAG annotation, taxonomy, metabolism, abundance
  envs/                # per-rule conda environment YAML files
  scripts/
    short_reads_processing.R      # AMR normalisation + summary
    contig_abundance.py           # contig cpg abundance
    contig_summary.py             # join all contig annotations
    detect_circular_contigs.py    # Yu et al. 2024 circularity algorithm
test/
  data/                # mock FASTQ pairs (mock1, mock2)
  dbs/                 # mini test databases
  scripts/
    generate_mock_data.py         # reproducible mock data generation
```
