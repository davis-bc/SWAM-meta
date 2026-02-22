# SWAM-meta

Snakemake workflow for end-to-end metagenomic analysis of antibiotic resistance in environmental bacterial communities.

## Overview

SWAM-meta takes paired-end FASTQ files and produces a fully normalised, multi-evidence AMR abundance table alongside per-contig and per-MAG annotation outputs:

- **Short-read summaries** — AMR gene abundance (copies per genome, cpg) via KMA alignment against AMRFinderPlus, metagenomic coverage (Nonpareil), and QC metrics (fastp)
- **Assembled contigs** — de novo assembly, contig classification (chromosome / plasmid / phage), plasmid system inference (MobMess), taxonomy (MMseqs2/UniRef50), AMR (AMRFinderPlus), and MGE (MobileElementFinder) annotations joined into one TSV
- **Unified AMR table** — blends short-read sensitivity with contig-level context; cpg-normalised; evidence column tracks detection source
- **Metagenome-Assembled Genomes (MAGs)** — binning, AMR/MGE annotation, taxonomy (GTDB-tk, optional), metabolic potential (METABOLIC, optional), quality (CheckM2, optional), and relative abundance per MAG

---

## Pipeline stages

```
raw reads
  └─ fastp QC + minimap2 host removal
       ├─ KMA → AMR gene alignment (AMRFinderPlus nucleotide DB)
       ├─ DIAMOND → 40 single-copy genes (genome count estimation)
       ├─ KMA → anthropogenic markers (pBI143, crAss001)
       ├─ Nonpareil → metagenomic coverage
       └─ short_reads_output.csv + fastp_summary.csv + markers_cpg.csv

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

contig_summary.tsv + short_reads_output.csv
  └─ amr_unified.py → AMR_unified.csv + AMR_abundance_summary.csv

contigs + contig BAM
  └─ MetaBAT2 → MAG binning (checkpoint: bins discovered dynamically)
       ├─ Prodigal → per-bin ORFs
       ├─ AMRFinderPlus → MAG AMR
       ├─ MobileElementFinder → MAG MGE
       ├─ CheckM2 → MAG completeness / contamination (optional)
       ├─ GTDB-tk classify_wf → MAG taxonomy (optional)
       ├─ METABOLIC-G → MAG metabolic potential (optional)
       └─ CoverM genome → MAG abundance (trimmed mean, ≥10% coverage)
```

---

## Requirements

- [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/) (Mamba recommended for faster env solving)
- [Snakemake](https://snakemake.readthedocs.io/) ≥ 7

All other tools are installed automatically into isolated conda environments on first run (`--use-conda`).

> **Note:** The default Snakemake MILP scheduler requires the `cbc` solver. If it is not installed, add `--scheduler greedy` to all `snakemake` commands.

---

## Installation

```bash
git clone https://github.com/<org>/SWAM-meta.git
cd SWAM-meta
```

No additional build step is needed; `--use-conda` handles all tool installation.

---

## Quick start: test mode

Test mode runs the full pipeline on two synthetic mock samples using pre-built mini databases. It requires ≤ 16 GB RAM and ≤ 4 cores. No configuration changes are needed.

```bash
# 1. Dry run — validate the DAG without executing any jobs
snakemake -n --use-conda --cores 4 --scheduler greedy --config test=True

# 2. Full end-to-end test run
snakemake --use-conda --cores 4 --scheduler greedy --config test=True
```

Expected runtime: ~30–60 minutes on a laptop. Outputs are written to `test/output/`.

---

## Production setup

### 1. Prepare required databases

#### AMRFinderPlus + human genome
Downloaded automatically on first run by the `initiate_dbs` rule (~1.2 GB total). No manual action needed.

#### 40 single-copy gene database (SCGs)
Provide the path to `SCGs_40_All.fasta` in `config/config.yaml` under `scg_db`.

#### Anthropogenic markers (pBI143, crAss001)
Provide the path to a directory containing `pBI143.fasta` and `crAss001.fasta` under `markers_db`.

#### UniRef50 MMseqs2 taxonomy database (for contig taxonomy)

```bash
# Download UniRef50 FASTA from https://www.uniprot.org/uniref/
mmseqs createdb uniref50.fasta /path/to/uniref50_mmseqs
mmseqs createtaxdb /path/to/uniref50_mmseqs tmp --ncbi-tax-dump /path/to/taxdump
mmseqs createindex /path/to/uniref50_mmseqs tmp --search-type 2
```

Set `uniref50_db: /path/to/uniref50_mmseqs` in `config/config.yaml`.

#### GTDB-tk reference data *(optional)*

```bash
download-db.sh /path/to/gtdbtk_db
```

Set `gtdbtk_db: /path/to/gtdbtk_db`. To skip GTDB-tk entirely, set `skip_gtdbtk: True`.

#### METABOLIC *(optional)*

Clone the [METABOLIC repository](https://github.com/AnantharamanLab/METABOLIC) and set `metabolic_dir` to the cloned path. To skip, set `skip_metabolic: True`.

#### CheckM2 *(optional)*

```bash
checkm2 database --download --path /path/to/checkm2_db
```

Set `checkm2_db: /path/to/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd`. To skip, set `skip_checkm2: True`.

---

### 2. Configure `config/config.yaml`

```yaml
# Required
in_dir:      /path/to/fastq_files        # paired-end FASTQs (*R1*.fastq* or *_1.fastq*)
out_dir:     /path/to/output
scg_db:      /path/to/SCGs_40_All.fasta
markers_db:  /path/to/markers_directory  # must contain pBI143.fasta + crAss001.fasta
uniref50_db: /path/to/uniref50_mmseqs    # MMseqs2 DB prefix

# Optional large databases
gtdbtk_db:    /path/to/gtdbtk_db
metabolic_dir: /path/to/METABOLIC
checkm2_db:   /path/to/CheckM2_database/uniref100.KO.1.dmnd

# Performance tuning
genomad_splits: 1     # increase (e.g. 4–8) to reduce geNomad peak memory

# Skip optional stages
skip_gtdbtk:   False
skip_metabolic: False
skip_checkm2:  False
```

---

### 3. Run the workflow

```bash
# Dry run (validate DAG + config)
snakemake -n --use-conda --scheduler greedy

# Local run
snakemake --use-conda --cores <N> --scheduler greedy

# SLURM cluster
snakemake --use-conda --slurm --jobs <N> --scheduler greedy

# Resume after failure
snakemake --use-conda --cores <N> --scheduler greedy --rerun-incomplete

# Force re-run a specific rule
snakemake --use-conda --cores <N> --scheduler greedy --forcerun <rule_name>
```

The first run downloads and indexes the AMRFinderPlus database (~300 MB) and human reference genome (~900 MB) automatically.

---

## Outputs

| Path | Description |
|------|-------------|
| `{out_dir}/fastp_summary.csv` | Per-sample fastp QC metrics + Nonpareil coverage |
| `{out_dir}/short_reads_output.csv` | AMR gene abundances (cpg) from short-read KMA alignment |
| `{out_dir}/markers_cpg.csv` | pBI143 and crAss001 cpg per sample |
| `{out_dir}/contig_summary.tsv` | Per-contig: abundance, molecule type, taxonomy, AMR genes, MGEs |
| `{out_dir}/AMR_unified.csv` | **Unified AMR table** — short-reads + contigs merged (see below) |
| `{out_dir}/AMR_abundance_summary.csv` | Per-sample: total AMR cpg (unified) + pBI143 cpg + crAss001 cpg |
| `{out_dir}/data/megahit/{sample}.contigs.fa` | Assembled contigs (headers prefixed `{sample}-`) |
| `{out_dir}/data/genomad/{sample}/` | geNomad classification output |
| `{out_dir}/data/mobmess/{sample}-mobmess_contigs.txt` | MobMess plasmid systems |
| `{out_dir}/data/bins/{sample}/.binning.done` | MetaBAT2 binning sentinel |
| `{out_dir}/data/bins/{sample}/mag_amr.tsv` | AMR genes per MAG |
| `{out_dir}/data/bins/{sample}/mag_mge.tsv` | MGEs per MAG |
| `{out_dir}/data/bins/{sample}/gtdbtk.summary.tsv` | GTDB-tk taxonomy per MAG *(if enabled)* |
| `{out_dir}/data/bins/{sample}/mag_abundance.tsv` | Trimmed-mean relative abundance per MAG |

### AMR_unified.csv columns

| Column | Description |
|--------|-------------|
| `sample` | Sample identifier |
| `gene_symbol` | AMRFinderPlus allele (join key across both streams) |
| `gene_family` | Gene family grouping |
| `product_name` | Full gene product name |
| `class` / `subclass` | AMR drug class / subclass |
| `type` / `subtype` | Feature type (AMR, STRESS, etc.) |
| `cpg` | Preferred abundance: contig cpg if assembled, else short-reads cpg |
| `evidence` | Detection source: `both`, `short_reads_only`, or `contigs_only` |
| `short_reads_cpg` | Raw KMA-derived cpg (0 if not detected) |
| `contig_cpg` | Contig-derived cpg (0 if not assembled) |
| `molecule_type` | Comma-separated molecule types from geNomad (e.g. `chromosome,plasmid`) |
| `taxonomy` | Modal contig taxonomy from MMseqs2 |
| `n_contigs` | Number of distinct contigs carrying this gene (0 for short_reads_only) |

---

## Abundance normalisation

All AMR and contig abundances are expressed in **copies per genome (cpg)**:

```
n_genomes = Σ(alignment_length / gene_length) / 40    # across 40 SCGs via DIAMOND
cpg       = mean_depth / n_genomes
```

The same `n_genomes` estimate from the SCG DIAMOND alignment is shared between short-read AMR normalisation and contig-level abundance.

**AMR_unified.csv cpg priority:** when a gene is detected by both streams, the contig cpg is used (higher specificity — derived from assembled sequence); otherwise the short-reads cpg is used (detects genes present at too low depth to assemble).

---

## Repository structure

```
config/
  config.yaml                   # user configuration (edit before production run)
workflow/
  Snakefile                     # entry point; defines rule all
  rules/
    common.smk                  # sample discovery, test-mode overrides, resource helper
    short_reads.smk             # QC, host filtering, AMR + SCG + marker alignment
    contigs.smk                 # assembly, geNomad, MobMess, Prodigal, taxonomy, MGE, abundance
    mags.smk                    # binning, per-MAG annotation, taxonomy, metabolism, abundance
    summary.smk                 # cross-stage: AMR_unified + AMR_abundance_summary
  envs/                         # per-rule conda environment YAML files
  scripts/
    short_reads_processing.R    # AMR + marker normalisation; outputs short_reads_output.csv + markers_cpg.csv
    contig_abundance.py         # contig cpg abundance from BAM depth
    contig_summary.py           # joins geNomad, MMseqs2, AMR, MGE, abundance per sample
    amr_unified.py              # merges short-read + contig AMR; generates AMR_unified.csv + AMR_abundance_summary.csv
    detect_circular_contigs.py  # Yu et al. 2024 circularity algorithm (used by MobMess rule)
test/
  data/                         # mock FASTQ pairs (mock1, mock2) — biologically realistic
  dbs/                          # mini test databases (AMRFinderPlus subset, UniRef50 subset, etc.)
  scripts/
    generate_mock_data.py       # reproducible mock data generation script
docs/
  session-log.md                # Copilot session log
```

