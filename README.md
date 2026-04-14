# SWAM-meta

Snakemake workflow for end-to-end metagenomic analysis of antibiotic resistance in environmental bacterial communities.

## What it produces

SWAM-meta starts from paired-end FASTQs and generates:

- **`fastp_summary.csv`** - per-sample read QC and Nonpareil coverage
- **`assembly_qa.tsv`** - per-sample assembly size, N50, and read-mapping stats
- **`contig_summary.tsv`** - per-contig abundance, taxonomy, AMR, MGE, and molecule type
- **`AMR_abundance_summary.csv`** - per-sample AMR abundance plus additive and multiplicative risk scores
- **`mag_summary.tsv`** - per-MAG abundance and annotations

It also keeps lower-level outputs under `out_dir/data/`, including assemblies, geNomad results, MobMess output, MAG files, and the short-read tables `data/QAQC/short_reads_output.csv` and `data/QAQC/markers_cpg.csv`.

---

## Pipeline at a glance

1. **Short reads** - fastp QC, host filtering, KMA AMR alignment, SCG alignment, anthropogenic marker alignment, Nonpareil
2. **Contigs** - MEGAHIT assembly, geNomad classification, MobMess, Prodigal, AMRFinderPlus, MobileElementFinder, MMseqs2 taxonomy, contig abundance
3. **MAGs** - MetaBAT2 binning, per-bin AMR/MGE annotation, abundance, and optional GTDB-tk, CheckM2, and METABOLIC
4. **Summary** - sample-level AMR abundance and risk scoring from AMR, mobility, host context, and marker-derived exposure

[![SWAM-meta rulegraph](docs/rulegraph.png)](docs/rulegraph.png)

---

## Installation

```bash
git clone https://github.com/<org>/SWAM-meta.git
cd SWAM-meta
conda install -c conda-forge -c bioconda "snakemake>=8" snakemake-executor-plugin-slurm
```

All workflow tools are installed automatically into per-rule conda environments on first run with `--use-conda`.

> **Always use `--scheduler greedy`.** The default MILP scheduler requires the `cbc` solver, which is not installed.

---

## Quick start

### Test mode

This runs the full workflow on two bundled mock samples and mini databases.

```bash
snakemake -n --use-conda --cores 4 --scheduler greedy --config test=True
snakemake    --use-conda --cores 4 --scheduler greedy --config test=True
```

Expected runtime is about 30-60 minutes on a laptop. Outputs go to `test/output/`.

### Production mode

Edit `config/config.yaml`:

```yaml
in_dir:  /path/to/fastq_files
out_dir: /path/to/output
```

Then run:

```bash
snakemake -n --use-conda --scheduler greedy
snakemake    --use-conda --cores <N> --scheduler greedy
```

To resume after a failed run:

```bash
snakemake --use-conda --cores <N> --scheduler greedy --rerun-incomplete
```

---

## Databases

Most databases are downloaded and prepared automatically in `dbs/` on first use, including:

- AMRFinderPlus databases
- human reference for host filtering
- anthropogenic markers (`pBI143`, `crAss001`)
- geNomad database
- UniRef50 MMseqs2 taxonomy database
- CheckM2 database
- METABOLIC

**GTDB-tk is the only manual setup step.**

```bash
conda run -n gtdbtk download-db.sh /path/to/gtdbtk_db
```

Then set:

```yaml
gtdbtk_db: /path/to/gtdbtk_db
```

Optional stages can be disabled in `config/config.yaml`:

```yaml
skip_gtdbtk: False
skip_metabolic: False
skip_checkm2: False
genomad_splits: 1
```

Increase `genomad_splits` if geNomad needs less peak memory.

---

## Running on SLURM

Two profiles are included:

| Profile | Best for |
|---|---|
| `config/slurm/small-batch` | Fewer than 50 samples |
| `config/slurm/large-batch` | 50 or more samples |

Set your account and partition in each profile, then run one of:

```bash
snakemake --profile config/slurm/small-batch
snakemake --profile config/slurm/large-batch
```

The profiles already enable `use-conda` and the required greedy scheduler.

---

## Running only part of the workflow

Use these config flags:

```yaml
run_short_reads: True
run_contigs: True
run_mags: True
```

| `run_short_reads` | `run_contigs` | `run_mags` | Result | Extra input needed |
|---|---|---|---|---|
| `True` | `True` | `True` | Full pipeline | None |
| `True` | `True` | `False` | Short reads + contigs | None |
| `True` | `False` | `False` | Short reads only | None |
| `False` | `True` | `True` | Contigs + MAGs | `clean_reads_dir` |
| `False` | `True` | `False` | Contigs only | `clean_reads_dir` |
| `False` | `False` | `True` | MAGs only | `clean_reads_dir`, `contigs_dir` |

When skipping upstream stages:

- `clean_reads_dir` must contain host-filtered paired FASTQs
- `contigs_dir` must contain `{sample}.contigs.fa`

If `run_short_reads: False`, contig abundance falls back to mean depth rather than fully normalized cpg because the SCG-based genome estimate is unavailable.

---

## Main outputs

| Path | Description |
|---|---|
| `{out_dir}/fastp_summary.csv` | Per-sample fastp metrics and Nonpareil coverage |
| `{out_dir}/assembly_qa.tsv` | Per-sample assembly QC metrics |
| `{out_dir}/contig_summary.tsv` | Per-contig annotations and abundance |
| `{out_dir}/AMR_abundance_summary.csv` | Per-sample AMR abundance and risk scores |
| `{out_dir}/mag_summary.tsv` | Per-MAG summary table |
| `{out_dir}/data/QAQC/short_reads_output.csv` | Short-read AMR gene abundances |
| `{out_dir}/data/QAQC/markers_cpg.csv` | `pBI143` and `crAss001` copies per genome |

`AMR_abundance_summary.csv` contains:

| Column | Meaning |
|---|---|
| `sample` | Sample ID |
| `AMR_total_cpg` | Total AMR abundance across detected genes |
| `pBI143_cpg`, `crAss001_cpg` | Anthropogenic marker abundance |
| `E_exposure` | Exposure score from marker abundance |
| `R_mean` | Mean resistance-hazard score |
| `M_mean` | Mean mobility score |
| `H_mean` | Mean host/pathogenicity score |
| `amr_risk_additive_raw`, `amr_risk_multiplicative_raw` | Raw study-level risk scores |
| `amr_risk_additive`, `amr_risk_multiplicative` | Min-max normalized 0-100 scores |

---

## Abundance and risk scoring

All AMR and contig abundances are reported as **copies per genome (cpg)**:

```text
n_genomes = Σ(alignment_length / gene_length) / 40
cpg       = mean_depth / n_genomes
```

The 40-gene denominator comes from the bundled single-copy gene reference in `workflow/resources/SCGs_40_All.fasta`.

Risk scoring combines four components:

- **R** - resistance hazard from AMRFinderPlus class/subclass
- **M** - mobility from contig type and nearby MGE evidence
- **H** - host/pathogenicity from MMseqs2 taxonomy
- **E** - exposure from `pBI143` and `crAss001`

Both additive and multiplicative sample-level scores are reported.

---

## Repository layout

```text
config/                 user configuration and SLURM profiles
docs/                   rule graph and session log
test/                   mock data, mini databases, reproducible test scripts
workflow/
  Snakefile             workflow entry point
  envs/                 per-rule conda environments
  resources/            bundled reference files
  rules/                short-read, contig, MAG, and summary rules
  scripts/              Python/R scripts used by rules
dbs/                    auto-managed databases (gitignored)
```
