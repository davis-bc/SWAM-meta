# SWAM-meta Copilot Instructions

## Overview

SWAM-meta is a Snakemake workflow for end-to-end metagenomic analysis of antibiotic resistance in environmental bacterial communities. Full pipeline: read QC and host filtering → AMR gene alignment (short-read) → metagenomic coverage estimation → de novo assembly → contig classification (chromosome/plasmid/phage) → taxonomic assignment → MGE annotation → contig abundance → contig summary → MAG binning → MAG characterisation (AMR, MGE, taxonomy, metabolism, abundance).

## Running the Workflow

```bash
# Dry run (validate DAG without executing)
snakemake -n --use-conda

# Full run (local)
snakemake --use-conda --cores <N>

# Run on SLURM cluster
snakemake --use-conda --slurm --jobs <N>

# Run a specific rule for a specific sample
snakemake --use-conda output/<sample>.contigs.fa

# Force re-run a specific rule
snakemake --use-conda --forcerun <rule_name>
```

No test suite exists yet; dry runs (`-n`) serve as validation. A `--config test=True` mode is planned (see **Testing** section).

```bash
# Run in test mode with mock samples (≤16 GB RAM)
snakemake --use-conda --cores 4 --config test=True
```

## Configuration

`config/config.yaml` is the only user-facing configuration file. Required keys:
- `in_dir`: path to directory containing paired-end FASTQ files
- `out_dir`: path for all outputs
- `scg_db`: path to the 40 single-copy genes FASTA (`SCGs_40_All.fasta`)

Keys to be added for planned stages:
- `uniref50_db`: path to pre-built MMseqs2 UniRef50 taxonomy database (`mmseqs createdb` + `createtaxdb`)
- `gtdbtk_db`: path to GTDB-tk reference data directory (set `GTDBTK_DATA_PATH`)
- `metabolic_db`: path to METABOLIC database directory

## Architecture

### Pipeline stages (in dependency order)

1. **`initiate_dbs`** (`short_reads.smk`) — one-time setup: downloads AMRFinderPlus DB + human reference genome, indexes them with KMA and DIAMOND. Outputs sentinel `.done.txt` files.
2. **`short_reads`** (`short_reads.smk`) — per-sample: fastp QC → minimap2 human-read removal → KMA alignment to AMRFinderPlus → DIAMOND alignment to SCGs → nonpareil coverage estimation.
3. **`short_reads_summary`** (`short_reads.smk`) — aggregates all samples into two CSVs: `fastp_summary.csv` and `short_reads_output.csv` (AMR genes normalized to copies-per-genome). Runs an R script via Snakemake's `script:` directive.
4. **`init_genomad`** + **`contigs`** (`contigs.smk`) — MEGAHIT assembly (per sample, `--presets meta-large`, min contig 1000 bp); contig headers are renamed `{sample}-<original>` in-place.
5. **`genomad`** (`contigs.smk`) — classifies contigs into chromosome/plasmid/phage using geNomad (`--relaxed --cleanup`).
6. **`mobmess`** (`contigs.smk`) — maps reads back to plasmid contigs, calls `detect_circular_contigs.py` to flag circular sequences, then runs MobMess to infer plasmid systems.

**Planned stages (not yet implemented):**

7. **`mmseqs_taxonomy`** (`contigs.smk`) — per-sample: MMseqs2 `easy-taxonomy` against a pre-built UniRef50 taxonomy database to assign taxonomy to each contig. Output: `{sample}_lca.tsv` (contig → lineage). Use `--lca-ranks superkingdom,phylum,class,order,family,genus,species`. The UniRef50 MMseqs2 DB must be pre-built externally and pointed to via `config["uniref50_db"]`.
8. **`mge_annotation`** (`contigs.smk`) — per-sample: MobileElementFinder (`mobilefinder`) to annotate MGEs (insertion sequences, transposons, etc.) on assembled contigs. Output: `{sample}_mge.tsv` with contig, start, stop, strand, MGE type.
9. **`contig_abundance`** (`contigs.smk`) — per-sample: map clean reads back to the assembled contigs with minimap2 (`-ax sr`), sort/index with samtools, then compute per-contig depth with `samtools depth`. Normalise to cpg using the same `n_genomes` estimate from the SCG DIAMOND alignment already produced by the `short_reads` rule. Output: `{sample}_contig_abundance.tsv`.
10. **`contig_summary`** (`contigs.smk`) — aggregation rule (all samples): joins per-sample outputs from steps 5–9 into a single TSV with columns: `sample | contig_id | abundance_cpg | molecule_type | taxonomy | feature_type | gene | start | stop | strand`. `feature_type` distinguishes AMR (from `prodigal` + AMRFinderPlus protein search on contigs) from MGE (from MobileElementFinder). Implemented as a Python/R script called via Snakemake `script:` directive.
11. **`bin`** (`mags.smk`, new file) — per-sample: index the contig-to-read BAM from step 9, run MetaBAT2 (`jgi_summarize_bam_contig_depths` → `metabat2`) to produce draft bins. Min contig size for binning: 2500 bp. Output: `output/data/bins/{sample}/`.
12. **`mag_amr`** (`mags.smk`) — per-bin: AMRFinderPlus protein mode on Prodigal-predicted ORFs. Requires Prodigal run per bin first.
13. **`mag_mge`** (`mags.smk`) — per-bin: MobileElementFinder on bin contigs.
14. **`mag_taxonomy`** (`mags.smk`) — per-sample batch: GTDB-tk `classify_wf` on all bins for a sample. Requires `GTDBTK_DATA_PATH` set from `config["gtdbtk_db"]`. Output: `gtdbtk.summary.tsv`.
15. **`mag_metabolism`** (`mags.smk`) — per-sample batch: METABOLIC-G.pl on all bin FASTAs. Output: `METABOLIC_result/`.
16. **`mag_abundance`** (`mags.smk`) — per-sample: CoverM `genome` mode with `--min-covered-fraction 0.1`, trimmed-mean method, against the contig-to-read BAM from step 9. Output: `{sample}_mag_abundance.tsv`.

The `rule all` target in `Snakefile` will need updating to include the new terminal outputs (contig summary TSV and per-sample MAG abundance TSVs).

### Sample discovery (`rules/common.smk`)

Samples are discovered automatically at DAG construction time by globbing `in_dir` for `*R1*.fastq*` / `*_1.fastq*` patterns and pairing with their R2 counterparts. The `sample` wildcard is constrained to `[^/]+`.

### Conda environments

Each rule specifies a conda env in `workflow/envs/`:
- `shortreads.yaml` — fastp, kma, nonpareil, diamond, minimap2, samtools, pigz
- `contigs.yaml` — megahit, minimap2, samtools, prodigal, mmseqs2, blast, pysam, MobMess + PlasX (pip from GitHub)
- `genomad.yaml` — genomad
- `Renv.yaml` — R, tidyverse, Nonpareil (for summary script)

Environments to be added for planned stages:
- `mge.yaml` — mobilefinder (MobileElementFinder)
- `mags.yaml` — metabat2, coverm, prodigal, ncbi-amrfinderplus, mobilefinder, gtdbtk
- `metabolic.yaml` — METABOLIC (separate env due to Perl dependency complexity)

### Key scripts

- `workflow/scripts/short_reads_processing.R` — called via `snakemake@input`/`snakemake@output` (Snakemake R script interface, not standalone). Normalizes AMR gene depth to copies-per-genome (cpg) using 40 SCGs as a proxy for genome count.
- `workflow/scripts/detect_circular_contigs.py` — standalone Python script (pysam). Implements the Yu et al. 2024 circularity algorithm: reads mapped in reverse-forward orientation spanning ≥ `contig_length − 3×median_insert` indicate circularity.

Scripts to be added for planned stages:
- `workflow/scripts/contig_summary.py` (or `.R`) — joins geNomad molecule type, MMseqs2 taxonomy, contig cpg abundance, AMR annotations, and MGE annotations into a single long-format TSV. One row per genomic feature (AMR gene or MGE); contig-level fields (abundance, molecule type, taxonomy) are repeated per feature row.

## Testing

The pipeline needs built-in test data to validate end-to-end on a laptop/workstation (≤16 GB RAM, ≤4 cores).

### Design

- Test FASTQ pairs live in `test/data/` (two mock samples: `mock1_R1.fastq.gz` + `mock1_R2.fastq.gz`, `mock2_R1.fastq.gz` + `mock2_R2.fastq.gz`).
- Each mock sample: ~50,000 read pairs, 150 bp, subsampled from a real metagenome or fully synthetic.
- Test databases live in `test/dbs/` (tiny versions: e.g. 100-sequence AMRFinderPlus subset, minimal UniRef50 MMseqs2 db, SCG subset).
- Activated via `--config test=True` in `common.smk`: when `config.get("test")` is truthy, override `input_dir`, `output_dir`, and all db paths with the `test/` equivalents.
- Resource overrides for test mode: set `mem_mb=4000`, `threads=2`, `time="0-01:00:00"` via a conditional in each rule's `resources:` block, or via a Snakemake `--resources` flag.
- `localrules` list must include any db-init rules that are skipped/mocked in test mode.
- Test run command: `snakemake --use-conda --cores 4 --config test=True -n` (dry) or without `-n` for a real end-to-end test.
- Mock FASTQ generation script: `test/scripts/generate_mock_reads.py` — uses `seqtk sample` or pure Python to subsample/synthesize reads; document the exact commands used so the test data is reproducible.

### Efficiency amendments to existing rules

When implementing the new stages, also apply these improvements to current rules:
- **`short_reads`**: the concatenated single-end temp file for DIAMOND (`TMP_SE`) can be replaced with `--query <(cat r1 r2)` process substitution to avoid writing ~doubled data to disk.
- **`contigs`**: the `--continue` flag on MEGAHIT is already set; ensure the BAM produced for `mobmess` is reused as the contig-abundance BAM (step 9) to avoid a redundant minimap2 mapping. The mobmess rule currently maps reads to *plasmid* contigs only — a separate mapping to *all* contigs is needed for abundance; keep them as distinct outputs.
- **`genomad`**: pass `--splits 8` on low-memory systems; make this configurable via `config.get("genomad_splits", 1)`.
- **Prodigal**: run once per sample on all contigs (not per-rule) and share the output `.faa` / `.gff` between AMR annotation and the contig summary. Add a dedicated `prodigal` rule in `contigs.smk` whose output is consumed by both `contig_amr` and `contig_summary`.
- **MAG binning reuse**: the contig-depth file from `jgi_summarize_bam_contig_depths` (MetaBAT2 step) is the same depth information needed for CoverM — write it once and reference it in both rules.

## Key Conventions

- **Protected outputs**: clean reads (`.clean.fastq.gz`), nonpareil `.npo` files, and assembled contigs are marked `protected()` to prevent accidental deletion.
- **Resource declarations**: all compute-intensive rules declare `mem_mb`, `threads`, and `time` under `resources:` for cluster schedulers.
- **Temporary files**: rules use `$TMPDIR` set to the output directory to keep temp files on the same filesystem, enabling atomic moves. Temp files are explicitly `rm -f`'d in shell blocks.
- **`localrules`**: `initiate_dbs`, `short_reads_summary`, and `init_genomad` are declared local (run on the head node, not submitted to the cluster).
- **`.gitignore`**: `input/`, `output/`, `slurm/`, `dbs/`, `.snakemake/`, and `run_snakemake.sh` are excluded — these are environment-specific and should not be committed.
- **MobMess/PlasX** are installed via pip from GitHub inside the conda env — internet access is required on first `--use-conda` setup.
- **cpg normalisation**: copies-per-genome is computed as `alignment_depth / n_genomes` where `n_genomes = sum(alignment_length / gene_length) / 40` across the 40 SCGs. Apply this same formula to contig-level depth for the contig abundance step.
- **New rule file convention**: MAG-related rules go in a new `workflow/rules/mags.smk`, included in `Snakefile`. Add its terminal outputs to `rule all`.
- **`prodigal` shared output**: add a `prodigal` rule in `contigs.smk` that outputs `{sample}.faa` and `{sample}.gff`; both the contig AMR step and the contig summary script consume these files. Do not run Prodigal redundantly per tool.
- **Test mode pattern**: in `common.smk`, check `config.get("test", False)` and override `input_dir`, `output_dir`, and all database path variables with `test/` paths. Resource values in each rule should be wrapped: `resources: mem_mb = 4000 if config.get("test") else 150000` etc.
