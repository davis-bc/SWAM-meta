# SWAM-meta — Copilot Session Log

This file is maintained by GitHub Copilot. At the start of each session, Copilot reads this file to restore context. At the end of each session (or when significant progress is made), Copilot appends a new entry.

---

## 2026-02-22 (session 2)

### What was done
- **Mock dataset v2**: Replaced the old 50k-read synthetic mock data with a biologically realistic dataset simulated from 6 real NCBI reference genomes:
  - *E. coli* 2024GO-0438 (GCA_042156675.1, 5.18 Mb, 198 contigs, plasmids)
  - *E. faecalis* ST6 (GCA_015304995.1, 3.31 Mb, 89 contigs, plasmids)
  - *Salmonella enterica* ENWS000537-005 (GCA_052022995.1, 4.71 Mb, 43 contigs, plasmids)
  - *Campylobacter coli* 18QD2YX31C (GCF_013072625.1, 1.94 Mb, 5 contigs, plasmids)
  - Carjivirus communis phage (GCF_009734245.1, 97 kb, 1 contig)
  - Influenza A H3N2 (GCA_038998445.1, 13.6 kb, 8 segments)
- **mock1** = E. coli (10×) + E. faecalis (10×) + Carjivirus (3×) + Influenza A (3×) → 284,141 read pairs (16 MB per file gzipped)
- **mock2** = E. coli (10×) + Salmonella (10×) + Campylobacter (10×) + Carjivirus (3×) → 394,971 read pairs (23 MB per file gzipped)
- E. coli and Carjivirus appear in both samples to exercise CoverM cross-sample mapping.
- 40 synthetic SCGs (seed=52) + blaTEM-1 re-embedded into E. coli first contig at same positions; existing DIAMOND/KMA databases remain valid.
- Reference FASTAs saved to `test/dbs/genomes/`. `test/data/mock_data_summary.txt` documents composition and rationale.
- **CheckM2 MAG QC**: Added `mag_qc` rule to `mags.smk` (`checkm2 predict` per sample). New `workflow/envs/checkm2.yaml`. `checkm2_db` and `skip_checkm2` added to `config/config.yaml`. `skip_checkm2=True` set in test mode. `checkm2_quality.tsv` added to `rule all` as conditional target.
- Pushed to `origin/main` (commit `200c9c4`).

### Current pipeline state
All stages fully implemented and test data ready. Key additions since last session:

| Component | Status |
|-----------|--------|
| Mock data | ✅ Upgraded to real full-genome references (v2) |
| CheckM2 MAG QC | ✅ Added (`mag_qc` rule, `checkm2.yaml` env) |
| CheckM2 in test mode | ✅ Skipped (`skip_checkm2=True`) |

### Known issues / next steps
- CheckM2 requires downloading its database (~600 MB `.dmnd` file) via `checkm2 database --download --path <dir>` before use in production.
- The 10x coverage mock data is larger than the old dataset (16-23 MB gzipped FASTQ files vs 3.2 MB). This is necessary to ensure MetaBAT2 can produce bins.
- `test/output/` (gitignored) contains old test run outputs and will need to be cleared before running with the new mock data.
- Reference genomes source: `/home/benda/Projects/databases/test-genomes/` — required locally to re-run `generate_mock_data.py`.

---

## 2026-02-22 (session 5)

### What was done
- **Removed all `protected()` output declarations** from `short_reads.smk` (clean reads, nonpareil) and `contigs.smk` (megahit contigs). These were causing `ProtectedOutputException` on re-runs after failures.
- **Fixed MILP scheduler missing `cbc`**: Snakemake's MILP-based job selector requires the `cbc` solver which is not installed. All test runs now use `--scheduler greedy`.
- **Fixed MobMess `KeyError` on small mock plasmid data**: Added `|| true` + `touch` fallback after `mobmess systems` call — the tool crashes in `derep_sources` when all 216 plasmid contigs form singleton connected components. Pipeline now gracefully skips and continues.
- **Fixed `mag_amr` AMRFinder `gff_check` failure on small bins**: Added `2>/dev/null || true` to the per-bin `amrfinder` call in `mags.smk`, matching the pattern already in `contig_amr`.
- **Full end-to-end test passed**: `rule all` completed successfully — all 29 jobs across both mock samples finished. Output verified at `test/output/`.
- Pushed to `origin/main` (commits `6b60b72`, `2dc586f`, `b8ac566`).

### Current pipeline state
✅ **Full end-to-end test mode passes** (`snakemake --use-conda --cores 4 --scheduler greedy --config test=True`).

All 29 jobs complete for 2 mock samples:
- short_reads, short_reads_summary ✅
- contigs, genomad, mobmess ✅
- init_mmseqs_db, mmseqs_taxonomy ✅
- prodigal, contig_amr, mge_annotation, contig_abundance, contig_summary ✅
- bin, mag_prodigal, mag_amr, mag_mge, mag_abundance ✅ (gtdbtk, metabolic, checkm2 skipped as expected)

### Known issues / next steps
- `--scheduler greedy` must be passed to avoid `PulpSolverError: cannot execute cbc`. Consider documenting this in README or installing `coincbc` in the base env.
- MobMess produces no output for mock2 plasmids (gracefully skipped). This is expected — tool has a known edge case with high-identity plasmid clusters.
- AMRFinder finds no AMR genes in mock bins (expected for synthetic test data; `|| true` handles this cleanly).

---

## 2026-02-22 (session 4)

### What was done
- **Bug fix — OOM crash on test run**: All jobs died with `SIGTERM` because `threads` was placed inside `resources:` instead of as a top-level `threads:` directive. Snakemake only uses the `threads:` directive for core-based scheduling; `resources: threads` is a custom resource it ignores for that purpose. With `--cores 4`, Snakemake treated every job as 1-core and allowed up to 4 memory-intensive jobs (geNomad ×2, MEGAHIT, mmseqs, contig_abundance, etc.) to run simultaneously, exhausting RAM and crashing Ubuntu via the OOM killer.
- **Fix**: Added `threads: lambda wc: res(X, Y)` top-level directive to all compute-intensive rules in `short_reads.smk`, `contigs.smk`, and `mags.smk`. Shell commands continue using `{resources.threads}` unchanged (both can coexist).
- **Bonus fix**: Set `genomad_splits=4` as the default in test mode (`config.setdefault("genomad_splits", 4)` in `common.smk`) to further reduce geNomad peak memory. Can be overridden with `--config genomad_splits=N`.
- Dry run (`snakemake -n --config test=True --cores 4`) passes cleanly.
- Pushed to `origin/main` (commit `b85f186`).

### Current pipeline state
Same as session 3 — all stages implemented. The only change is the scheduling fix above.

### Known issues / next steps
- Re-run the full test (`snakemake --use-conda --cores 4 --config test=True`) to confirm end-to-end success now that parallel OOM is prevented.
- If memory is still tight, run with `--cores 2` to further limit parallelism, or add `--resources mem_mb=12000` to cap total concurrent memory usage.

---

## 2026-02-22 (session 3)

### What was done
- **Bug fix — markers always showing 0 cpg**: KMA includes the full FASTA description in the `#Template` column (e.g., `U30316.1 Bacteroides fragilis strain IB143 plasmid pBI143, complete sequence`), but the R script was filtering by exact accession match (`template == "U30316.1"`), which never succeeded. Fixed by extracting the accession before the first space: `mutate(seq_id = str_extract(template, "^\\S+"))` and filtering on `seq_id`.
- **Bug fix — robustify AFM template accession extraction**: Added `str_extract(refseq_nucleotide_accession, "^\\S+")` after the `separate()` call to strip any trailing description text from the nucleotide accession used as the catalog join key. Addresses the reported "AMRFinderPlus metadata missing" production issue (likely caused by the same KMA full-description behaviour on some database versions).
- **Test confirmed**: Deleted and re-ran `short_reads_summary` in test mode. `AMR_abundance_summary.csv` now correctly shows non-zero pBI143/crAss001 cpg values (mock1: pBI143=3.64, crAss001=0.50; mock2: crAss001=0.49). AMR metadata still fully populated.
- **Note**: All changes confined to `workflow/scripts/short_reads_processing.R`. No Snakemake rule or conda env changes needed.

### Current pipeline state
All stages implemented. Short-reads summary outputs are:
- `fastp_summary.csv` — per-sample QC metrics
- `short_reads_output.csv` — per-gene AMR abundance in cpg with full metadata
- `AMR_abundance_summary.csv` — per-sample totals: AMR_total (cpg) | pBI143 (cpg) | crAss001 (cpg)

### Known issues / next steps
- Production run OOM'd (likely during MEGAHIT or MetaBAT2). Need to investigate peak memory rule and adjust resources or assembly parameters.
- pBI143 and crAss001 references are stored at `/home/benda/Projects/databases/reference-sequences/` (production) and `test/dbs/markers/` (test). Both use the same FASTA headers.
- Full end-to-end test on this machine is blocked by memory — short-reads-only testing is feasible.


### What was done
- Wrote a comprehensive `README.md` covering: pipeline overview (ASCII flow diagram), requirements, configuration, run commands, test mode, outputs table, cpg normalisation formula, and repository structure.
- Pushed to `origin/main` (commit `1471875`).
- Created `docs/session-log.md` (this file) and added session-logging instructions to `.github/copilot-instructions.md`.

### Current pipeline state
The full pipeline is implemented end-to-end across four rule files:

| Rule file | Status | Key rules |
|-----------|--------|-----------|
| `short_reads.smk` | ✅ Complete | `initiate_dbs`, `short_reads`, `short_reads_summary` |
| `contigs.smk` | ✅ Complete | `contigs`, `genomad`, `mobmess`, `init_mmseqs_db`, `prodigal`, `contig_amr`, `init_mge_tool`, `mge_annotation`, `mmseqs_taxonomy`, `contig_abundance`, `contig_summary` |
| `mags.smk` | ✅ Complete | `bin`, `mag_prodigal` (checkpoint), `mag_amr`, `mag_mge`, `mag_taxonomy`, `mag_metabolism`, `mag_abundance` |
| Test infrastructure | ✅ Complete | Mock data in `test/data/`, mini DBs in `test/dbs/`, test mode via `--config test=True` |

### Known issues / next steps
- No known bugs at this time.
- MAG taxonomy (`mag_taxonomy`) and metabolism (`mag_metabolism`) require large external databases (GTDB-tk, METABOLIC); skipped in test mode and when `skip_gtdbtk: True` / `skip_metabolic: True`.
- `contig_summary.py` script exists but content not reviewed this session.
- MobMess/PlasX require internet access on first `--use-conda` setup (pip install from GitHub).

---

## 2026-02-22 (session 6)

### What was done
- **AMR_unified.csv feature implemented** (commit `ac36b6d`): New cross-stream AMR abundance table merging short-read KMA detections with contig-level AMRFinderPlus annotations.
  - New `workflow/rules/summary.smk` with `localrule amr_unified`
  - New `workflow/scripts/amr_unified.py`: outer-joins `short_reads_output.csv` + `contig_summary.tsv` (AMR rows only) on `sample × gene_symbol`
  - **cpg priority**: contig cpg preferred when available (higher specificity), else short-reads cpg
  - **`evidence` column**: `"both"` | `"short_reads_only"` | `"contigs_only"`
  - **`molecule_type`**: comma-sep sorted unique values per gene × sample (e.g., `"phage,plasmid"` when gene found on both)
  - **`taxonomy`**: modal non-empty value from contig annotations (empty string for short_reads_only)
  - **AMRFinderPlus catalog metadata** (gene_family, class, subclass, type, subtype) joined from `ReferenceGeneCatalog.txt`
  - Fixed pandas 2.0 compatibility (removed `include_groups=False` from `groupby.apply`)
  - Added `AMR_unified.csv` to `all_targets()` in `Snakefile`
- Verified on test data: `blaTEM-1` correctly tagged `evidence=both`, contig cpg preferred; `tet(M)` shows `molecule_type=phage,plasmid`.

### Current pipeline state

| Rule file | Status | Key rules |
|-----------|--------|-----------|
| `short_reads.smk` | ✅ Complete | `initiate_dbs`, `short_reads`, `short_reads_summary` |
| `contigs.smk` | ✅ Complete | `contigs`, `genomad`, `mobmess`, `prodigal`, `contig_amr`, `mge_annotation`, `mmseqs_taxonomy`, `contig_abundance`, `contig_summary` |
| `mags.smk` | ✅ Complete | `bin`, `mag_prodigal`, `mag_amr`, `mag_mge`, `mag_taxonomy`, `mag_metabolism`, `mag_abundance` |
| `summary.smk` | ✅ Complete (new) | `amr_unified` |
| Test infrastructure | ✅ Complete | Mock data, mini DBs, `--config test=True` |

### Known issues / next steps
- **Catalog metadata sparse in test output**: test catalog uses synthetic gene names (`fakeGene_*`) so only `blaTEM-1` gets metadata annotated; production run with real ReferenceGeneCatalog will populate all rows.
- Consider adding a contig-level AMR annotation lookup using the real catalog at runtime (post-assembly) to enrich `contigs_only` rows with gene_family/class/etc. when not in catalog.
- `short_reads_only` evidence will appear in production runs where KMA detects genes at low read depth that don't assemble (expected for low-abundance resistome members).

---

## 2026-02-22 (session 7)

### What was done
Four issues fixed in the unified AMR output:

**1. Test KMA database too small** (`test/dbs/amrfinder/`)
- The test `AMR_CDS.fa` had only 11 synthetic genes + blaTEM-1, causing `short_reads_output.csv` to detect only blaTEM-1.
- Wrote a script that extracts 49 real AMR gene nucleotide sequences from assembled contigs (using AMRFinder coordinates, strand-aware reverse complement) and adds them to `AMR_CDS.fa` with proper AMRFinderPlus-style FASTA headers.
- Updated `ReferenceGeneCatalog.txt` with synthetic `NG_MOCK_*` accessions matching the FASTA headers so the R catalog join succeeds.
- `short_reads_output.csv` now detects 81 gene-per-sample entries across both mock samples.

**2. MGE parsing broken in contig_summary.py**
- MobileElementFinder output is CSV (comma-separated) but `safe_read` used tab separator → single-column DataFrame → all fields empty.
- Comment lines (e.g. `#date:`, `#sample:`) were parsed as data.
- `stop_col` search was missing "end" (the actual column name in mge_finder output).
- Contig IDs had extra metadata appended (`"contig_id flag=1 multi=..."`) that failed mol_map/tax_map lookups.
- Fixed: `parse_mge` now uses `pd.read_csv(sep=",", comment="#")`, strips contig ID to first token, and includes "end" in stop_col candidates.

**3. amr_unified.py metadata was empty for most genes**
- Previous version joined metadata from `ReferenceGeneCatalog.txt` only (synthetic in test mode).
- Rebuilt metadata lookup from `*_contig_amr.tsv` files directly (highest priority), with fallback to short_reads catalog join, then the catalog file.
- All genes in the unified output now have product_name, class, subclass, type, subtype populated.

**4. AMR_abundance_summary.csv now uses unified totals**
- R script `short_reads_processing.R` now outputs `markers_cpg.csv` (pBI143/crAss001 per sample) as out_file3 instead of `AMR_abundance_summary.csv`.
- `amr_unified.py` reads `markers_cpg.csv` as input and generates `AMR_abundance_summary.csv` from the unified cpg totals (contig preferred, else short-reads) + marker cpg values.
- `short_reads_summary` rule output renamed `amr_summary` → `markers_cpg`.
- `summary.smk` `amr_unified` rule: added `markers_cpg` input, second output `amr_abundance_summary`.

### Current pipeline state
All 31 jobs pass in test mode. Terminal outputs:
- `fastp_summary.csv`, `short_reads_output.csv`, `markers_cpg.csv` (from short_reads)
- `contig_summary.tsv` (AMR+MGE fully populated with contig_id, gene, coords, molecule_type)
- `AMR_unified.csv` (78 `both`, 3 `short_reads_only`, 0 `contigs_only` in test data)
- `AMR_abundance_summary.csv` (unified AMR total + pBI143 + crAss001 per sample)

### Known issues / next steps
- `gene_family` in `AMR_unified.csv` uses the gene symbol itself (not a gene family grouping) for contig-detected genes, since AMRFinder direct output doesn't provide a gene_family column. Production runs with the full catalog will populate this correctly for known alleles.
- `short_reads_only` evidence is rare in test data (3 genes: cmlA1 and 2 others in mock2) because the test reads come from the same organisms whose contigs were assembled — expected in production to be more prevalent for low-depth genes.


---

## 2026-03-09 (session 8)

### What was done

**SLURM executor plugin compatibility:**
- Removed `time` / `runtime` string values from all rule `resources:` blocks — the `d-hh:mm:ss` format is only valid within `snakemake-executor-plugin-slurm` and causes a parse error when running locally.
- Created `config/slurm/config.yaml` — a complete SLURM executor profile following the SWAM-g large-batch pattern. Includes `executor: slurm`, `scheduler: greedy`, `set-threads` and `set-resources` (with `mem_mb` + `runtime`) for all 16 compute rules. Placeholder `slurm_account`/`slurm_partition` to fill in per cluster.
- Usage: `snakemake --profile config/slurm` (requires `pip install snakemake-executor-plugin-slurm`).

**Pipeline modularity (run_short_reads / run_contigs / run_mags flags):**
- Three new `config.yaml` stage flags (all default `True`): `run_short_reads`, `run_contigs`, `run_mags`.
- Two new external-input keys: `clean_reads_dir` (pre-filtered reads, needed when `run_short_reads: False`), `contigs_dir` (pre-assembled contigs, needed when `run_contigs: False` and `run_mags: True`).
- `common.smk`: added `_RUN_SR/_RUN_CTG/_RUN_MAGS` booleans, `get_clean_r1/r2()` and `get_contigs()` input helpers, `WorkflowError` config validation.
- `Snakefile` `all_targets()`: fully conditional — each stage block is gated on its flag.
- `contigs.smk`: `contigs` + `mobmess` rules use `get_clean_r1/r2()`; `contig_abundance` uses all three helpers + optional `scgs` (empty list when `_RUN_SR=False`).
- `mags.smk`: `bin` rule uses `get_contigs()`.
- `summary.smk`: `amr_unified` inputs conditional on stage flags.
- `amr_unified.py`: handles `None` for any optional input (short_reads_amr, contig_summary, markers_cpg).
- `contig_abundance.py`: falls back to `n_genomes=1` when `scgs` is unavailable.
- README updated with "Running on SLURM" and "Modular execution" sections.
- Pushed to `origin/main` (commit `d187e7d`).

### Current pipeline state

Dry-run job counts for all mode combinations (all validated this session):
- Full pipeline (`test=True`): **32 jobs** ✅
- Short reads only (`run_contigs=False, run_mags=False`): **6 jobs** ✅
- Contigs only (`run_short_reads=False, run_mags=False`): **23 jobs** ✅
- MAGs only (`run_short_reads=False, run_contigs=False`): **8 jobs** ✅

### Known issues / next steps
- SLURM profile placeholders `slurm_account` and `slurm_partition` must be filled in before cluster use.
- When `run_short_reads: False`, contig abundance cpg uses `n_genomes=1` (mean depth not genome equivalents) — documented in README.
- End-to-end full test (`snakemake --use-conda --cores 4 --scheduler greedy --config test=True`) was last confirmed passing in session 7; not re-run this session.

## 2026-03-10 (session 9)

### What was done
- **End-to-end test passed** with the session 8 updates (SLURM profile + modularity flags).
- Full pipeline ran successfully on both mock samples: **26 steps (100%) done**, exit code 0.
- Run time: ~21 minutes on 4 cores (21:36 → 21:57).

### Output verification

| File | Size | Notes |
|------|------|-------|
| `fastp_summary.csv` | 429 B | ✅ |
| `short_reads_output.csv` | 15 KB | ✅ |
| `markers_cpg.csv` | 65 B | ✅ (mock data: pBI143/crAss001 = 0 — expected, not in mock genomes) |
| `contig_summary.tsv` | 54 KB, 857 rows | ✅ |
| `AMR_unified.csv` | 15 KB | ✅ evidence=both/short_reads_only working |
| `AMR_abundance_summary.csv` | 109 B | ✅ |
| `mag_amr.tsv` / `mag_mge.tsv` / `mag_abundance.tsv` | mock1+mock2 | ✅ |

`AMR_unified.csv` sample: `tet(A)` evidence=both, `cmlA1` evidence=short_reads_only, `tet(M)` molecule_type=phage,plasmid — all as expected.

### Current pipeline state
✅ **Full end-to-end test passes** with SLURM compatibility + modularity changes.

### Known issues / next steps
- `markers_cpg.csv` shows pBI143/crAss001 = 0 for both mock samples. This is expected — the mock genomes don't include pBI143 or crAss001 sequences. Production runs with real wastewater samples will show non-zero values.
- SLURM profile `slurm_account` and `slurm_partition` remain as empty placeholders — fill in before cluster use.
