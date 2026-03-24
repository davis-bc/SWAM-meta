# workflow/rules/summary.smk
#
# Cross-stage summary rules that depend on outputs from both
# short_reads.smk and contigs.smk, and MAG-stage outputs.

# ---------------------------------------------------------------------------
#   Unified AMR abundance table
#   Merges KMA short-read detections with contig-level AMRFinderPlus
#   annotations using cpg as the common normalization unit.
# ---------------------------------------------------------------------------

rule amr_unified:
    input:
        short_reads_amr = lambda wc: (
            os.path.join(output_dir, "data", "QAQC", "short_reads_output.csv") if _RUN_SR else []
        ),
        contig_summary  = lambda wc: (
            os.path.join(output_dir, "contig_summary.tsv") if _RUN_CTG else []
        ),
        catalog         = os.path.join(_SR_DBS_DIR, "ReferenceGeneCatalog.txt"),
        markers_cpg     = lambda wc: (
            os.path.join(output_dir, "data", "QAQC", "markers_cpg.csv") if _RUN_SR else []
        )
    output:
        csv                    = os.path.join(output_dir, "AMR_unified.csv"),
        amr_abundance_summary  = os.path.join(output_dir, "AMR_abundance_summary.csv")
    params:
        contig_amr_dir = os.path.join(output_dir, "data", "amr_contigs")
    conda: "../envs/contigs.yaml"
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "amr_unified.log")
    script:
        "../scripts/amr_unified.py"

# ---------------------------------------------------------------------------
#   MAG summary table
#   Joins per-sample mag_amr, mag_mge, mag_abundance (and optionally
#   checkm2 quality + GTDB-Tk taxonomy) into a single TSV.
#   One row per (sample, bin).
# ---------------------------------------------------------------------------

rule mag_summary:
    input:
        amr_files   = expand(
            os.path.join(output_dir, "data", "bins", "{sample}", "mag_amr.tsv"),
            sample=samples
        ),
        mge_files   = expand(
            os.path.join(output_dir, "data", "bins", "{sample}", "mag_mge.tsv"),
            sample=samples
        ),
        abund_files = expand(
            os.path.join(output_dir, "data", "bins", "{sample}", "mag_abundance.tsv"),
            sample=samples
        ),
    output:
        tsv = os.path.join(output_dir, "mag_summary.tsv")
    params:
        samples      = samples,
        output_dir   = output_dir,
        skip_checkm2 = lambda wc: config.get("skip_checkm2", False),
        skip_gtdbtk  = lambda wc: config.get("skip_gtdbtk",  False),
    conda: "../envs/contigs.yaml"
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "mag_summary.log")
    script:
        "../scripts/mag_summary.py"

