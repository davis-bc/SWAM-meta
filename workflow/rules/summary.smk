# workflow/rules/summary.smk
#
# Cross-stage summary rules for MAG-stage outputs.

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
