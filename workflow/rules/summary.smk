# workflow/rules/summary.smk
#
# Cross-stage summary rules for MAG-stage outputs.

# ---------------------------------------------------------------------------
#   MAG summary table
#   Joins per-sample mag_amr, mag_mge, mag_abundance (and optionally
#   checkm2 quality + GTDB-Tk taxonomy) into a single TSV.
#   One row per (sample, bin).
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
#   AMR risk score table
#   Computes additive and multiplicative AMR risk scores per sample.
#   Sources: short_reads_output.csv (gene cpg + class/subclass),
#            markers_cpg.csv (pBI143 / crAss001 — Exposure component),
#            contig_summary.tsv (mobility + host taxonomy; optional).
#   Produces AMR_abundance_summary.csv with per-sample risk columns.
# ---------------------------------------------------------------------------

rule amr_risk_score:
    input:
        short_reads    = SHORT_READS_OUTPUT,
        markers_cpg    = MARKERS_CPG_OUTPUT,
        contig_summary = ([CONTIG_SUMMARY_OUTPUT] if _RUN_CTG else []),
    output:
        AMR_ABUNDANCE_SUMMARY_OUTPUT
    params:
        samples = samples,
    conda: "../envs/contigs.yaml"
    log:
        rule_log("amr_risk_score")
    script:
        "../scripts/amr_risk_score.py"


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
        tsv = MAG_SUMMARY_OUTPUT
    params:
        samples      = samples,
        output_dir   = output_dir,
        skip_checkm2 = lambda wc: config.get("skip_checkm2", False),
        skip_gtdbtk  = lambda wc: config.get("skip_gtdbtk",  False),
    conda: "../envs/contigs.yaml"
    log:
        rule_log("mag_summary")
    script:
        "../scripts/mag_summary.py"
