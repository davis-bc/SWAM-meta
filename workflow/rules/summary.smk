# workflow/rules/summary.smk
#
# Cross-stage summary rules that depend on outputs from both
# short_reads.smk and contigs.smk.

# ---------------------------------------------------------------------------
#   Unified AMR abundance table
#   Merges KMA short-read detections with contig-level AMRFinderPlus
#   annotations using cpg as the common normalization unit.
# ---------------------------------------------------------------------------

rule amr_unified:
    input:
        short_reads_amr = os.path.join(output_dir, "short_reads_output.csv"),
        contig_summary  = os.path.join(output_dir, "contig_summary.tsv"),
        catalog         = os.path.join(output_dir, "data", "alignments", "dbs", "ReferenceGeneCatalog.txt"),
        markers_cpg     = os.path.join(output_dir, "markers_cpg.csv")
    output:
        csv                    = os.path.join(output_dir, "AMR_unified.csv"),
        amr_abundance_summary  = os.path.join(output_dir, "AMR_abundance_summary.csv")
    params:
        contig_amr_dir = os.path.join(output_dir, "data", "amr_contigs")
    conda: "../envs/contigs.yaml"
    script:
        "../scripts/amr_unified.py"
