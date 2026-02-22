# workflow/rules/mags.smk
#
# MAG (Metagenome-Assembled Genome) binning and characterisation.
#
# Rules in dependency order:
#   bin            – MetaBAT2 binning from contig BAM
#   mag_prodigal   – Prodigal gene prediction per bin
#   mag_amr        – AMRFinderPlus on each bin's predicted proteins
#   mag_mge        – MobileElementFinder on each bin's contigs
#   mag_taxonomy   – GTDB-tk classify_wf (skipped when config skip_gtdbtk=True)
#   mag_metabolism – METABOLIC (skipped when config skip_metabolic=True)
#   mag_abundance  – CoverM trimmed_mean per bin

import os

# ---------------------------------------------------------------------------
#   Binning with MetaBAT2
# ---------------------------------------------------------------------------

rule bin:
    input:
        contigs  = os.path.join(output_dir, "data", "megahit", "{sample}.contigs.fa"),
        bam      = os.path.join(output_dir, "data", "contig_abundance", "{sample}_contigs.bam"),
        bai      = os.path.join(output_dir, "data", "contig_abundance", "{sample}_contigs.bam.bai")
    output:
        depth    = os.path.join(output_dir, "data", "bins", "{sample}", "contig_depth.txt"),
        done     = os.path.join(output_dir, "data", "bins", "{sample}", ".binning.done")
    params:
        bin_dir    = os.path.join(output_dir, "data", "bins", "{sample}", "bins"),
        min_contig = lambda wc: res(2500, 1500),
        min_cls    = lambda wc: res(200000, 50000)
    resources:
        mem_mb  = lambda wc: res(150000, 6000),
        threads = lambda wc: res(32, 4),
        time    = "1-00:00:00"
    conda: "../envs/mags.yaml"
    shell:
        """
        mkdir -p {params.bin_dir}

        # Generate contig depth profile (reused by CoverM via the BAM,
        # but we write the MetaBAT2 depth file here for binning)
        jgi_summarize_bam_contig_depths \
            --outputDepth {output.depth} \
            {input.bam}

        metabat2 \
            -i {input.contigs} \
            -a {output.depth} \
            -o {params.bin_dir}/bin \
            -m {params.min_contig} \
            --minClsSize {params.min_cls} \
            -t {resources.threads} \
            --seed 42

        # Strip MetaBAT2's tab-separated depth info from bin FASTA headers.
        # CoverM genome mode matches sequence names from the BAM (which are
        # truncated at the first whitespace) against FASTA header names.
        # MetaBAT2 appends "\ttotal_depth=...\tsample_depths=..." to each
        # header; removing it ensures names match the BAM reference names.
        for fa in {params.bin_dir}/*.fa; do
            [ -f "$fa" ] || continue
            sed -i 's/\t.*//' "$fa"
        done

        # Write an empty bin placeholder if no bins were produced
        # (common with small test assemblies)
        if [ -z "$(ls -A {params.bin_dir})" ]; then
            echo "No bins produced for {wildcards.sample}" >&2
        fi

        touch {output.done}
        """

# ---------------------------------------------------------------------------
#   Per-bin: Prodigal gene calling  (needed by AMR + MGE rules below)
# ---------------------------------------------------------------------------

checkpoint mag_prodigal:
    """
    A checkpoint is used here because the number of bins is unknown until
    MetaBAT2 finishes. The downstream rules use aggregate() to collect
    per-bin outputs after this checkpoint resolves.
    """
    input:
        done = os.path.join(output_dir, "data", "bins", "{sample}", ".binning.done")
    output:
        done = os.path.join(output_dir, "data", "bins", "{sample}", ".prodigal.done")
    params:
        bin_dir = os.path.join(output_dir, "data", "bins", "{sample}", "bins"),
        prd_dir = os.path.join(output_dir, "data", "bins", "{sample}", "prodigal")
    resources:
        mem_mb  = lambda wc: res(16000, 4000),
        threads = 1,
        time    = "0-06:00:00"
    conda: "../envs/mags.yaml"
    shell:
        """
        mkdir -p {params.prd_dir}
        for fa in {params.bin_dir}/*.fa; do
            [ -f "$fa" ] || continue
            bin=$(basename "$fa" .fa)
            prodigal \
                -i "$fa" \
                -a {params.prd_dir}/$bin.faa \
                -f gff \
                -o {params.prd_dir}/$bin.gff \
                -p meta -q
        done
        touch {output.done}
        """

# ---------------------------------------------------------------------------
#   Per-sample: AMRFinderPlus on all bins (called after mag_prodigal checkpoint)
# ---------------------------------------------------------------------------

rule mag_amr:
    input:
        done = os.path.join(output_dir, "data", "bins", "{sample}", ".prodigal.done")
    output:
        tsv  = os.path.join(output_dir, "data", "bins", "{sample}", "mag_amr.tsv")
    params:
        bin_dir = os.path.join(output_dir, "data", "bins", "{sample}", "bins"),
        prd_dir = os.path.join(output_dir, "data", "bins", "{sample}", "prodigal")
    resources:
        mem_mb  = lambda wc: res(16000, 4000),
        threads = lambda wc: res(16, 4),
        time    = "0-06:00:00"
    conda: "../envs/mags.yaml"
    shell:
        """
        # Write a header-only TSV then append per-bin results
        TMP_AMR=$(mktemp)
        HEADER_WRITTEN=0
        for fa in {params.bin_dir}/*.fa; do
            [ -f "$fa" ] || continue
            bin=$(basename "$fa" .fa)
            faa="{params.prd_dir}/$bin.faa"
            gff="{params.prd_dir}/$bin.gff"
            [ -f "$faa" ] || continue
            amrfinder \
                -n "$fa" \
                -p "$faa" \
                -g "$gff" \
                --threads {resources.threads} \
                --annotation_format prodigal \
                -o "$TMP_AMR.bin"
            if [ -s "$TMP_AMR.bin" ]; then
                if [ $HEADER_WRITTEN -eq 0 ]; then
                    head -1 "$TMP_AMR.bin" > "$TMP_AMR"
                    HEADER_WRITTEN=1
                fi
                tail -n +2 "$TMP_AMR.bin" | awk -v bin="$bin" 'OFS="\\t"{{print bin,$0}}' >> "$TMP_AMR"
            fi
        done
        mv "$TMP_AMR" {output.tsv}
        [ -f {output.tsv} ] || printf "bin_id\\n" > {output.tsv}
        """

# ---------------------------------------------------------------------------
#   Per-sample: MobileElementFinder on all bins
# ---------------------------------------------------------------------------

rule mag_mge:
    input:
        done     = os.path.join(output_dir, "data", "bins", "{sample}", ".binning.done"),
        mef_done = os.path.join(output_dir, "data", "mge_contigs", ".mef_install.done")
    output:
        tsv  = os.path.join(output_dir, "data", "bins", "{sample}", "mag_mge.tsv")
    params:
        bin_dir = os.path.join(output_dir, "data", "bins", "{sample}", "bins"),
        tmp_dir = os.path.join(output_dir, "data", "bins", "{sample}", "mge_tmp")
    resources:
        mem_mb  = lambda wc: res(16000, 4000),
        threads = lambda wc: res(16, 4),
        time    = "0-06:00:00"
    conda: "../envs/mge.yaml"
    shell:
        """
        mkdir -p {params.tmp_dir}
        TMP_MGE=$(mktemp)
        HEADER_WRITTEN=0
        for fa in {params.bin_dir}/*.fa; do
            [ -f "$fa" ] || continue
            bin=$(basename "$fa" .fa)
            PREFIX={params.tmp_dir}/$bin
            mkdir -p {params.tmp_dir}/$bin
            mefinder find \
                --contig "$fa" \
                --temp-dir {params.tmp_dir}/$bin \
                --threads {resources.threads} \
                "${{PREFIX}}_out" || true
            csv="${{PREFIX}}_out.csv"
            if [ -s "$csv" ]; then
                if [ $HEADER_WRITTEN -eq 0 ]; then
                    head -1 "$csv" > "$TMP_MGE"
                    HEADER_WRITTEN=1
                fi
                tail -n +2 "$csv" | awk -v bin="$bin" 'OFS="\\t"{{print bin,$0}}' >> "$TMP_MGE"
            fi
        done
        rm -rf {params.tmp_dir}
        mv "$TMP_MGE" {output.tsv}
        [ -f {output.tsv} ] || printf "bin_id\\n" > {output.tsv}
        """

# ---------------------------------------------------------------------------
#   Per-sample: GTDB-tk taxonomy (skipped when config skip_gtdbtk=True)
# ---------------------------------------------------------------------------

rule mag_taxonomy:
    input:
        done = os.path.join(output_dir, "data", "bins", "{sample}", ".binning.done")
    output:
        tsv  = os.path.join(output_dir, "data", "bins", "{sample}", "gtdbtk.summary.tsv")
    params:
        bin_dir  = os.path.join(output_dir, "data", "bins", "{sample}", "bins"),
        out_dir  = os.path.join(output_dir, "data", "bins", "{sample}", "gtdbtk"),
        gtdb_db  = config.get("gtdbtk_db", ""),
        skip     = config.get("skip_gtdbtk", False)
    resources:
        mem_mb  = lambda wc: res(150000, 8000),
        threads = lambda wc: res(64, 4),
        time    = "1-12:00:00"
    conda: "../envs/gtdbtk.yaml"
    shell:
        """
        if [ "{params.skip}" = "True" ]; then
            echo "GTDB-tk skipped (skip_gtdbtk=True)" >&2
            printf "user_genome\\tclassification\\n" > {output.tsv}
        else
            export GTDBTK_DATA_PATH="{params.gtdb_db}"
            mkdir -p {params.out_dir}
            gtdbtk classify_wf \
                --genome_dir {params.bin_dir} \
                --extension fa \
                --out_dir {params.out_dir} \
                --cpus {resources.threads} \
                --skip_ani_screen
            # Combine archaea + bacteria summaries
            TMP=$(mktemp)
            for f in {params.out_dir}/classify/gtdbtk.*.summary.tsv; do
                [ -f "$f" ] || continue
                head -1 "$f" > "$TMP" && tail -n +2 "$f" >> "$TMP"
            done
            mv "$TMP" {output.tsv}
            [ -f {output.tsv} ] || printf "user_genome\\tclassification\\n" > {output.tsv}
        fi
        """

# ---------------------------------------------------------------------------
#   Per-sample: METABOLIC (skipped when config skip_metabolic=True)
# ---------------------------------------------------------------------------

rule mag_metabolism:
    input:
        done = os.path.join(output_dir, "data", "bins", "{sample}", ".binning.done")
    output:
        done = os.path.join(output_dir, "data", "bins", "{sample}", ".metabolic.done")
    params:
        bin_dir      = os.path.join(output_dir, "data", "bins", "{sample}", "bins"),
        out_dir      = os.path.join(output_dir, "data", "bins", "{sample}", "metabolic"),
        metabolic_dir = config.get("metabolic_dir", ""),
        skip         = config.get("skip_metabolic", False)
    resources:
        mem_mb  = lambda wc: res(150000, 8000),
        threads = lambda wc: res(64, 4),
        time    = "1-00:00:00"
    conda: "../envs/metabolic.yaml"
    shell:
        """
        if [ "{params.skip}" = "True" ]; then
            echo "METABOLIC skipped (skip_metabolic=True)" >&2
        else
            mkdir -p {params.out_dir}
            perl {params.metabolic_dir}/METABOLIC-G.pl \
                -in {params.bin_dir} \
                -r {params.out_dir}/reads_list.txt \
                -o {params.out_dir} \
                -t {resources.threads}
        fi
        touch {output.done}
        """

# ---------------------------------------------------------------------------
#   Per-sample: CheckM2 MAG quality assessment (completeness + contamination)
# ---------------------------------------------------------------------------

rule mag_qc:
    input:
        done = os.path.join(output_dir, "data", "bins", "{sample}", ".binning.done")
    output:
        tsv  = os.path.join(output_dir, "data", "bins", "{sample}", "checkm2_quality.tsv")
    params:
        bin_dir = os.path.join(output_dir, "data", "bins", "{sample}", "bins"),
        out_dir = os.path.join(output_dir, "data", "bins", "{sample}", "checkm2"),
        db_path = config.get("checkm2_db", ""),
        skip    = config.get("skip_checkm2", False)
    resources:
        mem_mb  = lambda wc: res(32000, 4000),
        threads = lambda wc: res(16, 4),
        time    = "0-08:00:00"
    conda: "../envs/checkm2.yaml"
    shell:
        """
        if [ "{params.skip}" = "True" ]; then
            echo "CheckM2 skipped (skip_checkm2=True)" >&2
            printf "Name\\tCompleteness\\tContamination\\n" > {output.tsv}
        else
            n_bins=$(find {params.bin_dir} -name "*.fa" 2>/dev/null | wc -l)
            if [ "$n_bins" -eq 0 ]; then
                printf "Name\\tCompleteness\\tContamination\\n" > {output.tsv}
            else
                mkdir -p {params.out_dir}
                checkm2 predict \
                    --input {params.bin_dir} \
                    --output-directory {params.out_dir} \
                    --database_path {params.db_path} \
                    --threads {resources.threads} \
                    -x fa
                cp {params.out_dir}/quality_report.tsv {output.tsv}
            fi
        fi
        """


rule mag_abundance:
    input:
        bam      = os.path.join(output_dir, "data", "contig_abundance", "{sample}_contigs.bam"),
        bai      = os.path.join(output_dir, "data", "contig_abundance", "{sample}_contigs.bam.bai"),
        done     = os.path.join(output_dir, "data", "bins", "{sample}", ".binning.done"),
        depth    = os.path.join(output_dir, "data", "bins", "{sample}", "contig_depth.txt")
    output:
        tsv = os.path.join(output_dir, "data", "bins", "{sample}", "mag_abundance.tsv")
    params:
        bin_dir = os.path.join(output_dir, "data", "bins", "{sample}", "bins")
    resources:
        mem_mb  = lambda wc: res(32000, 4000),
        threads = lambda wc: res(16, 4),
        time    = "0-04:00:00"
    conda: "../envs/mags.yaml"
    shell:
        """
        # Count bins; if none exist write empty output and exit cleanly
        n_bins=$(find {params.bin_dir} -name "*.fa" 2>/dev/null | wc -l)
        if [ "$n_bins" -eq 0 ]; then
            printf "Genome\\ttrimmed_mean\\n" > {output.tsv}
        else
            coverm genome \
                --bam-files {input.bam} \
                --genome-fasta-files {params.bin_dir}/*.fa \
                --methods trimmed_mean \
                --min-covered-fraction 0.1 \
                --threads {resources.threads} \
                --output-file {output.tsv}
        fi
        """
