# workflow/rules/mags.smk
#
# MAG (Metagenome-Assembled Genome) binning and characterisation.
#
# Rules in dependency order:
#   init_checkm2_db – one-time CheckM2 database download
#   bin            – MetaBAT2 binning from contig BAM
#   mag_prodigal   – Prodigal gene prediction per bin
#   mag_amr        – AMRFinderPlus on each bin's predicted proteins
#   mag_mge        – MobileElementFinder on each bin's contigs
#   mag_taxonomy   – GTDB-tk classify_wf (skipped when config skip_gtdbtk=True)
#   mag_metabolism – METABOLIC (skipped when config skip_metabolic=True)
#   mag_qc         – CheckM2 quality assessment (skipped when config skip_checkm2=True)
#   mag_abundance  – CoverM trimmed_mean per bin

import os

# ---------------------------------------------------------------------------
#   One-time: download CheckM2 reference database into _DBS_DIR
# ---------------------------------------------------------------------------

rule init_checkm2_db:
    output:
        done = os.path.join(_DBS_DIR, ".checkm2.done")
    params:
        db_dir    = os.path.join(_DBS_DIR, "checkm2"),
        db_path   = _CHECKM2_DB,
        skip      = config.get("skip_checkm2", False),
        test_mode = _TEST
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "init_checkm2_db.log")
    conda: "../envs/checkm2.yaml"
    shell:
        """
        mkdir -p {params.db_dir}
        if [ "{params.skip}" = "True" ]; then
            echo "CheckM2: database setup skipped (skip_checkm2=True)"
        elif [ "{params.test_mode}" = "True" ]; then
            echo "CheckM2: test mode — skipping database download"
        elif [ -f "{params.db_path}" ]; then
            echo "CheckM2: database already exists, skipping download"
        else
            echo "CheckM2: downloading database (~3 GB)..."
            checkm2 database --download --path {params.db_dir} >> {log} 2>&1
            echo "CheckM2: database ready"
        fi
        touch {output.done}
        """

# ---------------------------------------------------------------------------
#   Binning with MetaBAT2
# ---------------------------------------------------------------------------

rule bin:
    input:
        contigs  = get_contigs,
        bam      = os.path.join(output_dir, "data", "contig_abundance", "{sample}_contigs.bam"),
        bai      = os.path.join(output_dir, "data", "contig_abundance", "{sample}_contigs.bam.bai")
    output:
        depth    = os.path.join(output_dir, "data", "bins", "{sample}", "contig_depth.txt"),
        done     = os.path.join(output_dir, "data", "bins", "{sample}", ".binning.done")
    params:
        bin_dir    = os.path.join(output_dir, "data", "bins", "{sample}", "bins"),
        min_contig = lambda wc: res(2500, 1500),
        min_cls    = lambda wc: res(200000, 50000)
    threads: lambda wc: res(32, 4)
    resources:
        mem_mb  = lambda wc: res(150000, 6000),
        threads = lambda wc: res(32, 4),
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "{sample}.bin.log")
    conda: "../envs/mags.yaml"
    shell:
        """
        mkdir -p {params.bin_dir}

        echo "[{wildcards.sample}] MetaBAT2: profiling contig depths..."
        jgi_summarize_bam_contig_depths \
            --outputDepth {output.depth} \
            {input.bam} >> {log} 2>&1

        echo "[{wildcards.sample}] MetaBAT2: binning contigs..."
        metabat2 \
            -i {input.contigs} \
            -a {output.depth} \
            -o {params.bin_dir}/bin \
            -m {params.min_contig} \
            --minClsSize {params.min_cls} \
            -t {resources.threads} \
            --seed 42 >> {log} 2>&1

        # Strip MetaBAT2's tab-separated depth info from bin FASTA headers.
        for fa in {params.bin_dir}/*.fa; do
            [ -f "$fa" ] || continue
            sed -i 's/\t.*//' "$fa"
        done

        if [ -z "$(ls -A {params.bin_dir})" ]; then
            echo "[{wildcards.sample}] MetaBAT2: no bins produced"
        else
            echo "[{wildcards.sample}] MetaBAT2: done"
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
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "{sample}.mag_prodigal.log")
    conda: "../envs/mags.yaml"
    shell:
        """
        mkdir -p {params.prd_dir}
        echo "[{wildcards.sample}] Prodigal: predicting genes in MAG bins..."
        for fa in {params.bin_dir}/*.fa; do
            [ -f "$fa" ] || continue
            bin=$(basename "$fa" .fa)
            prodigal \
                -i "$fa" \
                -a {params.prd_dir}/$bin.faa \
                -f gff \
                -o {params.prd_dir}/$bin.gff \
                -p meta -q >> {log} 2>&1
        done
        touch {output.done}
        echo "[{wildcards.sample}] Prodigal: done"
        """

# ---------------------------------------------------------------------------
#   Per-sample: AMRFinderPlus on all bins (called after mag_prodigal checkpoint)
# ---------------------------------------------------------------------------

rule mag_amr:
    input:
        done        = os.path.join(output_dir, "data", "bins", "{sample}", ".prodigal.done"),
        afp_db_done = os.path.join(_DBS_DIR, ".amrfinder_db.done")
    output:
        tsv  = os.path.join(output_dir, "data", "bins", "{sample}", "mag_amr.tsv")
    params:
        bin_dir    = os.path.join(output_dir, "data", "bins", "{sample}", "bins"),
        prd_dir    = os.path.join(output_dir, "data", "bins", "{sample}", "prodigal"),
        afp_db_dir = _AFP_DB_DIR
    threads: lambda wc: res(16, 4)
    resources:
        mem_mb  = lambda wc: res(16000, 4000),
        threads = lambda wc: res(16, 4),
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "{sample}.mag_amr.log")
    conda: "../envs/mags.yaml"
    shell:
        """
        echo "[{wildcards.sample}] AMRFinderPlus: annotating MAG bins..."
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
                --database {params.afp_db_dir} \
                --threads {resources.threads} \
                --annotation_format prodigal \
                -o "$TMP_AMR.bin" >> {log} 2>&1 || true
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
        echo "[{wildcards.sample}] AMRFinderPlus: done"
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
    threads: lambda wc: res(16, 4)
    resources:
        mem_mb  = lambda wc: res(16000, 4000),
        threads = lambda wc: res(16, 4),
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "{sample}.mag_mge.log")
    conda: "../envs/contigs.yaml"
    shell:
        """
        mkdir -p {params.tmp_dir}
        echo "[{wildcards.sample}] MobileElementFinder: annotating MAG MGEs..."
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
                "${{PREFIX}}_out" >> {log} 2>&1 || true
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
        echo "[{wildcards.sample}] MobileElementFinder: done"
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
    threads: lambda wc: res(64, 4)
    resources:
        mem_mb  = lambda wc: res(150000, 8000),
        threads = lambda wc: res(64, 4),
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "{sample}.mag_taxonomy.log")
    conda: "../envs/gtdbtk.yaml"
    shell:
        """
        if [ "{params.skip}" = "True" ]; then
            echo "[{wildcards.sample}] GTDB-tk: skipped (skip_gtdbtk=True)"
            printf "user_genome\\tclassification\\n" > {output.tsv}
        else
            echo "[{wildcards.sample}] GTDB-tk: classifying MAG taxonomy..."
            export GTDBTK_DATA_PATH="{params.gtdb_db}"
            mkdir -p {params.out_dir}
            gtdbtk classify_wf \
                --genome_dir {params.bin_dir} \
                --extension fa \
                --out_dir {params.out_dir} \
                --cpus {resources.threads} \
                --skip_ani_screen >> {log} 2>&1
            # Combine archaea + bacteria summaries
            TMP=$(mktemp)
            for f in {params.out_dir}/classify/gtdbtk.*.summary.tsv; do
                [ -f "$f" ] || continue
                head -1 "$f" > "$TMP" && tail -n +2 "$f" >> "$TMP"
            done
            mv "$TMP" {output.tsv}
            [ -f {output.tsv} ] || printf "user_genome\\tclassification\\n" > {output.tsv}
            echo "[{wildcards.sample}] GTDB-tk: done"
        fi
        """

# ---------------------------------------------------------------------------
#   Per-sample: METABOLIC (skipped when config skip_metabolic=True)
# ---------------------------------------------------------------------------

rule mag_metabolism:
    input:
        done          = os.path.join(output_dir, "data", "bins", "{sample}", ".binning.done"),
        metabolic_done = os.path.join(_DBS_DIR, ".metabolic.done.txt")
    output:
        done = os.path.join(output_dir, "data", "bins", "{sample}", ".metabolic.done")
    params:
        bin_dir       = os.path.join(output_dir, "data", "bins", "{sample}", "bins"),
        out_dir       = os.path.join(output_dir, "data", "bins", "{sample}", "metabolic"),
        metabolic_dir = _METABOLIC_DIR,
        skip          = config.get("skip_metabolic", False)
    threads: lambda wc: res(64, 4)
    resources:
        mem_mb  = lambda wc: res(150000, 8000),
        threads = lambda wc: res(64, 4),
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "{sample}.mag_metabolism.log")
    conda: "../envs/metabolic.yaml"
    shell:
        """
        if [ "{params.skip}" = "True" ]; then
            echo "[{wildcards.sample}] METABOLIC: skipped (skip_metabolic=True)"
        else
            echo "[{wildcards.sample}] METABOLIC: profiling metabolism..."
            mkdir -p {params.out_dir}
            perl {params.metabolic_dir}/METABOLIC-G.pl \
                -in {params.bin_dir} \
                -r {params.out_dir}/reads_list.txt \
                -o {params.out_dir} \
                -t {resources.threads} >> {log} 2>&1
            echo "[{wildcards.sample}] METABOLIC: done"
        fi
        touch {output.done}
        """

# ---------------------------------------------------------------------------
#   Per-sample: CheckM2 MAG quality assessment (completeness + contamination)
# ---------------------------------------------------------------------------

rule mag_qc:
    input:
        done         = os.path.join(output_dir, "data", "bins", "{sample}", ".binning.done"),
        checkm2_done = os.path.join(_DBS_DIR, ".checkm2.done")
    output:
        tsv  = os.path.join(output_dir, "data", "bins", "{sample}", "checkm2_quality.tsv")
    params:
        bin_dir = os.path.join(output_dir, "data", "bins", "{sample}", "bins"),
        out_dir = os.path.join(output_dir, "data", "bins", "{sample}", "checkm2"),
        db_path = _CHECKM2_DB,
        skip    = config.get("skip_checkm2", False)
    threads: lambda wc: res(16, 4)
    resources:
        mem_mb  = lambda wc: res(32000, 4000),
        threads = lambda wc: res(16, 4),
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "{sample}.mag_qc.log")
    conda: "../envs/checkm2.yaml"
    shell:
        """
        if [ "{params.skip}" = "True" ]; then
            echo "[{wildcards.sample}] CheckM2: skipped (skip_checkm2=True)"
            printf "Name\\tCompleteness\\tContamination\\n" > {output.tsv}
        else
            n_bins=$(find {params.bin_dir} -name "*.fa" 2>/dev/null | wc -l)
            if [ "$n_bins" -eq 0 ]; then
                echo "[{wildcards.sample}] CheckM2: no bins to assess"
                printf "Name\\tCompleteness\\tContamination\\n" > {output.tsv}
            else
                echo "[{wildcards.sample}] CheckM2: assessing MAG quality..."
                mkdir -p {params.out_dir}
                checkm2 predict \
                    --input {params.bin_dir} \
                    --output-directory {params.out_dir} \
                    --database_path {params.db_path} \
                    --threads {resources.threads} \
                    -x fa >> {log} 2>&1
                cp {params.out_dir}/quality_report.tsv {output.tsv}
                echo "[{wildcards.sample}] CheckM2: done"
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
    threads: lambda wc: res(16, 4)
    resources:
        mem_mb  = lambda wc: res(32000, 4000),
        threads = lambda wc: res(16, 4),
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "{sample}.mag_abundance.log")
    conda: "../envs/shortreads.yaml"
    shell:
        """
        # Count bins; if none exist write empty output and exit cleanly
        n_bins=$(find {params.bin_dir} -name "*.fa" 2>/dev/null | wc -l)
        if [ "$n_bins" -eq 0 ]; then
            echo "[{wildcards.sample}] CoverM: no bins found, skipping"
            printf "Genome\\ttrimmed_mean\\n" > {output.tsv}
        else
            echo "[{wildcards.sample}] CoverM: computing MAG abundance..."
            coverm genome \
                --bam-files {input.bam} \
                --genome-fasta-files {params.bin_dir}/*.fa \
                --methods trimmed_mean \
                --min-covered-fraction 0.1 \
                --threads {resources.threads} \
                --output-file {output.tsv} >> {log} 2>&1
            echo "[{wildcards.sample}] CoverM: done"
        fi
        """
