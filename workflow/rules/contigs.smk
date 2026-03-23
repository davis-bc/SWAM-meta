# ------------------------------------------
#       prepare geNomad database 
# ------------------------------------------

rule init_genomad:
    output:
        genomad_db   = os.path.join(_DBS_DIR, ".genomad.db.done.txt")
    params:
        dbs_dir = _DBS_DIR
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "init_genomad.log")
    conda: "../envs/genomad.yaml"
    shell:
        """
        mkdir -p {params.dbs_dir}
        if [ ! -d {params.dbs_dir}/genomad_db ]; then
        
        echo "geNomad: downloading database..."
        cd {params.dbs_dir}
        curl -L https://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz 2>> {log} | tar -xz 2>> {log}
        touch {output.genomad_db}
        echo "geNomad: database ready"
        
        else 
        
        echo "geNomad: database already exists, skipping download"
        touch {output.genomad_db}
        
        fi
        """
# ------------------------------------------
#       de novo assembly + gene annotation 
# ------------------------------------------

rule contigs:
    input:
        r1       = get_clean_r1,
        r2       = get_clean_r2
    output:
        megahit  = os.path.join(output_dir, "data", "megahit", "{sample}.contigs.fa")
    threads: lambda wc: res(32, 4)
    resources:
        mem_mb  = lambda wc: res(150000, 8000),
        threads = lambda wc: res(32, 4),
    benchmark:
        os.path.join(output_dir, "data", "QAQC", "benchmarks", "{sample}.contigs.txt")
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "{sample}.megahit.log")
    conda: "../envs/contigs.yaml"
    shell:
        """
        set -euo pipefail
        echo "[{wildcards.sample}] MEGAHIT: assembling contigs..."

        # create temporary directory in same parent so moves are on same fs (atomic)
        tmpdir="$(mktemp -d "$(dirname {output.megahit})/megahit_{wildcards.sample}.tmp.XXXXXX")"

        # Run MEGAHIT into the temporary directory
        megahit -1 {input.r1} -2 {input.r2} \
                -t {resources.threads} \
                --presets meta-large \
                --min-contig-len 1000 \
                -o "$tmpdir" \
                --out-prefix {wildcards.sample} \
                --continue \
                -f >> {log} 2>&1

        # Move the contigs file into the final megahit directory.
        mv -f "$tmpdir/{wildcards.sample}.contigs.fa" {output.megahit}

        # Clean up the rest of temporary files
        rm -rf "$tmpdir"
        
        # rename fasta headers to include sample information
        sed -E -i 's/>/>{wildcards.sample}-/g' {output.megahit}
        echo "[{wildcards.sample}] MEGAHIT: done"
        
        """
    
# ------------------------------------------------------------------------
#       classify contigs into chromosome, plasmid, and phage with geNomad
# ------------------------------------------------------------------------

rule genomad:
    input:
        contigs      = os.path.join(output_dir, "data", "megahit", "{sample}.contigs.fa"),
        genomad_db   = os.path.join(_DBS_DIR, ".genomad.db.done.txt")
    output:
        genomad      = directory(os.path.join(output_dir, "data", "genomad", "{sample}")),
        plas_contigs = os.path.join(output_dir, "data", "genomad", "{sample}", "{sample}.contigs_summary", "{sample}.contigs_plasmid.fna"),
        plas_summary = os.path.join(output_dir, "data", "genomad", "{sample}", "{sample}.contigs_summary", "{sample}.contigs_plasmid_summary.tsv"),
        vir_summary  = os.path.join(output_dir, "data", "genomad", "{sample}", "{sample}.contigs_summary", "{sample}.contigs_virus_summary.tsv")
    params:
        splits = config.get("genomad_splits", 1)
    threads: lambda wc: res(32, 4)
    resources:
        mem_mb  = lambda wc: res(150000, 8000),
        threads = lambda wc: res(32, 4),
    benchmark:
        os.path.join(output_dir, "data", "QAQC", "benchmarks", "{sample}.genomad.txt")
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "{sample}.genomad.log")
    conda: "../envs/genomad.yaml"
    shell:
        """
        echo "[{wildcards.sample}] geNomad: classifying contigs (chromosome/plasmid/phage)..."
        genomad end-to-end {input.contigs} {output.genomad} $(dirname {input.genomad_db})/genomad_db \
            --relaxed --cleanup --splits {params.splits} >> {log} 2>&1
        # touch missing output files (geNomad omits them when no hits are found)
        touch {output.plas_contigs} {output.plas_summary} {output.vir_summary}
        echo "[{wildcards.sample}] geNomad: done"
        """

# ----------------------------------------------------------------------------------
#       determine cirularity of plasmid contigs, infer plasmid systems with MobMess 
# ----------------------------------------------------------------------------------

rule mobmess:
    input:
        r1           = get_clean_r1,
        r2           = get_clean_r2,
        plas_contigs = os.path.join(output_dir, "data", "genomad", "{sample}",  "{sample}.contigs_summary", "{sample}.contigs_plasmid.fna")
    output:
        circular_txt = os.path.join(output_dir, "data", "mobmess", "{sample}.contigs_circular.txt"),
        contig_bam   = os.path.join(output_dir, "data", "mobmess", "{sample}.contigs_plasmid.bam"),
        mobmess      = os.path.join(output_dir, "data", "mobmess", "{sample}-mobmess_contigs.txt")
    conda: "../envs/contigs.yaml"
    threads: lambda wc: res(32, 4)
    resources:
        mem_mb  = lambda wc: res(150000, 8000),
        threads = lambda wc: res(32, 4),
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "{sample}.mobmess.log")
    shell:
        """
        mkdir -p $(dirname {output.contig_bam})
        # Skip processing if plasmid FASTA is empty (no plasmids found by geNomad)
        if [ ! -s {input.plas_contigs} ]; then
            echo "[{wildcards.sample}] MobMess: no plasmids found, skipping"
            touch {output.contig_bam} {output.circular_txt} {output.mobmess}
        else
            echo "[{wildcards.sample}] minimap2: mapping reads to plasmid contigs..."
            minimap2 -t {resources.threads} -ax sr {input.plas_contigs} {input.r1} {input.r2} 2>> {log} \
            | samtools sort -@4 -o {output.contig_bam} - 2>> {log} && samtools index {output.contig_bam} 2>> {log}
            echo "[{wildcards.sample}] minimap2: done"
            
            echo "[{wildcards.sample}] MobMess: detecting circular plasmid contigs..."
            workflow/scripts/detect_circular_contigs.py -b {output.contig_bam} -o {output.circular_txt} >> {log} 2>&1
            
            echo "[{wildcards.sample}] MobMess: inferring plasmid systems..."
            mobmess systems --sequences {input.plas_contigs} --complete {output.circular_txt} --output $(dirname {output.mobmess})/{wildcards.sample}-mobmess --threads {resources.threads} >> {log} 2>&1 || true
            # ensure output file always exists
            [ -f {output.mobmess} ] || touch {output.mobmess}
            echo "[{wildcards.sample}] MobMess: done"
        fi
        """

# ---------------------------------------------------------------------------
#   Build / validate MMseqs2 UniRef50 taxonomy database
#   Test mode: build from mini FASTA in test/dbs/uniref50/
#   Production: validate user-provided db path
# ---------------------------------------------------------------------------

rule init_mmseqs_db:
    output:
        done = os.path.join(_DBS_DIR, ".mmseqs_db.done")
    params:
        test_mode    = _TEST,
        db_prefix    = _UNIREF50_DB,
        dbs_dir      = _DBS_DIR,
        test_fasta   = os.path.join(_REPO, "test", "dbs", "uniref50", "uniref50_mini.fasta"),
        test_tax_dir = os.path.join(_REPO, "test", "dbs", "uniref50"),
        test_mapping = os.path.join(_REPO, "test", "dbs", "uniref50", "tax_mapping.tsv"),
    threads: lambda wc: res(16, 4)
    resources:
        mem_mb  = lambda wc: res(32000, 4000),
        threads = lambda wc: res(16, 4),
    conda: "../envs/contigs.yaml"
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "init_mmseqs_db.log")
    shell:
        """
        mkdir -p "$(dirname {params.db_prefix})"

        if [ "{params.test_mode}" = "True" ]; then
            echo "MMseqs2: building test taxonomy database..."
            mkdir -p {params.test_tax_dir}/tmp
            touch {params.test_tax_dir}/merged.dmp
            mmseqs createdb {params.test_fasta} {params.db_prefix} >> {log} 2>&1
            mmseqs createtaxdb {params.db_prefix} {params.test_tax_dir}/tmp \
                --ncbi-tax-dump {params.test_tax_dir} \
                --tax-mapping-file {params.test_mapping} >> {log} 2>&1
            mmseqs createindex {params.db_prefix} {params.test_tax_dir}/tmp --search-type 2 >> {log} 2>&1
            rm -rf {params.test_tax_dir}/tmp
            echo "MMseqs2: test database ready"
        else
            if [ ! -f "{params.db_prefix}.dbtype" ]; then
                UNIREF_FASTA="{params.dbs_dir}/uniref50/uniref50.fasta"
                TAX_DIR="{params.dbs_dir}/uniref50/taxdump"
                TMP_DIR="{params.dbs_dir}/uniref50/tmp"
                mkdir -p "$TAX_DIR" "$TMP_DIR"

                echo "MMseqs2: downloading UniRef50 FASTA (~9 GB)..."
                wget -q -O "$UNIREF_FASTA.gz" \
                    https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz >> {log} 2>&1
                echo "MMseqs2: decompressing UniRef50..."
                gunzip -f "$UNIREF_FASTA.gz" >> {log} 2>&1

                echo "MMseqs2: downloading NCBI taxonomy..."
                wget -q -P "$TAX_DIR" \
                    https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz >> {log} 2>&1
                tar -xzf "$TAX_DIR/taxdump.tar.gz" -C "$TAX_DIR" >> {log} 2>&1

                echo "MMseqs2: building taxonomy database (this may take several hours)..."
                mmseqs createdb "$UNIREF_FASTA" {params.db_prefix} >> {log} 2>&1
                mmseqs createtaxdb {params.db_prefix} "$TMP_DIR" \
                    --ncbi-tax-dump "$TAX_DIR" \
                    --threads {resources.threads} >> {log} 2>&1
                mmseqs createindex {params.db_prefix} "$TMP_DIR" \
                    --search-type 2 \
                    --threads {resources.threads} >> {log} 2>&1

                rm -f "$UNIREF_FASTA"
                rm -rf "$TMP_DIR" "$TAX_DIR"
                echo "MMseqs2: production database ready"
            else
                echo "MMseqs2: database already exists at {params.db_prefix}, skipping"
            fi
        fi
        touch {output.done}
        """

# ---------------------------------------------------------------------------
#   Per-sample: gene calling with Prodigal (shared output for AMR + summary)
# ---------------------------------------------------------------------------

rule prodigal:
    input:
        contigs = os.path.join(output_dir, "data", "megahit", "{sample}.contigs.fa")
    output:
        faa = os.path.join(output_dir, "data", "prodigal", "{sample}.faa"),
        gff = os.path.join(output_dir, "data", "prodigal", "{sample}.gff")
    resources:
        mem_mb  = lambda wc: res(16000, 4000),
        threads = 1,
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "{sample}.prodigal.log")
    conda: "../envs/contigs.yaml"
    shell:
        """
        mkdir -p $(dirname {output.faa})
        echo "[{wildcards.sample}] Prodigal: predicting genes..."
        prodigal -i {input.contigs} -a {output.faa} -f gff -o {output.gff} -p meta -q >> {log} 2>&1
        echo "[{wildcards.sample}] Prodigal: done"
        """

# ---------------------------------------------------------------------------
#   One-time: download AMRFinderPlus database (skipped in test mode)
# ---------------------------------------------------------------------------

rule init_amrfinder_db:
    """Download the AMRFinderPlus database to dbs/amrfinderplus_db/. Skipped in test mode."""
    output:
        done = os.path.join(_DBS_DIR, ".amrfinder_db.done")
    params:
        db_dir    = _AFP_DB_DIR,
        test_mode = _TEST
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "init_amrfinder_db.log")
    conda: "../envs/mags.yaml"
    shell:
        """
        mkdir -p {params.db_dir}
        if [ "{params.test_mode}" = "True" ]; then
            echo "init_amrfinder_db: test mode — skipping database download"
        elif [ -f "{params.db_dir}/AMR_CDS.fa" ]; then
            echo "init_amrfinder_db: database already exists, skipping download"
        else
            echo "init_amrfinder_db: downloading AMRFinderPlus database (~600 MB)..."
            amrfinder --update >> {log} 2>&1
            # amrfinder --update downloads to $CONDA_PREFIX/share/amrfinderplus/data/latest/
            AFP_SRC="$CONDA_PREFIX/share/amrfinderplus/data/latest"
            if [ -d "$AFP_SRC" ]; then
                rsync -a "$AFP_SRC/" {params.db_dir}/ >> {log} 2>&1
                echo "init_amrfinder_db: database ready"
            else
                echo "ERROR: expected database at $AFP_SRC not found" >> {log}
                exit 1
            fi
        fi
        touch {output.done}
        """

# ---------------------------------------------------------------------------
#   Per-sample: AMR annotation on contigs (protein + nucleotide + GFF mode)
# ---------------------------------------------------------------------------

rule contig_amr:
    input:
        contigs     = os.path.join(output_dir, "data", "megahit", "{sample}.contigs.fa"),
        faa         = os.path.join(output_dir, "data", "prodigal", "{sample}.faa"),
        gff         = os.path.join(output_dir, "data", "prodigal", "{sample}.gff"),
        afp_db_done = os.path.join(_DBS_DIR, ".amrfinder_db.done")
    output:
        tsv = os.path.join(output_dir, "data", "amr_contigs", "{sample}_contig_amr.tsv")
    params:
        afp_db_dir = _AFP_DB_DIR
    threads: lambda wc: res(16, 4)
    resources:
        mem_mb  = lambda wc: res(16000, 4000),
        threads = lambda wc: res(16, 4),
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "{sample}.contig_amr.log")
    conda: "../envs/mags.yaml"
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        echo "[{wildcards.sample}] AMRFinderPlus: annotating contigs..."
        amrfinder \
            -n {input.contigs} \
            -p {input.faa} \
            -g {input.gff} \
            --database {params.afp_db_dir} \
            --threads {resources.threads} \
            -o {output.tsv} \
            --annotation_format prodigal >> {log} 2>&1 || true
        # Create empty file if no hits found or database not yet set up
        [ -f {output.tsv} ] || touch {output.tsv}
        echo "[{wildcards.sample}] AMRFinderPlus: done"
        """

# ---------------------------------------------------------------------------
#   Per-sample: MGE annotation on contigs with MobileElementFinder
# ---------------------------------------------------------------------------

rule init_mge_tool:
    """Install MobileElementFinder via pip --no-deps (conda provides all deps)."""
    output:
        done = os.path.join(output_dir, "data", "mge_contigs", ".mef_install.done")
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "init_mge_tool.log")
    conda: "../envs/contigs.yaml"
    shell:
        """
        mkdir -p $(dirname {output.done})
        echo "MobileElementFinder: installing..."
        pip install --no-deps --quiet mypy-extensions >> {log} 2>&1
        pip install --no-deps --quiet MGEdb==1.1.1 >> {log} 2>&1
        pip install --no-deps --quiet MobileElementFinder >> {log} 2>&1
        touch {output.done}
        echo "MobileElementFinder: installed"
        """

rule mge_annotation:
    input:
        contigs  = os.path.join(output_dir, "data", "megahit", "{sample}.contigs.fa"),
        mef_done = os.path.join(output_dir, "data", "mge_contigs", ".mef_install.done")
    output:
        tsv = os.path.join(output_dir, "data", "mge_contigs", "{sample}_mge.tsv"),
        tmp = temp(directory(os.path.join(output_dir, "data", "mge_contigs", "{sample}_tmp")))
    threads: lambda wc: res(16, 4)
    resources:
        mem_mb  = lambda wc: res(16000, 4000),
        threads = lambda wc: res(16, 4),
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "{sample}.mge_annotation.log")
    conda: "../envs/contigs.yaml"
    shell:
        """
        mkdir -p {output.tmp}
        PREFIX=$(dirname {output.tsv})/{wildcards.sample}_mge
        echo "[{wildcards.sample}] MobileElementFinder: annotating MGEs..."
        mefinder find \
            --contig {input.contigs} \
            --temp-dir {output.tmp} \
            --threads {resources.threads} \
            "$PREFIX" >> {log} 2>&1 || true
        # mefinder writes PREFIX.csv; rename to expected .tsv
        if [ -f "${{PREFIX}}.csv" ]; then
            mv "${{PREFIX}}.csv" {output.tsv}
        else
            touch {output.tsv}
        fi
        echo "[{wildcards.sample}] MobileElementFinder: done"
        """

# ---------------------------------------------------------------------------
#   Per-sample: MMseqs2 taxonomy assignment on assembled contigs
# ---------------------------------------------------------------------------

rule mmseqs_taxonomy:
    input:
        contigs  = os.path.join(output_dir, "data", "megahit", "{sample}.contigs.fa"),
        db_done  = os.path.join(_DBS_DIR, ".mmseqs_db.done")
    output:
        lca = os.path.join(output_dir, "data", "mmseqs", "{sample}_lca.tsv")
    params:
        db_prefix = _UNIREF50_DB,
        tmp       = lambda wc: os.path.join(output_dir, "data", "mmseqs", f"{wc.sample}_tmp")
    threads: lambda wc: res(32, 4)
    resources:
        mem_mb  = lambda wc: res(64000, 4000),
        threads = lambda wc: res(32, 4),
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "{sample}.mmseqs_taxonomy.log")
    conda: "../envs/contigs.yaml"
    shell:
        """
        mkdir -p $(dirname {output.lca}) {params.tmp}
        echo "[{wildcards.sample}] MMseqs2: taxonomy classification..."
        mmseqs easy-taxonomy \
            {input.contigs} \
            {params.db_prefix} \
            $(dirname {output.lca})/{wildcards.sample} \
            {params.tmp} \
            --threads {resources.threads} \
            --lca-ranks superkingdom,phylum,class,order,family,genus,species \
            --format-output "query,taxid,taxname,taxlineage" \
            --orf-filter 0 \
            -s 2 \
            >> {log} 2>&1 || true
        # easy-taxonomy writes {{prefix}}_lca.tsv; touch if absent (no hits or search failed)
        mv $(dirname {output.lca})/{wildcards.sample}_lca.tsv {output.lca} 2>/dev/null || true
        rm -rf {params.tmp}
        [ -f {output.lca} ] || touch {output.lca}
        echo "[{wildcards.sample}] MMseqs2: done"
        """

# ---------------------------------------------------------------------------
#   Per-sample: map clean reads to all contigs; compute cpg abundance
# ---------------------------------------------------------------------------

rule contig_abundance:
    input:
        r1      = get_clean_r1,
        r2      = get_clean_r2,
        contigs = get_contigs,
        scgs    = lambda wc: (
            os.path.join(output_dir, "data", "alignments", f"{wc.sample}.scgs")
            if _RUN_SR else []
        )
    output:
        bam = os.path.join(output_dir, "data", "contig_abundance", "{sample}_contigs.bam"),
        bai = os.path.join(output_dir, "data", "contig_abundance", "{sample}_contigs.bam.bai"),
        tsv = os.path.join(output_dir, "data", "contig_abundance", "{sample}_contig_abundance.tsv")
    threads: lambda wc: res(16, 4)
    resources:
        mem_mb  = lambda wc: res(20000, 4000),
        threads = lambda wc: res(16, 4),
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "{sample}.contig_abundance.log")
    conda: "../envs/contigs.yaml"
    script:
        "../scripts/contig_abundance.py"

# ---------------------------------------------------------------------------
#   Aggregation: join all per-sample annotations into one summary TSV
# ---------------------------------------------------------------------------

rule contig_summary:
    input:
        plas_summaries = expand(os.path.join(output_dir, "data", "genomad", "{sample}", "{sample}.contigs_summary", "{sample}.contigs_plasmid_summary.tsv"), sample=samples),
        vir_summaries  = expand(os.path.join(output_dir, "data", "genomad", "{sample}", "{sample}.contigs_summary", "{sample}.contigs_virus_summary.tsv"),  sample=samples),
        lca_files      = expand(os.path.join(output_dir, "data", "mmseqs", "{sample}_lca.tsv"), sample=samples),
        abund_files    = expand(os.path.join(output_dir, "data", "contig_abundance", "{sample}_contig_abundance.tsv"), sample=samples),
        amr_files      = expand(os.path.join(output_dir, "data", "amr_contigs", "{sample}_contig_amr.tsv"), sample=samples),
        mge_files      = expand(os.path.join(output_dir, "data", "mge_contigs", "{sample}_mge.tsv"), sample=samples),
        prodigal_gffs  = expand(os.path.join(output_dir, "data", "prodigal", "{sample}.gff"), sample=samples)
    output:
        tsv = os.path.join(output_dir, "contig_summary.tsv")
    params:
        samples    = samples,
        output_dir = output_dir
    conda: "../envs/contigs.yaml"
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "contig_summary.log")
    script:
        "../scripts/contig_summary.py"
        
        
        
        
        
        