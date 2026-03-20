# ------------------------------------------
#         initiate databases
# ------------------------------------------

rule initiate_dbs:
    output:
        afp_db       = os.path.join(output_dir, "data", "alignments", "dbs", ".indexing.done.txt"),
        scg_db       = os.path.join(output_dir, "data", "alignments", "dbs", ".dmnd.done.txt"),
        afp_metadata = os.path.join(output_dir, "data", "alignments", "dbs", "ReferenceGeneCatalog.txt"),
        h_genome     = os.path.join(output_dir, "data", "alignments", "dbs", "GCF_000001405.40_GRCh38.p14_genomic.fna.gz"),
        markers_db   = os.path.join(output_dir, "data", "alignments", "dbs", ".markers.done.txt")
    params:
        scg_db           = config["scg_db"],
        test_mode        = _TEST,
        test_amr_fa      = os.path.join(_REPO, "test", "dbs", "amrfinder", "AMR_CDS.fa"),
        test_meta        = os.path.join(_REPO, "test", "dbs", "amrfinder", "ReferenceGeneCatalog.txt"),
        test_human       = os.path.join(_REPO, "test", "dbs", "human", "human_mini.fna.gz"),
        test_markers_dir = os.path.join(_REPO, "test", "dbs", "markers"),
        markers_dir      = config.get("markers_db", ""),
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "initiate_dbs.log")
    conda: "../envs/shortreads.yaml"
    shell:
        """
        mkdir -p $(dirname {output.afp_db})

        if [ "{params.test_mode}" = "True" ]; then

            echo "initiate_dbs: indexing AMRFinderPlus database..."
            cp {params.test_amr_fa} $(dirname {output.afp_db})/AMR_CDS.fa
            kma index -i $(dirname {output.afp_db})/AMR_CDS.fa \
                      -o $(dirname {output.afp_db})/afp_db >> {log} 2>&1
            touch {output.afp_db}
            cp {params.test_meta} {output.afp_metadata}

            echo "initiate_dbs: building SCG DIAMOND database..."
            diamond makedb --in {params.scg_db} \
                           -d $(dirname {output.scg_db})/scg_db --quiet >> {log} 2>&1
            touch {output.scg_db}
            cp {params.test_human} {output.h_genome}

            echo "initiate_dbs: indexing anthropogenic markers..."
            cat {params.test_markers_dir}/pBI143.fasta {params.test_markers_dir}/crAss001.fasta \
                > $(dirname {output.afp_db})/markers.fa
            kma index -i $(dirname {output.afp_db})/markers.fa \
                      -o $(dirname {output.afp_db})/markers_db >> {log} 2>&1
            touch {output.markers_db}

        else

            echo "initiate_dbs: downloading AMRFinderPlus database..."
            wget -q -P $(dirname {output.afp_db}) https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMR_CDS.fa >> {log} 2>&1
            wget -q -P $(dirname {output.afp_db}) https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt >> {log} 2>&1

            echo "initiate_dbs: indexing AMRFinderPlus database..."
            kma index -i $(dirname {output.afp_db})/AMR_CDS.fa \
                      -o $(dirname {output.afp_db})/afp_db >> {log} 2>&1
            touch {output.afp_db}

            echo "initiate_dbs: building SCG DIAMOND database..."
            diamond makedb --in {params.scg_db} \
                           -d $(dirname {output.scg_db})/scg_db --quiet >> {log} 2>&1
            touch {output.scg_db}

            echo "initiate_dbs: downloading human reference genome..."
            wget -q -P $(dirname {output.afp_db}) https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz >> {log} 2>&1

            echo "initiate_dbs: indexing anthropogenic markers (pBI143, crAss001)..."
            cat {params.markers_dir}/pBI143.fasta {params.markers_dir}/crAss001.fasta \
                > $(dirname {output.afp_db})/markers.fa
            kma index -i $(dirname {output.afp_db})/markers.fa \
                      -o $(dirname {output.afp_db})/markers_db >> {log} 2>&1
            touch {output.markers_db}

        fi
        echo "initiate_dbs: done"
        """

# ----------------------------------------------------
#   read preprocessing, host filtering, and alignment
# ----------------------------------------------------

rule short_reads:
    input:
        r1          = get_r1,
        r2          = get_r2,
        afp_db      = os.path.join(output_dir, "data", "alignments", "dbs", ".indexing.done.txt"),
        scg_db      = os.path.join(output_dir, "data", "alignments", "dbs", ".dmnd.done.txt"),
        h_genome    = os.path.join(output_dir, "data", "alignments", "dbs", "GCF_000001405.40_GRCh38.p14_genomic.fna.gz"),
        markers_db  = os.path.join(output_dir, "data", "alignments", "dbs", ".markers.done.txt")
    output:
        json        = os.path.join(output_dir, "data", "QAQC", "fastp_reports", "{sample}.json"),
        r1_clean    = os.path.join(output_dir, "data", "clean_reads", "{sample}_R1.clean.fastq.gz"),
        r2_clean    = os.path.join(output_dir, "data", "clean_reads", "{sample}_R2.clean.fastq.gz"),
        nonpareil   = os.path.join(output_dir, "data", "QAQC", "nonpareil", "{sample}.npo"),
        afp         = os.path.join(output_dir, "data", "alignments", "{sample}.afp.res"),
        scgs        = os.path.join(output_dir, "data", "alignments", "{sample}.scgs"),
        markers     = os.path.join(output_dir, "data", "alignments", "{sample}.markers.res")
    threads: lambda wc: res(16, 4)
    resources:
        mem_mb  = lambda wc: res(20000, 4000),
        threads = lambda wc: res(16, 4),
    benchmark:
        os.path.join(output_dir, "data", "QAQC", "benchmarks", "{sample}.short_reads.txt")
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "{sample}.short_reads.log")
    conda: "../envs/shortreads.yaml"
    shell:
        """
        set -euo pipefail
        export TMPDIR="$(dirname {output.r1_clean})"
        
        TMP_R1="$TMPDIR/{wildcards.sample}_R1.fastq.gz"
        TMP_R2="$TMPDIR/{wildcards.sample}_R2.fastq.gz"
        
        echo "[{wildcards.sample}] fastp: quality filtering reads..."
        fastp -i {input.r1} -I {input.r2} -o "$TMP_R1" -O "$TMP_R2" --html /dev/null --json {output.json} >> {log} 2>&1
        echo "[{wildcards.sample}] fastp: done"

        echo "[{wildcards.sample}] minimap2: filtering host reads..."
        minimap2 -t {resources.threads} -ax sr {input.h_genome} "$TMP_R1" "$TMP_R2" 2>> {log} \
        | samtools view -u -f 12 -F 256 - 2>> {log} \
        | samtools fastq -n -1 >(pigz -p {resources.threads} > {output.r1_clean}) \
                 -2 >(pigz -p {resources.threads} > {output.r2_clean}) \
                 -0 /dev/null - 2>> {log}
        rm -f "$TMP_R1" "$TMP_R2"
        echo "[{wildcards.sample}] minimap2: done"
        
        echo "[{wildcards.sample}] KMA: aligning to AMRFinderPlus..."
        kma -ipe {output.r1_clean} {output.r2_clean} -o $(dirname {output.afp})/{wildcards.sample}.afp -t_db $(dirname {input.afp_db})/afp_db -1t1 -t {resources.threads} >> {log} 2>&1
        echo "[{wildcards.sample}] KMA: done (AMR)"
        
        echo "[{wildcards.sample}] KMA: aligning to anthropogenic markers..."
        kma -ipe {output.r1_clean} {output.r2_clean} -o $(dirname {output.markers})/{wildcards.sample}.markers -t_db $(dirname {input.markers_db})/markers_db -1t1 -t {resources.threads} >> {log} 2>&1
        echo "[{wildcards.sample}] KMA: done (markers)"
        
        echo "[{wildcards.sample}] DIAMOND: aligning to single-copy genes..."
        TMP_SE="$TMPDIR/{wildcards.sample}_SE.fastq.gz"
        cat {output.r1_clean} {output.r2_clean} > "$TMP_SE"
        diamond blastx --db $(dirname {input.scg_db})/scg_db \
            --query "$TMP_SE" \
            --out {output.scgs} \
            --outfmt 6 qseqid sseqid pident length evalue bitscore slen \
            --fast --max-target-seqs 1 --threads {resources.threads} --quiet >> {log} 2>&1
        rm -f "$TMP_SE"
        echo "[{wildcards.sample}] DIAMOND: done"
        
        echo "[{wildcards.sample}] nonpareil: estimating metagenomic coverage..."
        base=$(basename {output.r1_clean} .gz)
        gunzip -c {output.r1_clean} > $(dirname {output.nonpareil})/$base
        nonpareil -s $(dirname {output.nonpareil})/$base -T kmer -f fastq -b $(dirname {output.nonpareil})/{wildcards.sample} >> {log} 2>&1
        rm $(dirname {output.nonpareil})/$base
        echo "[{wildcards.sample}] nonpareil: done"
        
        """
        
# ------------------------------------------
#           Summarize results 
# ------------------------------------------        

rule short_reads_summary:
    input:
        afp_files       = expand(os.path.join(output_dir, "data", "alignments", "{sample}.afp.res"), sample=samples),
        scgs            = expand(os.path.join(output_dir, "data", "alignments", "{sample}.scgs"), sample=samples),
        npos            = expand(os.path.join(output_dir, "data", "QAQC", "nonpareil", "{sample}.npo"), sample=samples),
        jsons           = expand(os.path.join(output_dir, "data", "QAQC", "fastp_reports", "{sample}.json"), sample=samples),
        afp_metadata    = os.path.join(output_dir, "data", "alignments", "dbs", "ReferenceGeneCatalog.txt"),
        markers_files   = expand(os.path.join(output_dir, "data", "alignments", "{sample}.markers.res"), sample=samples)
    output:
        fastp       = os.path.join(output_dir, "fastp_summary.csv"),
        afp         = os.path.join(output_dir, "short_reads_output.csv"),
        amr_summary = os.path.join(output_dir, "markers_cpg.csv")
    conda: "../envs/Renv.yaml"
    log:
        os.path.join(output_dir, "data", "QAQC", "logs", "short_reads_summary.log")
    script:
        "../scripts/short_reads_processing.R"
