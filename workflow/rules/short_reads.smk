# ------------------------------------------
#         initiate databases
# ------------------------------------------

rule initiate_dbs:
    output:
        afp_db       = os.path.join(output_dir, "data", "alignments", "dbs", ".indexing.done.txt"),
        scg_db       = os.path.join(output_dir, "data", "alignments", "dbs", ".dmnd.done.txt"),
        afp_metadata = os.path.join(output_dir, "data", "alignments", "dbs", "ReferenceGeneCatalog.txt"),
        h_genome     = os.path.join(output_dir, "data", "alignments", "dbs", "GCF_000001405.40_GRCh38.p14_genomic.fna.gz")
    params:
        scg_db=config["scg_db"]
    conda: "../envs/shortreads.yaml"
    shell:
        """
        
        # download latest AMRFinderPlus database and metadata
        echo "downloading latest AMRFinderPlus database and metadata"
        wget -P $(dirname {output.afp_db}) https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMR_CDS.fa > /dev/null 2>&1
        wget -P $(dirname {output.afp_db}) https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt > /dev/null 2>&1
        
        # index AMRFinderPlus database
        echo "indexing AMRFinderPlus database"
        kma index -i $(dirname {output.afp_db})/AMR_CDS.fa -o $(dirname {output.afp_db})/afp_db > /dev/null 2>&1
        touch {output.afp_db}
        
        # make diamond database
        echo "making diamond database for single copy genes"
        diamond makedb --in {params.scg_db} -d $(dirname {output.scg_db})/scg_db > /dev/null 2>&1
        touch {output.scg_db}
        
        # download reference human genome
        echo "downloading human reference genome for host filtering"
        wget -P $(dirname {output.afp_db}) https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz > /dev/null 2>&1
        
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
        h_genome    = os.path.join(output_dir, "data", "alignments", "dbs", "GCF_000001405.40_GRCh38.p14_genomic.fna.gz")
    output:
        json        = os.path.join(output_dir, "data", "QAQC", "fastp_reports", "{sample}.json"),
        r1_clean    = protected(os.path.join(output_dir, "data", "clean_reads", "{sample}_R1.clean.fastq.gz")),
        r2_clean    = protected(os.path.join(output_dir, "data", "clean_reads", "{sample}_R2.clean.fastq.gz")),
        nonpareil   = protected(os.path.join(output_dir, "data", "QAQC", "nonpareil", "{sample}.npo")),
        afp         = os.path.join(output_dir, "data", "alignments", "{sample}.afp.res"),
        scgs        = os.path.join(output_dir, "data", "alignments", "{sample}.scgs")
    resources:
        mem_mb = 20000,
        threads = 16,
        time = "0-10:00:00"
    benchmark:
        os.path.join(output_dir, "data", "QAQC", "benchmarks", "{sample}.short_reads.txt")
    conda: "../envs/shortreads.yaml"
    shell:
        """
        set -euo pipefail
        export TMPDIR="$(dirname {output.r1_clean})"
        
        TMP_R1="$TMPDIR/{wildcards.sample}_R1.fastq.gz"
        TMP_R2="$TMPDIR/{wildcards.sample}_R2.fastq.gz"
        
        # run fastp using temporary files
        echo "cleaning reads with fastp"
        fastp -i {input.r1} -I {input.r2} -o "$TMP_R1" -O "$TMP_R2" --html /dev/null/ --json {output.json}

        # filter out human reads, ensure final output has equal number of reads
        echo "filtering out human DNA"
        
        minimap2 -t {resources.threads} -ax sr {input.h_genome} "$TMP_R1" "$TMP_R2" | samtools view -u -f 12 -F 256 - \
        | samtools fastq -n -1 >(pigz -p {resources.threads} > {output.r1_clean}) \
                 -2 >(pigz -p {resources.threads} > {output.r2_clean}) \
                 -0 /dev/null -
        
        rm -f "$TMP_R1" "$TMP_R2"
        
        # run KMA against AMRFinderPlus
        echo "aligning reads to AMRFinderPlus with KMA"
        kma -ipe {output.r1_clean} {output.r2_clean} -o $(dirname {output.afp})/{wildcards.sample}.afp -t_db $(dirname {input.afp_db})/afp_db -1t1 -t {resources.threads}
        
        # create temp single-end file for diamond
        TMP_SE="$TMPDIR/{wildcards.sample}_SE.fastq.gz"
        cat {output.r1_clean} {output.r2_clean} > "$TMP_SE"
        
        
        # Run diamond against single-copy genes for cell normalization
        echo "aligning reads to single copy genes with diamond"
        diamond blastx --db $(dirname {input.scg_db})/scg_db --query "$TMP_SE" --out {output.scgs} --outfmt 6 qseqid sseqid pident length evalue bitscore slen --fast --max-target-seqs 1 --threads {resources.threads} --quiet
        
        rm -f "$TMP_SE"
        
        # run nonpareil to calcualte metagenomic coverage
        echo "estimating metagenomic coverage with nonpareil"
        base=$(basename {output.r1_clean} .gz)
        gunzip -c {output.r1_clean} > $(dirname {output.nonpareil})/$base
        nonpareil -s $(dirname {output.nonpareil})/$base -T kmer -f fastq -b $(dirname {output.nonpareil})/{wildcards.sample}
        rm $(dirname {output.nonpareil})/$base
        
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
        afp_metadata    = os.path.join(output_dir, "data", "alignments", "dbs", "ReferenceGeneCatalog.txt")
    output:
        fastp   = os.path.join(output_dir, "fastp_summary.csv"),
        afp     = os.path.join(output_dir, "short_reads_output.csv")
    conda: "../envs/Renv.yaml"
    script:
        "../scripts/short_reads_processing.R"   
