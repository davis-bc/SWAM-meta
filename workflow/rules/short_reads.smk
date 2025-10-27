# ------------------------------------------
#         initiate databases
# ------------------------------------------

rule initiate_dbs:
    output:
        pr_db   = os.path.join(output_dir, "data", "alignments", "dbs", ".indexing.done.txt"),
        scg_db  = os.path.join(output_dir, "data", "alignments", "dbs", ".dmnd.done.txt")
    params:
        panres_db=config["panres_db"],
        scg_db=config["scg_db"]
    conda: "../envs/shortreads.yaml"
    shell:
        """
        # index PanRes
        mkdir -p $(dirname {output.pr_db})
        echo "indexing PanRes database"
        kma index -i {params.panres_db} -o $(dirname {output.pr_db})/pr_db > /dev/null 2>&1
        touch {output.pr_db}
        
        # make diamond database
        echo "making diamond database"
        diamond makedb --in {params.scg_db} -d $(dirname {output.scg_db})/scg_db > /dev/null 2>&1
        touch {output.scg_db}
        
        """

# ------------------------------------------
#       QAQC + short-read gene counting 
# ------------------------------------------

rule short_reads:
    input:
        r1 = get_r1,
        r2 = get_r2,
        pr_db = os.path.join(output_dir, "data", "alignments", "dbs", ".indexing.done.txt"),
        scg_db  = os.path.join(output_dir, "data", "alignments", "dbs", ".dmnd.done.txt")
    output:
        clean       = os.path.join(output_dir, "data", "clean_reads", "{sample}_merged.clean.fastq.gz"),
        json        = os.path.join(output_dir, "data", "QAQC", "fastp_reports", "{sample}.json"),
        nonpareil   = os.path.join(output_dir, "data", "QAQC", "nonpareil", "{sample}.npo"),
        panres      = os.path.join(output_dir, "data", "alignments", "{sample}.panres.res"),
        scgs        = os.path.join(output_dir, "data", "alignments", "{sample}.scgs")
    resources:
        mem_mb = 20000,
        threads = 16,
        time = "0-10:00:00"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.short_reads.txt")
    conda: "../envs/shortreads.yaml"
    shell:
        """
        # run fastp
        mkdir -p $(dirname {output.json}) 
        echo "running fastp"
        fastp -i {input.r1} -I {input.r2} --merge --include_unmerged --merged_out {output.clean} --html /dev/null/ --json {output.json}

        # run KMA against PanRes
        mkdir -p $(dirname {output.panres})
        echo "running KMA"
        kma -i {output.clean} -o $(dirname {output.panres})/{wildcards.sample}.panres -t_db $(dirname {input.pr_db})/pr_db  -1t1 -t {resources.threads}
        
        # Run diamond against single-copy genes for cell normalization
        mkdir -p $(dirname {output.scgs})
        diamond blastx --db $(dirname {input.scg_db})/scg_db --query {output.clean} --out {output.scgs} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen --mid-sensitive --max-target-seqs 1 --threads {resources.threads} --quiet 
        
        # run nonpareil to calcualte metagenomic coverage
        mkdir -p $(dirname {output.nonpareil})
        echo "running nonpareil"
        base=$(basename {output.clean} .gz)
        gunzip -c {output.clean} > $(dirname {output.nonpareil})/$base
        nonpareil -s $(dirname {output.nonpareil})/$base -T kmer -f fastq -b $(dirname {output.nonpareil})/{wildcards.sample}
        rm $(dirname {output.nonpareil})/$base
        """
        
# ------------------------------------------
#           Summarize results 
# ------------------------------------------        

rule short_reads_summary:
    input:
        panres_files  = expand(os.path.join(output_dir, "data", "alignments", "{sample}.panres.res"), sample=samples),
        scgs          = expand(os.path.join(output_dir, "data", "alignments", "{sample}.scgs"), sample=samples),
        npos          = expand(os.path.join(output_dir, "data", "QAQC", "nonpareil", "{sample}.npo"), sample=samples),
        jsons         = expand(os.path.join(output_dir, "data", "QAQC", "fastp_reports", "{sample}.json"), sample=samples)
    output:
        fastp   = os.path.join(output_dir, "fastp_summary.csv"),
        panres  = os.path.join(output_dir, "panres.csv")
    conda: "../envs/Renv.yaml"
    script:
        "../scripts/short_reads_processing.R"   
