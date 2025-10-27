# ------------------------------------------
#         index PanRes
# ------------------------------------------

rule index_panres:
    output:
        pr_db = os.path.join(output_dir, "data", "panres", "database", ".indexing.done.txt")
    params:
        panres_db=config["panres_db"]
    conda: "../envs/shortreads.yaml"
    shell:
        """
        # index PanRes
        mkdir -p $(dirname {output.pr_db})
        kma index -i {params.panres_db} -o $(dirname {output.pr_db})/pr_db > $(dirname {output.pr_db})/std.out
        touch {output.pr_db}
        """

# ------------------------------------------
#           Clean, merge, and map 
# ------------------------------------------

rule short_reads:
    input:
        r1 = get_r1,
        r2 = get_r2,
        pr_db = os.path.join(output_dir, "data", "panres", "database", ".indexing.done.txt")
    output:
        clean       = protected(os.path.join(output_dir, "data", "clean_reads", "{sample}_merged.clean.fastq.gz")),
        html        = os.path.join(output_dir, "data", "QAQC", "fastp_reports", "{sample}.html"),
        mcensus     = os.path.join(output_dir, "data", "QAQC", "MicrobeCensus", "{sample}.txt"),
        nonpareil   = os.path.join(output_dir, "data", "QAQC", "nonpareil", "{sample}.nonpareil"),
        panres      = os.path.join(output_dir, "data", "panres", "{sample}.panres")
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
        mkdir -p $(dirname {output.html}) 
        echo "running fastp"
        fastp -i {input.r1} -I {input.r2} --merge --include_unmerged --merged_out {output.clean} --compression 4 --html {output.html} --json /dev/null/

        # run KMA against PanRes
        mkdir -p $(dirname {output.panres})
        echo "running KMA"
        kma -i {output.clean} -o {output.panres} -t_db $(dirname {input.pr_db})/pr_db  -1t1 -t {resources.threads}
        
        # Perform rarefaction analysis on KMA alignments
        
        # run MicrobeCensus to calculate genome equivalents (GEQ)
        export TMPDIR=$(dirname {output.mcensus})
        mkdir -p $(dirname {output.mcensus})
        echo "running MicrobeCensus"
        run_microbe_census.py {output.clean} -n 100000000 -t {resources.threads} {output.mcensus}
        
        # run nonpareil to calcualte metagenomic coverage
        mkdir -p $(dirname {output.nonpareil})
        echo "running nonpareil"
        base=$(basename {output.clean} .gz)
        gunzip -c {output.clean} > $(dirname {output.nonpareil})/$base
        nonpareil -s $(dirname {output.nonpareil})/$base -T kmer -f fastq -b {output.nonpareil}
        rm $(dirname {output.nonpareil})/$base
        
        """