# ------------------------------------------
#       de novo assembly + gene annotation 
# ------------------------------------------

rule contigs:
    input:
        reads = os.path.join(output_dir, "data", "clean_reads", "{sample}_merged.clean.fastq.gz")
    output:
        megahit     = os.path.join(output_dir, "data", "megahit", "{sample}", "final.contigs.fa")
    resources:
        mem_mb = 150000,
        threads = 32,
        time = "3-00:00:00"
    benchmark:
        os.path.join(output_dir, "data", "benchmarks", "{sample}.contigs.txt")
    conda: "../envs/contigs.yaml"
    shell:
        """
        # run megahit
        mkdir -p $(dirname {output.megahit}) 
        echo "running megahit"
        megahit --12 {input.reads} -t 32 --presets meta-large --min-contig-len 1000 -o $(dirname {output.megahit})

        # map and calculate coverage depth
        minimiap
        
        # classify contigs into chromosome/phage/plasmid
        genomad 
        
        # annotate ARGs using KMA against PanRes
        kma
        
        # annotate MGEs using MobileElementFinder
        
        
        # annotate taxonomy using mmseqs2 against GTDB
        mmseqs
        
        
        """
    