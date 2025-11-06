# ------------------------------------------
#       prepare geNomad database 
# ------------------------------------------

rule init_genomad:
    output:
        genomad_db   = os.path.join(output_dir, "data", "genomad", ".genomad.db.done.txt") 
    conda: "../envs/genomad.yaml"
    shell:
        """
        if [ ! -d $(dirname {output.genomad_db})/genomad_db ]; then
        
        cd $(dirname {output.genomad_db})
        genomad download-database .
        touch {output.genomad_db}
        
        else 
        
        echo "geNomad database exists, skipping setup"
        touch {output.genomad_db}
        
        fi
        """
# ------------------------------------------
#       de novo assembly + gene annotation 
# ------------------------------------------

rule contigs:
    input:
        r1       = os.path.join(output_dir, "data", "clean_reads", "{sample}_R1.clean.fastq.gz"),
        r2       = os.path.join(output_dir, "data", "clean_reads", "{sample}_R2.clean.fastq.gz")
    output:
        megahit     = protected(os.path.join(output_dir, "data", "megahit", "{sample}.contigs.fa"))
    resources:
        mem_mb = 150000,
        threads = 32,
        time = "3-00:00:00"
    benchmark:
        os.path.join(output_dir, "data", "QAQC", "benchmarks", "{sample}.contigs.txt")
    conda: "../envs/contigs.yaml"
    shell:
        """
        set -euo pipefail
        echo "assembling with megahit"

        # create temporary directory in same parent so moves are on same fs (atomic)
        tmpdir="$(mktemp -d "$(dirname {output.megahit})/.megahit_{wildcards.sample}.tmp.XXXXXX")"

        # Run MEGAHIT into the temporary directory
        megahit -1 {input.r1} -2 {input.r2} \
                -t {resources.threads} \
                --presets meta-large \
                --min-contig-len 1000 \
                -o "$tmpdir" \
                --out-prefix {wildcards.sample} \
                --continue \
                -f 

        # Move the contigs file into the final megahit directory.
        # Overwrite if exists (protected() in Snakemake prevents accidental deletion by Snakemake, but we still overwrite intentionally).
        mv -f "$tmpdir/{wildcards.sample}.contigs.fa" {output.megahit}

        # Clean up the rest of temporary files
        rm -rf "$tmpdir"
        
        # rename all contigs to include sample information
        awk -v sample="{wildcards.sample}" ' /^>/ { header = substr($0, 2) # drop leading ">" match(header, /^[^[:space:]]+/) # find first token (contig id) id = substr(header, RSTART, RLENGTH) rest = substr(header, RLENGTH + 1) # the rest (may be "") printf(">%s_%s%s\n", sample, id, rest) next } { print } ' {output.megahit} > {output.megahit}.tmp && mv -f {output.megahit}.tmp {output.megahit}

        """
    
# ------------------------------------------------------------------------
#       classify contigs into chromosome, plasmid, and phage with geNomad
# ------------------------------------------------------------------------

rule genomad:
    input:
        r1           = os.path.join(output_dir, "data", "clean_reads", "{sample}_R1.clean.fastq.gz"),
        r2           = os.path.join(output_dir, "data", "clean_reads", "{sample}_R2.clean.fastq.gz"),
        contigs      = os.path.join(output_dir, "data", "megahit", "{sample}.contigs.fa"),
        genomad_db   = os.path.join(output_dir, "data", "genomad", ".genomad.db.done.txt")
    output:
        genomad      = directory(os.path.join(output_dir, "data", "genomad", "{sample}")),
        plas_contigs = os.path.join(output_dir, "data", "genomad", "{sample}", "{sample}_plasmid.fna"),
        contig_bam   = os.path.join(output_dir, "data", "genomad", "{sample}", "{sample}.plasmid.bam"),
        circular_txt = os.path.join(output_dir, "data", "genomad", "{sample}", "{sample}.contigs-circular.txt")
    resources:
        mem_mb = 150000,
        threads = 32,
        time = "3-00:00:00"
    benchmark:
        os.path.join(output_dir, "data", "QAQC", "benchmarks", "{sample}.genomad.txt")
    conda: "../envs/genomad.yaml"
    shell:
        """
        # classify contigs into chromosome/phage/plasmid
        echo "running geNomad"
        genomad end-to-end {input.contigs} {output.genomad} $(dirname {input.genomad_db})/genomad_db --relaxed --cleanup
        
        # map reads to plasmid contigs
        minimap2 -t {resources.threads} -ax sr {output.plas_contigs} {input.r1} {input.r2} | samtools sort -o {output.contig_bam}
        
        # detect circularity of plasmid contigs
        workflow/scripts/detect_circular_contigs.py -b {output.contig_bam} -o {output.circular_txt}
        
        """

# ----------------------------------------------------------------------------------
#       determine cirularity of plasmid contigs, infer plasmid systems with MobMess 
# ----------------------------------------------------------------------------------

rule mobmess:
    input:
        circular_txt    = os.path.join(output_dir, "data", "genomad", "{sample}", "{sample}.contigs-circular.txt"),
        plas_contigs    = os.path.join(output_dir, "data", "genomad", "{sample}", "{sample}_plasmid.fna")
    output:
        mobmess         = os.path.join(output_dir, "data", "mobmess", "{sample}")
    conda: "../envs/contigs.yaml"
    resources:
        mem_mb = 150000,
        threads = 32,
        time = "1-00:00:00"
    shell:
        """
        # run mobmess to infer plasmid systems 
        mobmess systems \
        --sequences {input.plas_contigs} \
        --complete {input.circular_txt} \
        --output {output.mobmess}-mobmess \
        --threads {resources.threads}
        """
        
        
        
        
        
        