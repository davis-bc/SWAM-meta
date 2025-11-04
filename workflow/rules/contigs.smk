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
        
        echo "Preparing geNomad database"
        cd $(dirname {output.genomad_db})
        genomad download-database . > /dev/null 2>&1
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
        reads       = os.path.join(output_dir, "data", "clean_reads", "{sample}_merged.clean.fastq.gz")
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
        megahit --12 {input.reads} \
            -t {resources.threads} \
            --presets meta-large \
            --min-contig-len 1000 \
            -o "$tmpdir" \
            --out-prefix {wildcards.sample} \
            --continue

        # Move the contigs file into the final megahit directory.
        # Overwrite if exists (protected() in Snakemake prevents accidental deletion by Snakemake, but we still overwrite intentionally).
        mv -f "$tmpdir/{wildcards.sample}.contigs.fa" {output.megahit}

        # Clean up the rest of temporary files
        rm -rf "$tmpdir"

        """
    
# ------------------------------------------------------------
#       classify contigs into chromosome, plasmid, and phage 
# ------------------------------------------------------------

rule genomad:
    input:
        contigs     = os.path.join(output_dir, "data", "megahit", "{sample}.contigs.fa"),
        genomd_db   = os.path.join(output_dir, "data", "genomad", ".genomad.db.done.txt")
    output:
        genomad     = directory(os.path.join(output_dir, "data", "genomad", "{sample}"))
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
        genomad end-to-end {output.megahit} {output.genomad} $(dirname {input.genomad_db})/genomad_db --relaxed --cleanup
        
        """
        
        
        
        
        