#!/usr/bin/env Rscript

libraries <- c("tidyverse", "Nonpareil")

invisible(lapply(libraries, function(x) {
  suppressMessages(suppressWarnings(library(x, character.only = T)))
  }))
  
  
# Input files

afp_files       <- snakemake@input[["afp_files"]]
scgs            <- snakemake@input[["scgs"]]
npos            <- snakemake@input[["npos"]]
jsons           <- snakemake@input[["jsons"]]
afp_path        <- snakemake@input[["afp_metadata"]]

# Output files

out_file1 <- snakemake@output[[1]]
out_file2 <- snakemake@output[[2]]


###############################################
###      summarize metagenomic coverage data
###############################################

nps <- as.data.frame(summary(Nonpareil.set(npos)))
nps <- nps %>% rownames_to_column(var="sample")

#####################################################
###    generate summary table from fastp json files
#####################################################

# Robust pluck helper with default
pluck_default <- function(x, ..., .default = NA) {
  tryCatch(
    purrr::pluck(x, ..., .default = .default),
    error = function(e) .default
  )
}

# Summarize a single fastp JSON file into a one-row tibble
summarize_fastp <- function(json_path) {
  j <- jsonlite::fromJSON(json_path, simplifyVector = TRUE)
  
  # fastp's main summary is typically in j$summary$before_filtering and j$summary$after_filtering
  before <- pluck_default(j, "summary", "before_filtering", .default = list())
  after  <- pluck_default(j, "summary", "after_filtering",  .default = list())
  filtering <- pluck_default(j, "summary", "filtering_result", .default = list())
  
  # helper to safely get numeric fields (returns NA if missing)
  get_num <- function(lst, name) {
    val <- pluck_default(lst, name, .default = NA)
    if (is.null(val)) NA_real_ else as.numeric(val)
  }
  
  # Basic fields from before/after
  total_reads_before <- get_num(before, "total_reads")
  total_reads_after  <- get_num(after,  "total_reads")
  total_bases_before <- get_num(before, "total_bases")
  total_bases_after  <- get_num(after,  "total_bases")
  
  q30_before  <- get_num(before,  "q30_rate")
  q30_after   <- get_num(after,  "q30_rate")
  
  read1_len_after  <- get_num(after, "read1_mean_length")
  read2_len_after  <- get_num(after, "read2_mean_length")

  # Derived metrics: retention, percent changes, adapter trimmed percent
  pct_reads_retained <- if (!is.na(total_reads_before) && total_reads_before > 0) total_reads_after / total_reads_before * 100 else NA_real_
  pct_bases_retained <- if (!is.na(total_bases_before) && total_bases_before > 0) total_bases_after / total_bases_before * 100 else NA_real_
  avg_read_length    <-  (read1_len_after + read2_len_after) / 2
  
  tibble(
    sample = tools::file_path_sans_ext(basename(json_path)),
    
    total_reads_before = total_reads_before,
    total_reads_after  = total_reads_after,
    pct_reads_retained = pct_reads_retained,
    
    total_bases_before = total_bases_before,
    total_bases_after  = total_bases_after,
    pct_bases_retained = pct_bases_retained,
    
    q30_before = q30_before,
    q30_after  = q30_after,
    
    avg_read_length = avg_read_length
    
  )
}

fastp_summary <- purrr::map_dfr(jsons, summarize_fastp)

fastp_summary <- fastp_summary %>% left_join(nps, by="sample") %>%
                 rename(metagenomic_coverage=C, sequence_diversity=diversity) %>% 
                 select(sample, metagenomic_coverage, sequence_diversity, 
                        total_reads_before, total_reads_after, pct_reads_retained, 
                        total_bases_before, total_bases_after, pct_bases_retained,
                        q30_before, q30_after, 
                        avg_read_length)


####################################################
###      read in alignment files, normalize to cpg
####################################################

afp_metadata <- read_tsv(afp_path, show_col_types = FALSE)

catalog_long <- afp_metadata %>% 
                      pivot_longer(cols = c("refseq_nucleotide_accession", "genbank_nucleotide_accession"), values_to = "key") %>%
                      distinct(key, .keep_all = TRUE)

afp <- do.call(bind_rows, lapply(afp_files, function(f) {
  x <- tryCatch(read_tsv(f, col_types = cols()),
                error=function(e) NULL)
  if (is.null(x)) return(NULL)
  x$sample <- sub(".afp.res", "", basename(f))
  x <- x %>% relocate(sample)
  x <- x %>% separate(`#Template`, into=c("GI", "refseq_protein_accession", "refseq_nucleotide_accession",
                      "fg1", "fg2", "node_id", "parent_node_id", "res_mech_type"), sep="\\|")
  x <- x %>% separate(res_mech_type, into=c("res_mech_type", "protein_name"), sep=" ")
  x
}))


scg <- do.call(bind_rows, lapply(scgs, function(f) {
  x <- tryCatch(read.table(f),
                error=function(e) NULL)
  if (is.null(x)) return(NULL)
  colnames(x) <- c("qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "slen")
  x$sample <- sub(".scgs", "", basename(f))
  x <- x %>% relocate(sample)
  x
}))

# estimate genome counts as average coverage depth of 40 SCGs
genome.counts <- scg %>% 
                group_by(sample) %>% 
                distinct(qseqid, .keep_all = T) %>% 
                summarise(n_genomes = sum(length/slen)/40)

# estimate copies per genome (cpg) as gene coverage depth / estimated genome counts
afp_long <-  afp %>% 
                  left_join(genome.counts, by="sample") %>% 
                  left_join(fastp_summary %>% select(sample, avg_read_length), by="sample") %>%
                  mutate(read_count = round(Depth * Template_length / avg_read_length),
                         cpg = Depth / n_genomes) %>%
                  filter(read_count >= 1) %>% 
                  left_join(catalog_long, by=c("refseq_nucleotide_accession" = "key")) %>%
                  mutate(allele_pass = ifelse(Template_Identity >= 80 & Query_Identity >= 99 & !is.na(allele), "yes", "no")) %>%
                  select(sample, read_count, Depth, cpg, allele, allele_pass, gene_family, product_name, 
                         scope, type, subtype, class, subclass, refseq_nucleotide_accession) 


###############################################
###       write csv outputs
###############################################

write.csv(fastp_summary, out_file1, row.names=F)
write.csv(afp_long, out_file2, row.names=F)



