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
  
  # merge information can be stored in different keys depending on fastp version/flags
  merge_info <- pluck_default(j, "merge_result", .default = NULL)
  if (is.null(merge_info) || identical(merge_info, NA)) {
    merge_info <- pluck_default(j, "merged_result", .default = NULL)
  }
  if (is.null(merge_info) || identical(merge_info, NA)) {
    merge_info <- pluck_default(j, "merge", .default = list())
  }
  if (is.null(merge_info)) merge_info <- list()
  
  # Some fastp versions (or when merge+filter is enabled) use "merged_and_filtered"
  merged_and_filtered <- pluck_default(j, "merged_and_filtered", .default = NULL)
  if (is.null(merged_and_filtered) || identical(merged_and_filtered, NA)) {
    merged_and_filtered <- NULL
  }
  
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
  
  q20_before <- get_num(before, "q20_rate")
  q30_before <- get_num(before, "q30_rate")
  q20_after  <- get_num(after,  "q20_rate")
  q30_after  <- get_num(after,  "q30_rate")
  
  gc_before <- get_num(before, "gc_content")
  gc_after  <- get_num(after,  "gc_content")
  
  read1_len_before <- get_num(before, "read1_mean_length")
  read2_len_before <- get_num(before, "read2_mean_length")
 
  merged_reads <- get_num(merge_info, "merged_reads")
  merged_bases <- get_num(merge_info, "merged_bases")
  avg_merged_len <- get_num(merge_info, "avg_merged_len")
  merged_rate <- get_num(merge_info, "merged_rate")  # often fraction between 0 and 1
  
  # If merged_and_filtered block exists, prefer its fields (it contains merged+filtered merged-read stats)
  if (!is.null(merged_and_filtered)) {
    maf_reads  <- get_num(merged_and_filtered, "total_reads")
    maf_bases  <- get_num(merged_and_filtered, "total_bases")
    maf_q20    <- get_num(merged_and_filtered, "q20_bases")
    maf_q30    <- get_num(merged_and_filtered, "q30_bases")
    maf_cycles <- get_num(merged_and_filtered, "total_cycles")
    
    # prefer these values if they look valid
    if (!is.na(maf_reads) && maf_reads > 0) merged_reads <- maf_reads
    if (!is.na(maf_bases)) merged_bases <- maf_bases
    # compute average merged length if possible
    if (is.na(avg_merged_len) && !is.na(merged_reads) && merged_reads > 0 && !is.na(merged_bases)) {
      avg_merged_len <- merged_bases / merged_reads
    } else if (!is.na(maf_bases) && !is.na(maf_reads) && maf_reads > 0) {
      avg_merged_len <- maf_bases / maf_reads
    }
    
    # compute Q20/Q30 rates on merged reads (percent)
    merged_q20_rate <- if (!is.na(maf_q20) && !is.na(merged_bases) && merged_bases > 0) (maf_q20 / merged_bases) * 100 else NA_real_
    merged_q30_rate <- if (!is.na(maf_q30) && !is.na(merged_bases) && merged_bases > 0) (maf_q30 / merged_bases) * 100 else NA_real_
    
    # If merged_rate not provided, we can compute pct_reads_merged from merged_reads / total_reads_before
    if (!is.na(merged_reads) && !is.na(total_reads_before) && total_reads_before > 0) {
      pct_reads_merged <- merged_reads / total_reads_before * 100
    } else {
      pct_reads_merged <- NA_real_
    }
    
  } else {
    # No merged_and_filtered; fallback to older/other structures
    if (is.na(avg_merged_len) && !is.na(merged_reads) && merged_reads > 0 && !is.na(merged_bases)) {
      avg_merged_len <- merged_bases / merged_reads
    }
    merged_q20_rate <- NA_real_
    merged_q30_rate <- NA_real_
    pct_reads_merged <- NA_real_
    if (!is.na(merged_reads) && !is.na(total_reads_before) && total_reads_before > 0) {
      pct_reads_merged <- merged_reads / total_reads_before * 100
    } else if (!is.na(merged_rate)) {
      if (merged_rate <= 1) pct_reads_merged <- merged_rate * 100 else pct_reads_merged <- merged_rate
    }
  }
  
  # Derived metrics: retention, percent changes, adapter trimmed percent
  pct_reads_retained <- if (!is.na(total_reads_before) && total_reads_before > 0) total_reads_after / total_reads_before * 100 else NA_real_
  pct_bases_retained <- if (!is.na(total_bases_before) && total_bases_before > 0) total_bases_after / total_bases_before * 100 else NA_real_
  
  
  tibble(
    sample = tools::file_path_sans_ext(basename(json_path)),
    
    total_reads_before = total_reads_before,
    total_reads_after  = total_reads_after,
    pct_reads_retained = pct_reads_retained,
    
    total_bases_before = total_bases_before,
    total_bases_after  = total_bases_after,
    pct_bases_retained = pct_bases_retained,
    
    q20_before = q20_before,
    q30_before = q30_before,
    q20_after  = q20_after,
    q30_after  = q30_after,
    
    gc_before = gc_before,
    gc_after  = gc_after,
    
    read1_len_before = read1_len_before,
    read2_len_before = read2_len_before,
    
    merged_reads = merged_reads,
    merged_bases = merged_bases,
    avg_merged_len = maf_cycles,
    pct_reads_merged = pct_reads_merged
  )
}

fastp_summary <- purrr::map_dfr(jsons, summarize_fastp)
fastp_summary <- fastp_summary %>% 
                 mutate(avg_read_length = (avg_merged_len * pct_reads_merged / 100) + 
                                          (read1_len_before * (1 - (pct_reads_merged / 100)))
                                          )

fastp_summary <- fastp_summary %>% left_join(nps, by="sample") %>%
                 rename(meta_coverage=C, seq_div=diversity) %>% 
                 select(sample, meta_coverage, seq_div, 
                        total_reads_before, total_reads_after, 
                        total_bases_before, total_bases_after,
                        q30_before, q30_after, 
                        read1_len_before, read2_len_before,
                        merged_reads, merged_bases, avg_merged_len, 
                        pct_reads_merged, avg_read_length)


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



