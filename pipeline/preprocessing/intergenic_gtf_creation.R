#!/usr/bin/env Rscript
#
# Build a custom GTF with intergenic regions appended following the Bgee guidelines.
#
# Run this BEFORE `integration_pipeline.R` / `snakemake`.
#
# Steps performed:
#   1. Load the annotation and split it into transcript / exon subsets,
#      restricted to the assembled chromosomes.
#   2. Derive per-gene coordinates, strand and biotype by deduplicating
#      the transcript records by gene_id.
#   3. Build candidate intergenic regions (1500-500bp upstream/downstream
#      of each gene) in parallel across chromosomes.
#   4. Trim runs of N bases from the candidate regions and drop regions
#      that end up too short.
#   5. Write a combined GTF (`gtf_all`) with the original gene/exon records
#      plus the new intergenic "exon" records.
#
# Usage:
#   Rscript intergenic_gtf_creation.R
#

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
get_annot_value <- function(split_annotation, field_name){
  field_all <- split_annotation[grep(field_name, split_annotation, fixed = T)];
  field_value <- strsplit(field_all, ' ', fixed = T)[[1]][2];
  field_value <- sub(';', '', field_value,fixed = T)
  field_value <- gsub('"', '', field_value, fixed = T)

  return(field_value)
}


remove_Ns_from_intergenic <- function (chr_number, chr_intergenic_regions, max_block_size = 31, max_proportion_N = 0.05, min_intergenic_length = 999) {

  intergenic_regions_without_N <- matrix(ncol = 6, nrow = 0)
  colnames(intergenic_regions_without_N) <- c("chr", "start", "end", "upstream/downstream","strand", "sequence")

  
  for (line in 1 : nrow(chr_intergenic_regions)) {

    chr_start <- as.numeric(chr_intergenic_regions[line,"start"])
    chr_end <- as.numeric(chr_intergenic_regions[line,"end"])
    intergenic_start <- 1
    intergenic_sequence <- chr_intergenic_regions[line,"sequence"]

    # if less N than the minimum size of a block we add the intergenic sequence to corrected intergenic regions
    if (max_block_size >=  count_number_of_occurences("N", intergenic_sequence)) {
      intergenic_regions_without_N <- rbind(intergenic_regions_without_N, chr_intergenic_regions[line,])

      # if no block of N and proportion of N lower than threshold then we add the sequence to corrected intergenic regions
    } else if (count_number_of_occurences(strrep("N", max_block_size), intergenic_sequence) ==  0
               & proportion_of_N(intergenic_sequence) <=  max_proportion_N) {
      intergenic_regions_without_N <- rbind(intergenic_regions_without_N, chr_intergenic_regions[line,])

    } else {

      seq_position <- 0
      while (seq_position + 1 <=  nchar(intergenic_sequence)) {
        seq_position <- seq_position + 1
        if (substr(intergenic_sequence, seq_position, seq_position) ==  "N") {
          block_size = 1
          while (seq_position + 1 <=  nchar(intergenic_sequence) && substr(intergenic_sequence,seq_position + 1,seq_position + 1) ==  "N") {
            block_size <- block_size + 1
            seq_position <- seq_position + 1
          }

          if (block_size >=  max_block_size) {

            current_chr_stop <- as.numeric(chr_intergenic_regions[line,"start"]) + seq_position - (block_size + 1)
            current_intergenic_stop <- seq_position - block_size
            current_size <- (current_chr_stop - chr_start) + 1
            if (current_size >=  min_intergenic_length && proportion_of_N(substr(intergenic_sequence, intergenic_start, current_intergenic_stop)) <=  max_proportion_N) {
              intergenic_regions_without_N <- rbind(intergenic_regions_without_N,
                                                    cbind(chr_number, chr_start, current_chr_stop, substr(intergenic_sequence, intergenic_start, current_intergenic_stop))[1,])
            }
            intergenic_start <- 1 + seq_position
            chr_start <- as.numeric(chr_intergenic_regions[line,"start"]) + seq_position
          }
        }
      }
      last_portion_size <- (seq_position - intergenic_start) + 1
      if (last_portion_size >=  min_intergenic_length && proportion_of_N(substr(intergenic_sequence, intergenic_start, seq_position)) <=  max_proportion_N) {
        intergenic_regions_without_N <- rbind(intergenic_regions_without_N,
                                              cbind(chr_number, chr_start, chr_end, substr(intergenic_sequence, intergenic_start, seq_position))[1,])
      }
    }

  }
  return(intergenic_regions_without_N)
}

count_number_of_occurences <- function (pattern, string) {
  return (lengths(regmatches(string, gregexpr(pattern, string))))
}

proportion_of_N <- function (sequence) {
  return (countPattern("N", sequence) / nchar(sequence))
}

# ---------------------------------------------------------------------------
# 1. Load GTF and derive gene-level metadata
# ---------------------------------------------------------------------------

gene_gtf <- as.matrix(read.table(file = "Prudul26A.chromosomes.clean.gtf", sep= "\t",strip.white = TRUE, as.is = TRUE, colClasses = "character", comment.char = '#'))

# Only the 8 assembled pseudo-chromosomes are used; unplaced scaffolds are excluded.
chromosomes_all <- c("Pd01","Pd02","Pd03","Pd04","Pd05","Pd06", "Pd07", "Pd08")

gene_gtf_exon <- gene_gtf[gene_gtf[,3] == "exon",]
gene_gtf_exon <- gene_gtf_exon[gene_gtf_exon[,1] %in% chromosomes_all,]

# The annotation only carries "transcript"/"exon" rows (no explicit "gene"
# feature), so per-gene chr/strand/biotype come from the transcript records.
gene_gtf_transcript <- gene_gtf[gene_gtf[,3] == "transcript",]
gene_gtf_transcript <- gene_gtf_transcript[gene_gtf_transcript[,1] %in% chromosomes_all,]

# tx2gene map, derived straight from the exon records' attributes.
split_annotation_list <- strsplit(gene_gtf_exon[,9], "; ", fixed = T)
gene_ids <- sapply(split_annotation_list, function(x) {get_annot_value(x, 'gene_id') })
transcript_ids <- sapply(split_annotation_list, function(x) {get_annot_value(x, 'transcript_id') })
tx2gene_ids <- unique(cbind(transcript_ids,gene_ids), MARGIN = 1)
split_annotation_list <- strsplit(gene_gtf_transcript[,9], "; ",  fixed  = T)
transcript_gene_ids <- sapply(split_annotation_list, function(x) {get_annot_value(x, 'gene_id') })
gene_gtf <- gene_gtf_transcript[!duplicated(transcript_gene_ids), ]

split_annotation_list <- strsplit(gene_gtf[,9], "; ",  fixed  = T)
gene_biotypes <- cbind(
  sapply(split_annotation_list, function(x) {get_annot_value(x, 'gene_id') }),
  sapply(split_annotation_list, function(x) {get_annot_value(x, 'gene_type') }),
  "genic")
gene_start <- sapply(split(as.numeric(gene_gtf_exon[,4]), gene_ids), function(x){ sort(as.numeric(x))[1] })
gene_stop <- sapply(split(as.numeric(gene_gtf_exon[,5]), gene_ids), function(x){ rev(sort(as.numeric(x)))[1] })
gene_chr <- sapply(split(gene_gtf_exon[,1], gene_ids), function(x){ x[1] })

chromosomes <- unique(gene_gtf_exon[,1])

summary_N_removal <- matrix(ncol = 4, nrow = 2, data = 0)
colnames(summary_N_removal) <- c("intergenic_regions", "total_bp", "total_N", "proportion_N")
rownames(summary_N_removal)<- c("before","after")

scipen_initial_value <- getOption("scipen")
options(scipen = 999)
gene_ids <- sapply(split_annotation_list, function(x) {get_annot_value(x, 'gene_id') })

# ---------------------------------------------------------------------------
# 2. Build candidate intergenic windows (1500-500bp flanking each gene)
# ---------------------------------------------------------------------------
#
# For every gene, define an upstream and a downstream window offset by
# 500-1500bp from the gene boundary (flipped for genes on the "-" strand,
# so "upstream"/"downstream" stay relative to transcription direction, not
# to chromosome coordinates).

library(parallel)

num_cores <- 10
intergenic_chr = c()
intergenic_starts = c()
intergenics_end = c()
intergenic_name = c()
intergenic_strand = c()

cl <- makeCluster(num_cores)
clusterExport(cl, varlist = c("gene_gtf", "gene_ids", "gene_chr", "gene_start", "gene_stop", "intergenic_starts", "intergenics_end", "intergenic_name","intergenic_chr","intergenic_strand"))

result <- parLapply(cl, 1:nrow(gene_gtf), function(gene_nrow) {

  intergenic_chr =  append(intergenic_chr, as.character(gene_gtf[gene_nrow, 1]))
  intergenic_chr =  append(intergenic_chr, as.character(gene_gtf[gene_nrow, 1]))
  intergenic_strand = append(intergenic_strand, as.character(gene_gtf[gene_nrow, 7]))
  intergenic_strand = append(intergenic_strand, as.character(gene_gtf[gene_nrow, 7]))

  if(gene_gtf[gene_nrow, 7] == "+"){
    intergenic_starts = append(intergenic_starts, max(0, as.numeric(gene_gtf[gene_nrow, 4]) - 1500))
    intergenics_end = append(intergenics_end, max(0, as.numeric(gene_gtf[gene_nrow, 4]) - 500))
    intergenic_name = append(intergenic_name, paste0("upstream_", gene_ids[gene_nrow]))

    intergenic_starts = append(intergenic_starts, max(0, as.numeric(gene_gtf[gene_nrow, 5]) + 500))
    intergenics_end = append(intergenics_end, max(0, as.numeric(gene_gtf[gene_nrow, 5]) + 1500))
    intergenic_name = append(intergenic_name, paste0("downstream_", gene_ids[gene_nrow]))
  } else {
    intergenic_starts = append(intergenic_starts, max(0, as.numeric(gene_gtf[gene_nrow,4]) - 1500))
    intergenics_end = append(intergenics_end, max(0, as.numeric(gene_gtf[gene_nrow, 4]) - 500))
    intergenic_name = append(intergenic_name, paste0("downstream_", gene_ids[gene_nrow]))

    intergenic_starts = append(intergenic_starts, max(0, as.numeric(gene_gtf[gene_nrow, 5]) + 500))
    intergenics_end = append(intergenics_end, max(0, as.numeric(gene_gtf[gene_nrow, 5]) + 1500))
    intergenic_name = append(intergenic_name, paste0("upstream_", gene_ids[gene_nrow]))
  }

  return(list(intergenic_chr, intergenic_starts, intergenics_end, intergenic_name,intergenic_strand))
})

stopCluster(cl)

intergenic_chr <- unlist(lapply(result, function(x) x[[1]]))
intergenic_starts <- unlist(lapply(result, function(x) x[[2]]))
intergenics_end <- unlist(lapply(result, function(x) x[[3]]))
intergenic_name <- unlist(lapply(result, function(x) x[[4]]))
intergenic_strand <- unlist(lapply(result, function(x) x[[5]]))

final_intergenic_regions <- matrix(c(intergenic_chr, intergenic_starts, intergenics_end, intergenic_name, intergenic_strand), ncol = 5)
colnames(final_intergenic_regions) <- c("chr", "start", "end", "upstream/downstream", "strand")

Reference_intergenic <- matrix(ncol = 6, nrow = 0)
colnames(Reference_intergenic) <- c("chr", "start", "end", "upstream/downstream", "strand", "sequence")

# Reference genome, used to pull the actual sequence for each candidate window.
ref_fasta_genome <- readDNAStringSet("./data/pdulcis26.chromosomes.fasta")

# ---------------------------------------------------------------------------
# 3. Restrict windows to true intergenic space, extract sequence, filter Ns
# ---------------------------------------------------------------------------
#
# Per chromosome (parallelized, since each chromosome is independent):
#   - compute the actual gene-free space on the chromosome (coverage() gap
#     regions between genes), so windows that overlap a neighbouring gene
#     get clipped out;
#   - intersect that with the 500-1500bp candidate windows from step 2;
#   - drop anything left shorter than 1001bp;
#   - pull the sequence from the reference FASTA and remove N-heavy blocks.

library(parallel)
library(IRanges)
library(Biostrings)

num_cores <- 10
cl <- makeCluster(num_cores)

# Defaults set by Bgee for N-block filtering.
N_block_size = 31
N_proportion = 0.05
clusterExport(cl, varlist = c("chromosomes", "gene_chr", "gene_start", "gene_stop",
                              "intergenic_starts", "intergenics_end", "final_intergenic_regions",
                              "ref_fasta_genome", "N_block_size", "N_proportion", "remove_Ns_from_intergenic",
                              "count_number_of_occurences", "intergenic_chr", "proportion_of_N"))

results_N <- parLapply(cl, chromosomes, function(chr) {
  library(IRanges)
  library(Biostrings)
  print(chr)
  if(( sum(gene_chr == as.character(chr)) == 0 )){ return(NULL) }
  gene_IR <- IRanges(start = gene_start[gene_chr == as.character(chr)], end = gene_stop[gene_chr == as.character(chr)])

  intergenic_all_IR <- slice(coverage(gene_IR), lower = 0, upper = 0, rangesOnly = TRUE)
  Max_chromosomal_distance = max(end(intergenic_all_IR)[length(intergenic_all_IR)], end(gene_IR)[length(gene_IR)])
  intergenic_gene1k_IR <- IRanges(start = intergenic_starts[intergenic_chr == as.character(chr)], end = intergenics_end[intergenic_chr == as.character(chr)])
  intergenic_ref_IR <- intersect(intergenic_all_IR, intergenic_gene1k_IR)
  intergenic_ref_IR <- intergenic_ref_IR[intergenic_ref_IR@width >= 1001]
  chr_intergenic_regions = final_intergenic_regions[final_intergenic_regions[,"chr"] == chr & final_intergenic_regions[, "start"] %in% intergenic_ref_IR@start & as.numeric(final_intergenic_regions[, "end"]) <= as.numeric(Max_chromosomal_distance), ]

  if(is.null(nrow(chr_intergenic_regions))){
    return(NULL)
  } else if(length(duplicated(chr_intergenic_regions[, "start"])) != 0){
    chr_intergenic_regions = chr_intergenic_regions[!duplicated(chr_intergenic_regions[, "start"]), ]
  }

  if(nrow(chr_intergenic_regions) == 0) return(NULL)

  chr_sequence <- as.character(ref_fasta_genome[[chr]])
  sequence <- apply(chr_intergenic_regions, 1, function(x) substr(chr_sequence, x["start"], x["end"]))
  chr_intergenic_regions <- cbind(chr_intergenic_regions, sequence)

  # N content before filtering, kept for the summary_N_removal report.
  total_N_before <- sum(apply(chr_intergenic_regions, 1, function(x) count_number_of_occurences("N", x["sequence"])))
  intergenic_regions_before <- nrow(chr_intergenic_regions)
  total_bp_before <- sum(as.numeric(chr_intergenic_regions[,"end"]) - as.numeric(chr_intergenic_regions[,"start"]) + 1)

  chr_intergenic_regions_after_N_removal <- remove_Ns_from_intergenic(chr, chr_intergenic_regions, as.numeric(N_block_size), as.numeric(N_proportion))

  # N content after filtering, for the same report.
  total_N_after <- sum(apply(chr_intergenic_regions_after_N_removal, 1, function(x) count_number_of_occurences("N", x["sequence"])))
  intergenic_regions_after <- nrow(chr_intergenic_regions_after_N_removal)
  total_bp_after <- sum(as.numeric(chr_intergenic_regions_after_N_removal[,"end"]) - as.numeric(chr_intergenic_regions_after_N_removal[,"start"]) + 1)

  return(list(
    regions = chr_intergenic_regions_after_N_removal[,1:6],
    before = c(total_N_before, intergenic_regions_before, total_bp_before),
    after = c(total_N_after, intergenic_regions_after, total_bp_after)
  ))
})

stopCluster(cl)

# Chromosomes with no surviving intergenic regions return NULL; drop those before rbind.
results_N_clean <- Filter(Negate(is.null), results_N)
Reference_intergenic <- do.call(rbind, lapply(results_N_clean, function(x) x$regions))

summary_N_removal["before", "total_N"] <- sum(sapply(results_N, function(x) x$before[1]))
summary_N_removal["before", "intergenic_regions"] <- sum(sapply(results_N, function(x) x$before[2]))
summary_N_removal["before", "total_bp"] <- sum(sapply(results_N, function(x) x$before[3]))

summary_N_removal["after", "total_N"] <- sum(sapply(results_N, function(x) x$after[1]))
summary_N_removal["after", "intergenic_regions"] <- sum(sapply(results_N, function(x) x$after[2]))
summary_N_removal["after", "total_bp"] <- sum(sapply(results_N, function(x) x$after[3]))

options(scipen = scipen_initial_value)

# ---------------------------------------------------------------------------
# 4. Assemble the combined GTF and write it out
# ---------------------------------------------------------------------------

cat("\nPreparing intergenic GTF data...\n")
intergenic_regions_gtf <- matrix(ncol = 9, nrow = nrow(Reference_intergenic))
intergenic_regions_gtf[,1] <- Reference_intergenic[,1]
intergenic_regions_gtf[,2] <- "intergenic"
intergenic_regions_gtf[,3] <- "exon"
intergenic_regions_gtf[,c(4,5)] <- Reference_intergenic[,c(2,3)]
intergenic_regions_gtf[,c(6,8)] <- "."
intergenic_regions_gtf[,7] <- Reference_intergenic[, 5]
intergenic_id <- Reference_intergenic[, 4]
gene_id <- paste0('gene_id "', intergenic_id, '";')
transcript_id <- paste0('transcript_id "', intergenic_id, '";')
intergenic_regions_gtf[,9] <- apply(cbind(gene_id, transcript_id), 1, function(x){ paste(x, collapse = " ") })

summary_N_removal[,"proportion_N"] <- summary_N_removal[,"total_N"] / summary_N_removal[,"total_bp"]


gtf_all <- rbind(gene_gtf_transcript, gene_gtf_exon, intergenic_regions_gtf)
write.table(x = gtf_all,
            file = paste(output_file_path, ".gtf_all", sep = ""),
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
