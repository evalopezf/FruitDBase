
#!/usr/bin/env Rscript
#
# Integrate per-sample Kallisto quantifications into a gene-level expression
# atlas and expression score matrices.
#
# Run this AFTER `snakemake` has completed:
#
# Steps performed:
#   1. Collapse technical replicates (same biological sample sequenced twice),
#      keeping the run with the most processed reads.
#   2. Split the transcriptome FASTA into coding and intergenic sequences .
#   3. Aggregate transcript-level Kallisto estimates to gene level (tximport).
#   4. Call presence/absence per sample with BgeeCall.
#   5. Build an expression score matrix (by tissue, and by tissue_developmental
#      stage)
#
# Usage:
#   Rscript integration_pipeline.R \
#     --metadata metadata.tsv \
#     --kallisto-report results/kallisto/kallisto_report.rds \
#     --kallisto-dir results/kallisto \
#     --fasta transcriptome.fasta \
#     --gtf annotation.gtf \
#     --species-id 3755 \
#     --output-dir results/integration
#


.libPaths(c("~/R-libraries", .libPaths()))

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(data.table)
  library(fs)
  library(tximport)
  library(BgeeCall)
  library(Biostrings)
  library(GenomicFeatures)
  library(txdbmaker)
  library(AnnotationDbi)
  library(rtracklayer)
})
# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

parse_cli_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  option_list <- list(
    make_option("--metadata", type = "character",
                help = "Curated sample metadata (.tsv). Requires columns: external_id_sample, IDBiologicalSampleFruitDBase, tissue_name, development_stage."),
    make_option("--kallisto-report", type = "character", dest = "kallisto_report",
                help = "kallisto_report.rds produced by the 'kallisto_report' Snakemake rule."),
    make_option("--kallisto-dir", type = "character", dest = "kallisto_dir",
                help = "Kallisto output directory from Snakemake (one subfolder per sample, each with abundance.tsv)."),
    make_option("--fasta", type = "character",
                help = "Transcriptome FASTA used to build the Kallisto index (coding + intergenic sequences; intergenic headers prefixed 'upstream_'/'downstream_')."),
    make_option("--gtf", type = "character",
                help = "GTF annotation matching --fasta. Used to derive the transcript-to-gene map and for BgeeCall presence/absence calling."),
    make_option("--species-id", type = "character", dest = "species_id",
                help = "NCBI Taxonomy ID for the species (passed to BgeeCall)."),
    make_option("--output-dir", type = "character", dest = "output_dir", default = "results",
                help = "Directory where result .rds files are written [default: %default]."),
    make_option("--work-dir", type = "character", dest = "work_dir", default = "tmp_bgeecall",
                help = "Scratch directory for BgeeCall's intermediate files [default: %default].")
  )

  parser <- OptionParser(
    option_list = option_list,
    description = "Integrate per-sample Kallisto quantifications into a gene-level expression atlas."
  )
  opt <- optparse::parse_args(parser, args = argv)

  required <- c("metadata", "kallisto_report", "kallisto_dir", "fasta", "gtf", "species_id")
  missing <- required[vapply(required, function(k) is.null(opt[[k]]), logical(1))]
  if (length(missing) > 0) {
    print_help(parser)
    stop("Missing required argument(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  opt
}

# ---------------------------------------------------------------------------
# 1. Technical replicate aggregation
# ---------------------------------------------------------------------------

# When the same biological sample was sequenced more than once (technical
# replicates), keep only the run with the most processed reads.
# We use IDBiologicalSampleFruitDBase but BioSample is the parameter provided
# by SRA to identify technical replicates. 
aggregate_technical_replicates <- function(metadata, kallisto_report) {
  kallisto_report$biosample <- metadata$IDBiologicalSampleFruitDBase[
    match(kallisto_report$sample, metadata$external_id_sample)
  ]

  duplicated_biosamples <- metadata %>%
    group_by(IDBiologicalSampleFruitDBase) %>%
    filter(n() > 1) %>%
    ungroup()

  duplicated_runs <- kallisto_report %>%
    filter(sample %in% duplicated_biosamples$external_id_sample)

  best_run_per_biosample <- duplicated_runs %>%
    group_by(biosample) %>%
    slice_max(n_processed, n = 1) %>%
    ungroup()

  non_duplicated <- metadata %>%
    filter(!IDBiologicalSampleFruitDBase %in% best_run_per_biosample$biosample) %>%
    pull(external_id_sample)

  c(non_duplicated, best_run_per_biosample$sample)
}

# ---------------------------------------------------------------------------
# 2. FASTA splitting (coding vs. intergenic)
# ---------------------------------------------------------------------------

split_transcriptome_fasta <- function(fasta_file, work_dir) {
  fasta_all <- readDNAStringSet(fasta_file)
  intergenic_idx <- grepl("^upstream_|^downstream_", names(fasta_all))

  intergenic_path <- file.path(work_dir, "transcripts_intergenic.fasta")
  genes_path <- file.path(work_dir, "transcripts_genes.fasta")

  writeXStringSet(fasta_all[intergenic_idx], intergenic_path)
  writeXStringSet(fasta_all[!intergenic_idx], genes_path)

  list(intergenic = intergenic_path, genes = genes_path)
}

# ---------------------------------------------------------------------------
# 3. Load per-sample Kallisto results and aggregate to gene level
# ---------------------------------------------------------------------------

# Read each sample's abundance.tsv straight from the Snakemake kallisto
# output directory.
load_kallisto_results <- function(kallisto_dir, samples) {
  result_list <- lapply(samples, function(s) {
    path <- file.path(kallisto_dir, s, "abundance.tsv")
    if (!file.exists(path)) {
      stop("Missing Kallisto output for sample '", s, "': ", path, call. = FALSE)
    }
    df <- as.data.frame(data.table::fread(path))
    rownames(df) <- df$target_id
    df
  })
  names(result_list) <- samples
  result_list
}

# Derive a transcript_id to gene_id map directly from the GTF's TxDb
derive_tx2gene <- function(txdb) {
  tx_ids <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
  AnnotationDbi::select(txdb, keys = tx_ids, keytype = "TXNAME", columns = "GENEID") %>%
    transmute(transcript_id = as.character(TXNAME), gene_id = as.character(GENEID)) %>%
    distinct()
}

# Aggregate transcript-level estimates to gene level with tximport.
aggregate_to_gene_level <- function(result_list, tx2gene) {
  tx_ids <- result_list[[1]]$target_id
  coding_ids <- tx_ids[!grepl("^(upstream_|downstream_)", tx_ids)]

  est_counts_mat <- matrix(NA, nrow = length(tx_ids), ncol = length(result_list),
                            dimnames = list(tx_ids, names(result_list)))
  tpm_mat <- est_counts_mat
  efflen_mat <- est_counts_mat

  for (sample_name in names(result_list)) {
    df <- result_list[[sample_name]]
    est_counts_mat[, sample_name] <- df$est_counts
    tpm_mat[, sample_name] <- df$tpm
    efflen_mat[, sample_name] <- df$eff_length
  }

  txi <- list(
    abundance = tpm_mat[coding_ids, , drop = FALSE],
    counts = est_counts_mat[coding_ids, , drop = FALSE],
    length = efflen_mat[coding_ids, , drop = FALSE]
  )

  summarizeToGene(txi, tx2gene, countsFromAbundance = "lengthScaledTPM")
}

# ---------------------------------------------------------------------------
# 4. Presence/absence calling with BgeeCall
# ---------------------------------------------------------------------------

build_txdb <- function(gtf_file, species_id) {
  txdbmaker::makeTxDbFromGFF(file = gtf_file, format = "gtf", taxonomyId = as.integer(species_id))
}


call_presence_absence <- function(samples, kallisto_results, fasta_paths, gtf_file,
                                  species_id, work_dir) {

  gtf_obj <- rtracklayer::import(gtf_file, format = "gtf",
                                  colnames = union(txdbmaker:::GTF_COLNAMES,
                                                    c("source", "gene_type", "gene_biotype")),
                                  feature.type = txdbmaker:::GFF_FEATURE_TYPES)
  bgee <- new("BgeeMetadata", intergenic_release = "custom")
  user_meta_base <- list(
    custom_intergenic_path = fasta_paths$intergenic,
    species_id             = species_id,
    transcriptome_name     = fasta_paths$genes,
    annotation_object      = gtf_obj,
    gtf_source             = "gencode",
    working_path           = work_dir,
    simple_arborescence    = TRUE,
    verbose                = TRUE
  )


  stale_cache <- list.files(work_dir, pattern = "tx2gene|gene2biotype|intergenic_ids",
                            recursive = TRUE, full.names = TRUE)
  if (length(stale_cache) > 0) file.remove(stale_cache)

  kallisto_meta <- new("KallistoMetadata")
  all_calls <- list()

  for (s in samples) {
    message("Calling presence/absence for sample: ", s)

    user_meta <- do.call(new, c("UserMetadata", user_meta_base,
                                list(run_ids = s, rnaseq_lib_path = s)))

    abundance_dest <- BgeeCall:::get_abundance_file_path(kallisto_meta, bgee, user_meta)
    fs::dir_create(dirname(abundance_dest))
    data.table::fwrite(kallisto_results[[s]], abundance_dest, sep = "\t")

    all_calls[[s]] <- tryCatch(
      BgeeCall:::generate_presence_absence(
        myAbundanceMetadata = kallisto_meta,
        myUserMetadata      = user_meta,
        myBgeeMetadata      = bgee
      ),
      error = function(e) {
        message("ERROR calling presence/absence for ", s, ": ", e$message)
        NULL
      }
    )
  }

  all_results <- list()
  for (s in samples) {
    calls_dir <- file.path(work_dir, "intergenic_custom", "all_results", paste0(s, "_", s))
    calls_file <- list.files(calls_dir, full.names = TRUE)

    if (length(calls_file) > 0) {
      df <- fread(calls_file[3])
      df$sample <- s
      all_results[[s]] <- df
    } else {
      message("No BgeeCall result found for sample: ", s)
    }
  }

  rbindlist(all_results, fill = TRUE)
}

# ---------------------------------------------------------------------------
# 5. Expression atlas and expression score
# ---------------------------------------------------------------------------

# Summarise presence/absence calls into an expression atlas for an arbitrary
# grouping (e.g. tissue, or tissue x developmental stage).
build_expression_atlas <- function(calls_df, metadata, group_vars) {
  calls_df %>%
    left_join(metadata, by = "sample") %>%
    group_by(across(all_of(c(group_vars, "id")))) %>%
    summarise(
      n_samples      = n(),
      n_present      = sum(call == "present"),
      prop_present   = mean(call == "present"),
      expressed_gold   = any(p.adjust(pValue, method = "BH") <= 0.01),
      expressed_silver = any(p.adjust(pValue, method = "BH") <= 0.05),
      expressed        = expressed_silver,
      .groups = "drop"
    )
}

# Fractional gene rank per sample (highest TPM = rank 1), used as input to
# the weighted expression score below.
gene_rank_matrix <- function(tpm_matrix) {
  apply(tpm_matrix, 2, function(x) rank(-x, ties.method = "average"))
}

# Weighted mean of ranks per gene and condition, weighted by the number of
# distinct ranks contributed by each sample.
weighted_mean_rank <- function(gene_ranks, condition, weights) {
  tapply(seq_along(condition), condition, function(idx) {
    sum(weights[idx] * gene_ranks[idx]) / sum(weights[idx])
  })
}

# Expression score, scaled from the weighted mean rank as BgeeCall's ExpressionScore does.
compute_expression_score <- function(tpm_rank_matrix, condition_labels, output_path) {
  distinct_counts <- apply(tpm_rank_matrix, 2, function(x) length(unique(x)))

  weighted_ranks <- t(apply(tpm_rank_matrix, 1, function(gene_ranks) {
    weighted_mean_rank(gene_ranks, condition_labels, distinct_counts)
  }))
  weighted_ranks <- as.data.frame(weighted_ranks)
  weighted_ranks <- cbind(Gene = rownames(weighted_ranks), weighted_ranks)

  max_rank_per_cond <- apply(weighted_ranks[, -1, drop = FALSE], 2, max)
  max_rank_global <- max(max_rank_per_cond)

  score <- weighted_ranks
  score[, -1] <- ((max_rank_global + 1 - weighted_ranks[, -1]) * 100) / max_rank_global
  colnames(score)[-1] <- paste0("Score_", colnames(score)[-1])
  rownames(score) <- score$Gene

  score_matrix <- as.matrix(score[, -1, drop = FALSE])
  saveRDS(score_matrix, output_path)
  score_matrix
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  opt <- parse_cli_args()

  fs::dir_create(opt$output_dir)
  fs::dir_create(opt$work_dir)

  message("Reading metadata and kallisto report ...")
  metadata <- read_delim(opt$metadata,delim='\t')
  metadata$sample <- metadata$external_id_sample
  kallisto_report <- readRDS(opt$kallisto_report)

  message("Collapsing technical replicates ...")
  samples <- aggregate_technical_replicates(metadata, kallisto_report)

  message("Splitting transcriptome FASTA into coding / intergenic ...")
  fasta_paths <- split_transcriptome_fasta(opt$fasta, opt$work_dir)

  message("Loading per-sample Kallisto results from ", opt$kallisto_dir, " ...")
  kallisto_results <- load_kallisto_results(opt$kallisto_dir, samples)

  message("Building TxDb from GTF and deriving transcript-to-gene map ...")
  txdb <- build_txdb(opt$gtf, opt$species_id)
  tx2gene <- derive_tx2gene(txdb)

  message("Aggregating to gene level (tximport) ...")
  txi_gene <- aggregate_to_gene_level(kallisto_results, tx2gene)
  tpm_gene <- txi_gene$abundance

  message("Calling presence/absence with BgeeCall ...")

  all_calls_df <- call_presence_absence(
    samples, kallisto_results, fasta_paths, opt$gtf, opt$species_id, opt$work_dir
  )
  saveRDS(all_calls_df, file.path(opt$output_dir, "presence_absence_df.rds"))
  message("Samples processed: ", length(unique(all_calls_df$sample)))

  message("Building expression atlas ...")
  atlas_by_tissue <- build_expression_atlas(all_calls_df, metadata, "tissue_name")
  saveRDS(atlas_by_tissue, file.path(opt$output_dir, "presence_by_tissue.rds"))

  atlas_by_dev <- build_expression_atlas(all_calls_df, metadata, c("tissue_name", "development_stage"))
  saveRDS(atlas_by_dev, file.path(opt$output_dir, "presence_by_dev.rds"))

  message("Computing expression score by tissue ...")
  genes_expressed <- atlas_by_tissue %>% filter(expressed) %>% pull(id) %>% unique()
  tpm_expressed <- tpm_gene[rownames(tpm_gene) %in% genes_expressed, , drop = FALSE]
  rank_matrix <- gene_rank_matrix(tpm_expressed)
  tissue_labels <- metadata$tissue_name[match(colnames(rank_matrix), metadata$external_id_sample)]
  compute_expression_score(rank_matrix, tissue_labels,
                            file.path(opt$output_dir, "ExpressionScorebyTissue.rds"))

  message("Computing expression score by tissue x developmental stage ...")
  genes_expressed <- atlas_by_dev %>% filter(expressed) %>% pull(id) %>% unique()
  tpm_expressed <- tpm_gene[rownames(tpm_gene) %in% genes_expressed, , drop = FALSE]
  rank_matrix <- gene_rank_matrix(tpm_expressed)
  rank_matrix <- rank_matrix[, colnames(rank_matrix) %in% metadata$external_id_sample, drop = FALSE]
  dev_labels <- paste0(
    metadata$tissue_name[match(colnames(rank_matrix), metadata$external_id_sample)], "_",
    metadata$development_stage[match(colnames(rank_matrix), metadata$external_id_sample)]
  )
  compute_expression_score(rank_matrix, dev_labels,
                            file.path(opt$output_dir, "ExpressionScoreByDev.rds"))

  message("Integration complete. Results written to ", opt$output_dir)
}

if (sys.nframe() == 0) {
  main()
}
