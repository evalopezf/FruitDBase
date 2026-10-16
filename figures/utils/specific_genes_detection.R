# =============================================================================
# Utils
# Identify tissue-specific genes and developmental stage-specific genes
# Calculate medain tpm by tissue and developmetnal stage
# Input:  metadata(tsv), expression score matrix (obtain by integration_pipeline_v2.R),
#         presence absence calling (obtain by integration_pipeline_v2.R), tpm matrix (obtain by integration_pipeline_v2.R)
# Output: highlySpecificGeneBgee.rds, highlySpecificGeneMeanTPM.rds, median_tpm_by_tissue.rds, median_tpm_by_tissue_dev.rds
# =============================================================================


# ---- Libraries --------------------------------------------------------------
library(readr)                 
library(dplyr)                 
library(tidyr)
library(ggplot2)               
library(patchwork)


# ---- Function -------------------------------------------------------------

classify_with_median <- function(expr_matrix, tau_scores) {

  results <- data.frame(
    gene = rownames(expr_matrix),
    tau = tau_scores,
    category = NA,
    n_above_mean = NA,
    tissues_above_mean = NA,
    max_tissue = NA,
    n_expressed = NA
  )
  
  for(i in 1:nrow(expr_matrix)) {
    expr <- as.numeric(expr_matrix[i, ])
    names(expr) <- colnames(expr_matrix)
    gene_id <- rownames(expr_matrix)[i] 
    
    tau <- tau_scores[i]
    
    # Clasificación
     if (!is.na(tau) && tau >= 0.8) {
      
      category <- "Tissue-Specific"
      
      
    }
    else{
      category <- 'nn'
    }

    results$category[i] <- category
    results$max_tissue[i] <- names(expr)[which.max(expr)]
  }
  
  return(results)
}


# ---- Parameters -------------------------------------------------------------

ATLAS_PATH <- 'fruitdbase/results/output_bgeecall/presence_by_tissue.rds'
EXPRESSION_SCORE_PATH <- 'fruitdbase/results/output_bgeecall/ExpressionScoreByTissue.rds'
TPM_GENE_PATH  <- 'fruitdbase/results/output_bgeecall/tpm_gene.rds'
METADATA_PATH <- 'fruitdbase/results/input_bgeecall/metadata_filtered_20260911.tsv'
EXPRESSION_SCORE_DEV_PATH <- 'fruitdbase/results/output_bgeecall/ExpressionScoreByDev.rds'
ATLAS_BYDEV_PATH <- 'fruitdbase/results/output_bgeecall/presence_by_dev.rds'
OUTPUT_DATA <- "C:/Users/Lab.Mejora VIII/Desktop/Eva/proyectos/fruitdbase_pipeline/fruitdbase/pipeline_post_processing/cleanForGiga/output_data/"

# ---- Load data -------------------------------------------------------------

# Presence absence file
atlas <- readRDS(ATLAS_PATH)
atlas_bydev <- readRDS(ATLAS_BYDEV_PATH)

# Expression Score Matrix 

expressionScore <- readRDS(EXPRESSION_SCORE_PATH)
expressionScoreByDev <- readRDS(EXPRESSION_SCORE_DEV_PATH)

# TPM matrix
tpm_gene <-  readRDS(TPM_GENE_PATH)

# Metadata
metadata <- read_delim(METADATA_PATH, delim='\t')

# ---- Pre-processing  -------------------------------------------------------------

# Filter by genes expressed at least in one sample

genesexpressedOneSample <- atlas %>% filter(expressed) %>% pull(id) %>% unique()
expressionScore <- expressionScore[rownames(expressionScore) %in% genesexpressedOneSample, ]

genesexpressedOneSampledev <- atlas_bydev %>% filter(expressed) %>% pull(id) %>% unique()
expressionScoreByDev <- expressionScoreByDev[rownames(expressionScoreByDev) %in% genesexpressedOneSampledev, ]

# Calculate median tpm per tissue and developmental stage

# Filter  by genes expressed at least in one sample
tpm_matrix <- tpm_gene[rownames(tpm_gene) %in% genesexpressedOneSample, ]

# Set the labels
tissue_vec <- metadata$tissue_name[match(colnames(tpm_matrix), metadata$BioSample)]
dev_stage_vec <- paste0(metadata$tissue_name[match(colnames(tpm_matrix), metadata$BioSample)],"_",metadata$development_stage[match(colnames(tpm_matrix), metadata$BioSample)])

median_tpm_by_tissue <- sapply(unique(tissue_vec), function(t) {
  apply(
    tpm_matrix[, tissue_vec == t, drop = FALSE],
    1,
    median,
    na.rm = TRUE
  )
})

median_tpm_by_tissue_dev <- sapply(unique(dev_stage_vec), function(t) {
  apply(
    tpm_matrix[, dev_stage_vec == t, drop = FALSE],
    1,
    median,
    na.rm = TRUE
  )
})


# ---- Calculation of tissue-specific genes  -------------------------------------------------------------

# By Bgee approach

# Obtain Tau Index Matrix 
TauBgee <- apply(expressionScore, 1, function(x) {
  if (max(x) == 0) return(NA)
  sum(1 - (x / max(x))) / (length(x) - 1)
})

# Clasify by Tau Index  
resultsBgee <- classify_with_median( expressionScore,TauBgee)
resultsBgee$max_tissue <- gsub("Score_","",resultsBgee$max_tissue)

highlySpecificGeneBgee <- list()

for(tissue in unique(resultsBgee$max_tissue)){
  
  highlyGene <- resultsBgee$gene[resultsBgee$category == "Tissue-Specific" & resultsBgee$max_tissue == tissue]
    
  highlySpecificGeneBgee[[tissue]] <- highlyGene
  
}

# By TPM Median approach

# Obtain Tau Index Matrix 

TauMedianTPM <- apply(median_tpm_by_tissue, 1, function(x) {
  if (max(x) == 0) return(NA)
  sum(1 - (x / max(x))) / (length(x) - 1)
})


# Clasify by Tau Index  
resultsMedianTPM <- classify_with_median( median_tpm_by_tissue,TauMedianTPM)

resultsMedianTPM$max_tissue <- gsub("Score_","",resultsMedianTPM$max_tissue)

highlySpecificGeneMeanTPM <- list()

for(tissue in unique(resultsMedianTPM$max_tissue)){
  
  highlyGene <- resultsMedianTPM$gene[resultsMedianTPM$category == "Tissue-Specific" & resultsMedianTPM$max_tissue == tissue]
    
  highlySpecificGeneMeanTPM[[tissue]] <- highlyGene
  
}


saveRDS(highlySpecificGeneBgee, paste0(OUTPUT_DATA,"highlySpecificGeneBgee.rds"))
saveRDS(highlySpecificGeneMeanTPM, paste0(OUTPUT_DATA,"highlySpecificGeneMeanTPM.rds"))
saveRDS(median_tpm_by_tissue, paste0(OUTPUT_DATA,"median_tpm_by_tissue.rds"))
saveRDS(median_tpm_by_tissue_dev, paste0(OUTPUT_DATA,"median_tpm_by_tissue_dev.rds"))
