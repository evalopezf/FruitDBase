# =============================================================================
# Utils
# Functional analysis
#
# Input:  highlySpecificGeneBgee and highlySpecificGeneMeanTPM files obtained 
#         by running utils/specific_genes_detection.R.
#         atlas obtained from `snakemake post_processing`
#         
# Output: tpm_enrichment.rds,bgee_enrichment.rds
# =============================================================================


# ---- Libraries --------------------------------------------------------------

library(clusterProfiler)
library(dplyr)
library(org.At.tair.db)

# ---- Functions --------------------------------------------------------------


functional_analysis <- function(specificTissueList){

# Gene universe
fdr_filtered <- atlas %>% dplyr::filter(expressed)
universe_query <- unique(fdr_filtered$id)

universe_at <- map_df_gene %>%
  dplyr::filter(query_gene %in% universe_query) %>%
  pull(at_gene) %>%
  unique()

# Result list
enrichment_results <- list()
ego_list <- list()

# Iterate over tissues
for (tissue in names(specificTissueList)) {
  genes_especie <- specificTissueList[[tissue]]
  
  if (length(genes_especie) == 0) next
  
  # Map P. dulcis genes to their Arabidpopsis ortholog
  genes_at <- map_df_gene %>%
    dplyr::filter(query_gene %in% genes_especie) %>%
    pull(at_gene) %>%
    unique()
  cat("Tissue:", tissue, "- Genes mapeados a TAIR:", length(genes_at), "\n")
  

  ego <- enrichGO(
    gene          = genes_at,
    universe      = universe_at,
    OrgDb         = org.At.tair.db,
    keyType       = "TAIR",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,        
    qvalueCutoff  = 1,        
    minGSSize     = 3,        
    maxGSSize     = 1000,    
    readable      = FALSE
  )
  


  if (!is.null(ego) && nrow(ego@result) > 0) {
      ego_list[[tissue]] <- ego  

    table <- ego@result %>%
      dplyr::filter(p.adjust < 0.05) %>%
      arrange(p.adjust) %>%
      dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, Count, geneID, qvalue) 
    
    enrichment_results[[tissue]] <- table
  }
}
  return(list(
  enrichment_results = enrichment_results,
  ego_list = ego_list
))

}



# ---- Parameters --------------------------------------------------------------

map_df_gene <- readRDS("D:/PhD/proyectos/atlasExpresion/pipelineBgee_v2/data/map_df_gene.rds")
OUTPUT_DATA <- "C:/Users/Lab.Mejora VIII/Desktop/Eva/proyectos/fruitdbase_pipeline/fruitdbase/pipeline_post_processing/cleanForGiga/output_data/"
ATLAS_PATH <- 'fruitdbase/results/output_bgeecall/presence_by_tissue.rds'


# ---- Load Data --------------------------------------------------------------
highlySpecificGeneMeanTPM <- readRDS( paste0(INPUT_DATA,"highlySpecificGeneMeanTPM.rds"))
highlySpecificGeneBgee <- readRDS(paste0(INPUT_DATA,"highlySpecificGeneBgee.rds"))
atlas <- readRDS(ATLAS_PATH)

# ---- Functional analysis based on Arabidopsis Thaliana --------------------------------------------------------------

# Cal the function 
tpm <- functional_analysis(highlySpecificGeneMeanTPM)
tpm_enrichment <- tpm$enrichment_results

ego_list<- tpm$ego_list

bgee <- functional_analysis(highlySpecificGeneBgee)
bgee_enrichment <- bgee$enrichment_results


tpm_dev <- functional_analysis(highlySpecificGeneMeanTPM_by_dev)
tpm_enrichment_dev <- tpm_dev$enrichment_results
ego_list<- tpm_dev$ego_list

saveRDS(tpm_enrichment,paste0(OUTPUT_DATA,'tpm_enrichment.rds'))
saveRDS(bgee_enrichment,paste0(OUTPUT_DATA,'bgee_enrichment.rds'))