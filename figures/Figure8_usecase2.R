# =============================================================================
# Plotting DAM genes expression acrros flower bud and wood developmental stages
#
# Input:  expression score matrix (obtained from `snakemake post_processing`),
#         median_tpm_by_tissue_dev (obtain by utils/specific_genes_detection.R)
# Output: dam_genes_expression_pattern.svg
# =============================================================================

# ---- Libraries --------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(grid)
# ---- Function -------------------------------------------------------------

a_largo <- function(mat, cols, tejido) {
  as.data.frame(mat[names(dams), cols, drop = FALSE]) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "stage_raw", values_to = "value") %>%
    mutate(developmental_stage  = factor(stage_raw, levels = cols, labels = etiquetas),
           gene   = factor(gene, levels = names(dams), labels = etiq_gen),
           tissue = tejido)
}


panel <- function(df, ylab) {
  ggplot(df, aes(x = developmental_stage, y = value, group = tissue, colour = tissue)) +
    geom_line(linewidth = 1) +                 # antes 0.6
    geom_point(size = 2.2, stroke = 0) +       # antes 1.3
    facet_wrap(~ gene, ncol = 4, scales = "free_y") +
    scale_colour_manual(values = color, name = NULL) +
scale_y_continuous(trans = scales::log1p_trans(),   # transform_log1p() en scales ≥ 1.3
                   breaks = c(0, 1, 10, 100, 1000),
                   labels = c("0", "1", "10", "100", "1000")) +
labs(x = "Developmental Stage", y = ylab)+
    theme_classic(base_size = 12) +            # antes 8
    theme(
      strip.background = element_blank(),
      strip.text       = element_text(size = 11, hjust = 0, face = "bold"),  # antes 7.5
      axis.line        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks       = element_line(colour = "black", linewidth = 0.4),
      axis.text        = element_text(colour = "black", size = 10),         # antes 6.5
      axis.text.x      = element_text(angle = 45, hjust = 1),
      legend.position  = "top",
      legend.text      = element_text(size = 11),
      legend.key.size  = unit(12, "pt"),
      panel.spacing    = unit(10, "pt")
    )
}

# ---- Parameters -------------------------------------------------------------
INPUT_DATA <-  "C:/Users/Lab.Mejora VIII/Desktop/Eva/proyectos/fruitdbase_pipeline/fruitdbase/pipeline_post_processing/cleanForGiga/output_data/"
OUTPUT_PATH <- "C:/Users/Lab.Mejora VIII/Desktop/Eva/proyectos/fruitdbase_pipeline/fruitdbase/pipeline_post_processing/cleanForGiga/output_figures/"
EXPRESSION_SCORE_DEV_PATH <- 'fruitdbase/results/output_bgeecall/ExpressionScoreByDev.rds'

# ---- Load Data -------------------------------------------------------------
median_tpm_by_tissue_dev<-readRDS(paste0(INPUT_DATA,"median_tpm_by_tissue_dev.rds"))
expressionScoreByDev <- readRDS(EXPRESSION_SCORE_DEV_PATH)

# ---- USE CASE II  -------------------------------------------------------------

# Genes identified by BLAST using DAM sequences from Bielenberg et al. (2008)
dams <- c(
  Prudul26A015462 = "DAM1",
  Prudul26A007437 = "DAM2",
  Prudul26A007027 = "DAM3/4",
  Prudul26A019427 = "DAM5/6"

)

# Set labels
etapas_bud  <- c("Flower_Bud_T0","Flower_Bud_T1","Flower_Bud_T2",
                 "Flower_Bud_T3","Flower_Bud_T4","Flower_flowering_stage")
etapas_wood <- c("Wood_T0","Wood_T1","Wood_T2",
                 "Wood_T3","Wood_T4","Wood_flowering_stage")
etiquetas   <- c("T0","T1","T2","T3","T4","Flowering")

etiq_gen <- setNames(paste0(dams, "\n(", names(dams), ")"), names(dams))

# Prepare the format to the plot
df_tpm <- bind_rows(
  a_largo(median_tpm_by_tissue_dev, etapas_bud,  "Flower bud"),
  a_largo(median_tpm_by_tissue_dev, etapas_wood, "Wood")) %>%
  mutate(value = log2(value + 1))          # <- TPM en log2

colnames(expressionScoreByDev) <- gsub('Score_','',colnames(expressionScoreByDev))
df_bgee <- bind_rows(
  a_largo(expressionScoreByDev, etapas_bud,  "Flower bud"),
  a_largo(expressionScoreByDev, etapas_wood, "Wood"))


color <- c("Flower bud" = "#1C578A", "Wood" = "#EBD152")

fig <- panel(df_tpm,  expression(log[2] *"(median TPM + 1)")) /
       panel(df_bgee, "Bgee expression score") +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "top",
        plot.tag = element_text(face = "bold", size = 12))

ggsave(paste0(OUTPUT_PATH,"dam_genes_expression_pattern.svg"), plot = fig, width = 15, height = 10)


