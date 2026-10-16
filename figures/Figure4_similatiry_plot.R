# =============================================================================
# Plotting overlap coefficient between bgee and median TPM approach
# Plotting correlation between number of samples and number of specific genes detected
#
# Input:  Files obtained by running utils/specific_genes_detection.R
# Output: overlap_coefficient_plot.svg, relation_num_sample_genes_detected.svg, combined_plot.svg
# =============================================================================

# ---- Libraries --------------------------------------------------------------
library(readr)                 
library(dplyr)                 
library(tidyr)
library(ggplot2)               
library(patchwork)


# ---- Parameters -------------------------------------------------------------

INPUT_DATA <-  "C:/Users/Lab.Mejora VIII/Desktop/Eva/proyectos/fruitdbase_pipeline/fruitdbase/pipeline_post_processing/cleanForGiga/output_data/"
OUTPUT_PATH <- "C:/Users/Lab.Mejora VIII/Desktop/Eva/proyectos/fruitdbase_pipeline/fruitdbase/pipeline_post_processing/cleanForGiga/output_figures/"

# ---- Load data -------------------------------------------------------------


highlySpecificGeneBgee <- readRDS(paste0(INPUT_DATA,"highlySpecificGeneBgee.rds"))
highlySpecificGeneMeanTPM<- readRDS( paste0(INPUT_DATA,"highlySpecificGeneMeanTPM.rds"))
median_tpm_by_tissue<-readRDS(paste0(INPUT_DATA,"median_tpm_by_tissue.rds"))
median_tpm_by_tissue_dev<-readRDS(paste0(INPUT_DATA,"median_tpm_by_tissue_dev.rds"))


# ---- Szymkiewicz–Simpson coefficient  -------------------------------------------------------------
overlap_coef <- function(A, B) {
  length(intersect(A, B)) / min(length(A), length(B))
}

tissues <- names(highlySpecificGeneBgee)

df <- data.frame(
  Tissue = tissues,
  Overlap = sapply(tissues, function(t) {
    overlap_coef(highlySpecificGeneBgee[[t]],
                 highlySpecificGeneMeanTPM[[t]])
  })
)

# ---- Similarity Plot  -------------------------------------------------------------

bar_simpson <- ggplot(df, aes(x = Tissue, y=Overlap)) +
  geom_col(fill = "#9F9FA9") +
  coord_cartesian(ylim = c(0.25, 1.0))+

  labs(x = "Tissues", y = "Szymkiewicz–Simpson coefficient") +
  theme_minimal()+
    theme(
    axis.line.x = element_line(),
    axis.line.y = element_line(),
        axis.text.x = element_text(angle = 45, hjust = 1, size=16),
        axis.text.y = element_text(
      size = 16
    ),
    axis.title.y = element_text(
      size = 17
    )
  )

# ---- Comparative betwen number of samples and tissue-specific genes detected  -------------------------------

# Collect data in a dataframe
summary_expression <- data.frame(
  Tissue  = unique(atlas$tissue_name),
  Specific_gene_bgee = NA_integer_,
  '% Specific_genes_bgee' = NA_integer_,
    Specific_gene_tpm = NA_integer_,
  '% Specific_genes_tpm' = NA_integer_,
  Overlap = NA_integer_
)


for ( tissue in unique(atlas$tissue_name)){
  n_samples <- length(metadata$BioSample[metadata$tissue_name == tissue])
  bgee <- length(highlySpecificGeneBgee[[tissue]])
  tpm <- length(highlySpecificGeneMeanTPM[[tissue]])
  bgee_norm <- bgee / n_samples
  tpm_norm  <- tpm / n_samples
  total <- summary_expression$Expressed_genes[summary_expression$Tissue == tissue]
  summary_expression$Specific_gene_bgee[summary_expression$Tissue == tissue] <- bgee
  summary_expression$Specific_gene_tpm[summary_expression$Tissue == tissue] <- tpm
  summary_expression$n_samples[summary_expression$Tissue == tissue] <- length(unique(metadata$BioSample[metadata$tissue_name == tissue]))

}


summary_long <- summary_expression %>%
  dplyr::select(Tissue, n_samples,
         Specific_gene_bgee,
         Specific_gene_tpm) %>%
  pivot_longer(
    cols = c(Specific_gene_bgee, Specific_gene_tpm),
    names_to = "Method",
    values_to = "Specific_genes"
  ) %>%
  mutate(
    Method = recode(Method,
                    Specific_gene_bgee = "BGEE method",
                    Specific_gene_tpm  = "TPM method")
  )

method_labels <- c(
  "BGEE method" = "BGEE method",
  "TPM method"  = "TPM method"
)


shape_values <- c(16, 17, 15, 18, 1, 2, 0, 5, 8, 6, 4,3,9)

# Spearman calculation
cor_labels <- summary_long %>%
  group_by(Method) %>%
  summarise(
    rho = cor(n_samples, Specific_genes, method = "spearman"),
    p   = cor.test(n_samples, Specific_genes, method = "spearman")$p.value
  ) %>%
  mutate(label = paste0("Spearman's ρ = ", round(rho, 3),
                        "\n", "p-value = ", signif(p, 3)))


# ---- Comparative Plot  -------------------------------

specific_cor <- ggplot(summary_long, aes(n_samples, Specific_genes)) +
  geom_smooth(method = "lm", se = TRUE,
              color = "black", fill = "grey80") +
  geom_point(aes(shape = Tissue),
             size = 3.5,
             color = "black",
             fill = "white",
             stroke = 0.8) +
  scale_shape_manual(values = shape_values) +
  facet_wrap(~ Method, ncol = 1) +

  geom_label(data = cor_labels,
             aes(x = Inf, y = Inf, label = label),
             hjust = 1.08,   
             vjust = 1.2,  
             size = 5.0,
             fill = "white",
             color = "black",
             label.size = 0.5) +

  theme_bw( ) + 

theme(
  plot.title      = element_text(size = 20, face = "bold"),
  axis.title      = element_text(size = 16),
  axis.text       = element_text(size = 14),
  strip.text      = element_text(size = 16, face = "bold"),  
  legend.title    = element_text(size = 14),
  legend.text     = element_text(size = 12)
) +
    
labs(
  x = "Number of samples",
  y = "Number of tissue-specific genes"
)


combined_plot <- bar_simpson + specific_cor 

ggsave(paste0(OUTPUT_PATH,"similarity_plot.svg"), plot = bar_simpson, width = 25, height = 15)
ggsave(paste0(OUTPUT_PATH,"comparative_detection_plot.svg"), plot = specific_cor, width = 25, height = 15)
ggsave(paste0(OUTPUT_PATH,"figure_4_combined_plot.svg"), plot = combined_plot, width = 25, height = 15)
