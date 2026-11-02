# =============================================================================
# Plotting Expression stability in almond of putative orthologs of Arabidopsis reference genes.
#
# Input:  
#         median_tpm_by_tissue_dev (obtain by utils/specific_genes_detection.R)
# Output: Figure8_reference_genes.svg
# =============================================================================

# ---- Libraries --------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(grid)
# ---- Function -------------------------------------------------------------


a_largo <- function(mat, genes, tissues) {
  as.data.frame(mat) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "tissue", values_to = "value") %>%
    mutate(tissue = factor(tissue, levels = tissues),
          gene   = factor(gene, levels = names(genes), labels = genes),
          group  = factor(groups[as.character(tissue)],
                          levels = c("Male", "Female", "Whole flower", "Vegetative")))
}
panel <- function(df, ylab, y_limits, by) {
  ggplot(df, aes(x = tissue, y = value, fill = group)) +
    geom_col(width = 0.7) +
    facet_wrap(~ gene, ncol = 4, axes = "all_y") +
    scale_fill_manual(values = color, name = NULL) +
    scale_y_continuous(limits = y_limits,
                      breaks = seq(y_limits[1], y_limits[2], by = by),
                      expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Tissue", y = ylab) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_blank(),
          strip.text = element_text(face = "italic"),
          legend.position = "bottom")
}

# ---- Parameters -------------------------------------------------------------
INPUT_DATA <-  " "
OUTPUT_PATH <- " "
MIN_TPM_BACKGROUND   <- 5  
MIN_TPM_RECOMMENDED <- 10 
PCT_RECOMMENDED     <- 10  
PCT_SUITABLE        <- 25  
# ---- Load Data -------------------------------------------------------------
median_tpm_by_tissue_dev<-readRDS(paste0(INPUT_DATA,"median_tpm_by_tissue_dev.rds"))
map_df_gene <- "\pipeline\data\map_df_gene.rds"
# ---- USE CASE II  -------------------------------------------------------------

# Genes identified as stable from Czechowski et al. 2005
traditional <- c(
  AT4G05320 = "UBQ10", AT5G60390 = "EF-1a", 
  AT1G13440 = "GAPDH",
  AT5G25760 = "UBC",  AT3G18780 = "ACT2",
  AT1G07920 = "EF-1a",  AT1G07930 = "EF-1a", AT1G07940 = "EF-1a"
)

novel <- c(
  AT5G46630 = "Clathrin adaptor", AT5G08290 = "YLS8",
  AT4G33380 = "Expressed",        AT2G32170 = "Expressed",
  AT4G34270 = "TIP41-like",      AT5G53300 = "UBC10",
  AT4G27960 = "UBC9",
  AT1G13320 = "PDF2",      AT2G28390 = "SAND",
  AT4G26410 = "Expressed",        AT3G01150 = "PTB",
  AT5G55840 = "PPR",              AT1G58050 = "Helicase",
  AT3G53090 = "UPL7",             AT5G15710 = "F-box",
  AT4G38070 = "bHLH",             AT5G12240 = "Expressed",
  AT1G62930 = "PPR"
)

# ---- Stability of all genes ---------------------------------------------------

m <- as.matrix(median_tpm_by_tissue)

stab <- tibble(
  query_gene = rownames(m),
  mean_tpm   = rowMeans(m),
  min_tpm    = apply(m, 1, min),
  cv         = apply(m, 1, sd) / rowMeans(m),
  mfc        = (apply(m, 1, max) + 1) / (apply(m, 1, min) + 1) 
) %>%
  filter(mean_tpm > 0)

# Background, genes with more than MIN_TPM_BACKGROUND in at least one tissue
background <- stab %>% filter(mean_tpm >= MIN_TPM_BACKGROUND)
pct_fun    <- ecdf(background$cv)

stab <- stab %>% mutate(cv_percentile = round(100 * pct_fun(cv), 1)) 

# ---- Almond orthologs of the reference genes ----------------------------------
ref <- bind_rows(
  tibble(at_gene = names(traditional), name = traditional, origin = "Traditional"),
  tibble(at_gene = names(novel),    name = novel,    origin = "Novel")
) %>%
  inner_join(map_df_gene, by = "at_gene") %>%      # all almond orthologs
  inner_join(stab, by = "query_gene") %>%
  mutate(
    category = case_when(
      cv_percentile <= PCT_RECOMMENDED & min_tpm >= MIN_TPM_RECOMMENDED ~ "Recommended",
      cv_percentile <= PCT_SUITABLE                                     ~ "Suitable",
      TRUE                                                              ~ "Not recommended"
    ),
    category = factor(category, levels = c("Recommended", "Suitable", "Not recommended")),
    gene_lab = ifelse(name == "Expressed", at_gene, name),   
    label    = paste0(gene_lab, " (", query_gene, ")")
  ) %>%
  distinct(at_gene, cv, mean_tpm, .keep_all = TRUE) %>%   
  arrange(cv_percentile)

  
# ---- Summary by Arabidopsis gene ----------------------------------------------
# A reference gene is supported in almond when at least one ortholog is in the
# 25% most stable genes
per_agi <- ref %>%
  group_by(at_gene, name, origin) %>%
  summarise(n_orthologs   = n(),
            orthologs     = paste(query_gene, collapse = ", "),
            best_ortholog = query_gene[which.min(cv_percentile)],
            best_pct      = min(cv_percentile),
            .groups = "drop") %>%
  mutate(candidate = best_pct <= PCT_SUITABLE) %>%
  arrange(best_pct)

print(per_agi, n = Inf)
per_agi %>% count(origin, candidate)


# ---- Figure: stability percentile ---------------------------------------------
cat_colors <- c("Recommended"     = "#1C578A",
                "Suitable"        = "#98C1D9",
                "Not recommended" = "#8B2E3C")

cut_labels <- tibble(x = c(PCT_RECOMMENDED, PCT_SUITABLE),
                    lab = paste0(c(PCT_RECOMMENDED, PCT_SUITABLE), "%"))

fig <- ggplot(ref, aes(x = cv_percentile, y = reorder(label, -cv_percentile),
                      colour = category)) +
  geom_vline(xintercept = cut_labels$x, linetype = "dashed", colour = "grey60") +
  geom_text(data = cut_labels, aes(x = x, y = Inf, label = lab), inherit.aes = FALSE,
            vjust = -0.5, size = 3, colour = "grey40") +
  geom_segment(aes(x = 0, xend = cv_percentile, yend = reorder(label, -cv_percentile))) +
  geom_point(aes(shape = origin), size = 2.5) +
  scale_colour_manual(values = cat_colors, name = NULL) +
  scale_shape_manual(values = c("Traditional" = 17, "Novel" = 16), name = NULL) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "Stability percentile (CV across tissues)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", legend.box = "vertical",
        plot.margin = margin(15, 10, 5, 5))


ggsave(paste0(OUTPUT_PATH, "Figure8_reference_genes.svg"), fig, width = 17, height = 14, units = "cm")
