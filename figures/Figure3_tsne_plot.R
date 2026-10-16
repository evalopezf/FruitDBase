# =============================================================================
# t-SNE of RNA-seq samples by BioProject and tissue
#
# Log-normalizes gene counts, selects the top 5000 HVGs and runs t-SNE.
#
# Input:  metadata(tsv), gene_counts (from `snakemake post_processing`)
# Output: tsne_5000_hvg.png
# =============================================================================

# ---- Libraries --------------------------------------------------------------
library(readr)                 
library(dplyr)                 
library(SingleCellExperiment)  
library(scater)                
library(scran)                 
library(ggplot2)               
library(ggforce)               

# ---- Parameters -------------------------------------------------------------

METADATA <- read_delim('metadata_filtered_20260911.tsv') #Change this to your metadata
OUTPUT_PATH  <- 'C:/Users/Lab.Mejora VIII/Desktop/Eva/proyectos/fruitdbase_pipeline/fruitdbase/pipeline_post_processing/cleanForGiga/output/'
COUNTS_PATH <- 'counts_gene.rds' #Change this to your counts_gene file obtained from `snakemake post_processing` (post_processing/integration/counts_gene.rds)
N_HVG         <- 5000
PERPLEXITY    <- 30
SEED          <- 678 

# ---- Load Data ---------------------------------------------------------


est_counts <- readRDS(COUNTS_PATH)
est_counts <- as.data.frame(est_counts)


METADATA$external_id[METADATA$external_id == "unpublished"] <- "PRJNA1540450"
# Sample IDs
sample_ids <- colnames(est_counts)

sorted_metadata <- METADATA[
  match(sample_ids, METADATA$BioSample),
  ,
  drop = FALSE
]

if (!all(sorted_metadata$BioSample == sample_ids)) {
  stop("The samples in the expression matrix do not match the metadata.")
}

# SingleCellExperiment
sce <- SingleCellExperiment(
  assays = list(
    counts = as.matrix(est_counts)
  )
)

# Add metadata
colData(sce) <- S4Vectors::DataFrame(sorted_metadata)

# ---- Normalization ---------------------------------------------------------


sce <- logNormCounts(sce)

# Model mean-variance
dec <- modelGeneVar(sce)

plot_data <- data.frame(
  mean = rep(dec$mean, 2),
  variance = c(dec$total, dec$tech),
  type = rep(c("Total variance", "Technical variance"),
             each = nrow(dec))
)

# Mean-variance plot
mean_variance_plot <- ggplot(
  plot_data,
  aes(x = mean, y = variance, color = type)
) +
  geom_point() +
  scale_color_manual(
    values = c(
      "Total variance" = "black",
      "Technical variance" = "red"
    )
  ) +
  labs(
    x = "Mean log-counts",
    y = "Variance",
    color = NULL
  ) +
  theme_classic()

# ---- Highly variable genes ---------------------------------------------------------


# Select top 5000 high variables genes
top_hvgs <- getTopHVGs(
  dec,
  n = N_HVG
)

# ---- t-SNE plot ---------------------------------------------------------


# Set seed
set.seed(SEED)

# Call tsne 
sce <- runTSNE(
  sce,
  subset_row = top_hvgs,
  name = "TSNE",
  perplexity = PERPLEXITY  #Try with different perplexity values
)

# Extraer coordenadas
tsne_df <- as.data.frame(
  reducedDim(sce, "TSNE")
)

colnames(tsne_df) <- c(
  "TSNE1",
  "TSNE2"
)

tsne_df$Tissue <- colData(sce)$tissue_name
tsne_df$BioProject <- colData(sce)$external_id
tsne_df$BioSample <- colData(sce)$BioSample


# ---- Plot ---------------------------------------------------------


# Set colors


proj_pal <- c(
  "#4477AA", "#EE6677", "#228833", "#CCBB44",
  "#66CCEE", "#AA3377", "#BBBBBB", "#44AA99",
  "#CC6677", "#332288", "#117733", "#999933",
  "#882255", "#88CCEE", "#E69F00", "#56B4E9","#7F77DD"
)


# Set the BioProject ID to identify each project
bioproj_levels <- sort(unique(tsne_df$BioProject))
n_proj <- length(bioproj_levels)

tsne_df <- tsne_df %>%
  mutate(
    proj_label = factor(
      BioProject,
      levels = bioproj_levels
    )
  )

# Assign each color to a project

proj_color_map <- setNames(
  proj_pal[seq_len(n_proj)],
  levels(tsne_df$proj_label)
)


p_hull <- ggplot(
  tsne_df,
  aes(
    x = TSNE1,
    y = TSNE2
  )
) +

  geom_mark_hull(
    aes(
      fill = Tissue,
      color = Tissue,
      label = Tissue
    ),
    concavity = 3,
    expand = unit(2.5, "mm"),
    alpha = 0.08,
    linewidth = 0.4,
    label.buffer = unit(0.58, "mm"),
    label.fontsize = 7,
    label.colour = "grey30",
    label.fill = NA,
    con.colour = "grey70",
    con.size = 0.3,
    show.legend = FALSE
  ) +

  ggnewscale::new_scale_fill() +

  geom_point(
    aes(
      fill = proj_label
    ),
    shape = 21,
    color = "black",
    stroke = 0.3,
    size = 2.2,
    alpha = 0.85
  ) +

  scale_fill_manual(
    name = "BioProject",
    values = proj_color_map
  )


# ---- Save output ---------------------------------------------------------

ggsave(paste0(OUTPUT_PATH,"mean_variance_plot.png"), plot = mean_variance_plot, width = 10, height = 6, dpi = 300)
ggsave(paste0(OUTPUT_PATH,"tsne_500_hvg.png"), plot = p_hull, width = 10, height = 6, dpi = 300)
