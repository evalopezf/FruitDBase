# =============================================================================
# Distribution Plot
#
# Geographic and biological characterization of 305 RNA-Seq samples from 
# Prunus dulcis available in NCBI SRA (A) Global distribution of projects 
# (B) Projects by country and number of samples. (C-D) Sample distribution 
# per tissue and cultivar.
#
# Input:  metadata(tsv)
# Output: plotDistributionMetadata.svg
# =============================================================================

# ---- Libraries --------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)
library(forcats)
library(openxlsx)            

# ---- Parameters -------------------------------------------------------------

# Fill the parameters with our own data

METADATA_PATH <- "curated_20260902.xlsx"  
OUTPUT_PATH <- "output/" 

# ---- Load Data -------------------------------------------------------------

metadata <-  read.xlsx(METADATA_PATH)
# Correct a USA name 
metadata <- metadata %>%
  mutate(origin_country = case_when(
    origin_country == "USA" ~ "United States of America",
    TRUE ~ origin_country
  ))


# ---- Plot A: World Map -------------------------------------------------------------

project_summary <- metadata %>%
  filter(!is.na(origin_country) & !is.na(external_id)) %>%
  group_by(origin_country) %>%
  summarise(
    n_projects = n_distinct(external_id),
    n_samples = n()
  ) %>%
  arrange(desc(n_projects))


world <- ne_countries(scale = "medium", returnclass = "sf")

world_data <- world %>%
  left_join(project_summary, by = c("name" = "origin_country"))

plot_map <- ggplot(data = world_data) +
  geom_sf(aes(fill = n_projects), color = "white", size = 0.2) +
  scale_fill_gradient(low = "#FFF3B0", high = "#274703", 
                      na.value = "grey95",
                      name = "Number of projects",
                      breaks = scales::pretty_breaks(n = 5)) +
  theme_minimal() +
  labs(title = "Global distribution of projects") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(1.5, "cm"),
    panel.grid = element_line(color = "grey90", size = 0.2),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# ---- Plot B: Porjects per country  -------------------------------------------------------------

plot_projects <- ggplot(project_summary, 
                        aes(x = reorder(origin_country, n_projects), 
                            y = n_projects)) +
  geom_bar(stat = "identity", fill = "#274703", alpha = 0.8) +
  geom_text(aes(label = paste0(n_projects, " (", n_samples, ")")), 
            hjust = -0.1, size = 3.5) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Samples and projects per country",
       subtitle = "Format: Projects (samples)",
       x = "Country",
       y = "Number of projects") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        panel.grid.minor = element_blank()) +
  expand_limits(y = max(project_summary$n_projects) * 1.3)

# ---- Plot C: Samples by tissue  -------------------------------------------------------------

tissue_count <- metadata %>%
  dplyr::count(tissue_name)

plot_tissue <- tissue_count %>%
  arrange(n) %>%
  ggplot(aes(y = fct_reorder(tissue_name, n), x = n)) +
  geom_bar(stat = "identity", fill = "#274703", alpha = 0.8) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.5, fontface = "bold") +
  theme_minimal() +
  labs(title = "Tissue distribution",
       y = "Tissue",
       x = "Number of samples") +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        panel.grid.minor = element_blank()) +
  expand_limits(x = max(tissue_count$n) * 1.1)

# ---- Plot D: Samples by cultivar  -------------------------------------------------------------

plot_cultivar <- metadata %>%
  dplyr::count(accession_name) %>%
  ggplot(
    aes(
      y = fct_reorder(accession_name, n),
      x = n,
      fill = accession_name == "Missing_pdu"
    )
  ) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.5, fontface = "bold") +
  scale_fill_manual(
    values = c(
      "FALSE" = "#274703",  
      "TRUE"  = "grey70"    # Missing values
    ),
    guide = "none"
  ) +
  theme_minimal() +
  labs(
    title = "Cultivar distribution",
    y = "Accession",
    x = "Number of samples"
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.grid.minor = element_blank()
  ) +
  expand_limits(x = max(dplyr::count(metadata, accession_name)$n) * 1.1)



# ---- Final Layout -------------------------------------------------------------
panel_superior <- plot_grid(
  plot_map,
  plot_projects,
  ncol = 2,
  labels = c("A", "B"),
  rel_widths = c(2, 1)  
)

panel_inferior <- plot_grid(
  plot_tissue,
  plot_cultivar,
  ncol = 2,
  labels = c("C", "D"),
  rel_widths = c(1, 1)
)

figura_final <- plot_grid(
  panel_superior,
  panel_inferior,
  ncol = 1,
  rel_heights = c(1.2, 1)
)

# ---- Output plot -------------------------------------------------------------


ggsave(
  paste0(OUTPUT_PATH,"plotDistributionMetadata.svg"),
  plot = figura_final,
  width = 37,
  height = 24,
  units = "cm",
  bg = "white"
)