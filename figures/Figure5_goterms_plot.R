# =============================================================================
# Plotting GO semantic similarity networks of the top 15 significant GO terms 
# resulting from the functional enrichment analysis of tissue-specific genes for
#  each tissue
#
# Input:  Files obtained by running utils/functional_analysis.R (tpm_enrichment,bgee_enrichment)
# Output: figure_5_goclusterplot.svg
# =============================================================================

# ---- Libraries --------------------------------------------------------------
library(scales)
library(igraph)
library(ggrepel)
library(GOSemSim)
library(org.At.tair.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(ggforce)
library(patchwork)

# ---- Parameters --------------------------------------------------------------
INPUT_DATA <-  "C:/Users/Lab.Mejora VIII/Desktop/Eva/proyectos/fruitdbase_pipeline/fruitdbase/pipeline_post_processing/cleanForGiga/output_data/"
tpm_enrichment <- readRDS(paste0(INPUT_DATA,'tpm_enrichment.rds'))
bgee_enrichment <- readRDS(paste0(INPUT_DATA,'bgee_enrichment.rds'))

# Maximum terms by method
MAX_TERMS_TPM        <- 15

SEM_METHOD           <- "Wang" #similarity semantic method
EDGE_THRESH          <- 0.35
CLUSTER_CUTOFF       <- 0.7 # Dendogram cut off
MIN_GENES            <- 2
NET_SCALE            <- 9
FACET_NCOL           <- 5   
NODE_SIZE_MIN        <- 3.5
NODE_SIZE_MAX        <- 7
NODE_ALPHA           <- 0.92
MAX_EDGES_PER_TISSUE <- 30 
EDGE_ALPHA_RANGE     <- c(0.12, 0.55)
HULL_EXPAND          <- 0.85
HULL_ALPHA_BG        <- 0.10
HALO_COLOR           <- "#FFD600"
HALO_MULT            <- 1.8
SHOW_TERM_LABELS_GLOBAL  <- FALSE
CLUSTER_LABEL_SIZE       <- 5.5
STRIP_LABEL_SIZE         <- 13   
OUTPUT_PATH <- "C:/Users/Lab.Mejora VIII/Desktop/Eva/proyectos/fruitdbase_pipeline/fruitdbase/pipeline_post_processing/cleanForGiga/output_figures/"


tissue_colors <- c(
  Wood           = "#8D6E63",
  Anther         = "#E91E63",
  Seed           = "#FF6F00",
  Branch         = "#33691E",
  Phloem         = "#00ACC1",
  Pistil         = "#7B1FA2",
  Flower         = "#D50000",
  Fruit          = "#2E7D32",
  Flower_Bud     = "#AD1457",
  Leaf           = "#558B2F",
  Ovary          = "#1565C0",
  pedicel        = "#C2185B",
  Fruit_SeedCoat = "#F9A825"
)

# ---- Functions --------------------------------------------------------------

semData_file <- "semData_AtBP.rds"
if (file.exists(semData_file)) {
  semData <- readRDS(semData_file)
} else {
  semData <- godata("org.At.tair.db", keytype="TAIR", ont="BP", computeIC=FALSE ) # computeIC not use by Wnag method
  saveRDS(semData, semData_file)
}

# Create similarity matrix
get_sem_matrix <- function(go_ids, semData, method="Wang") {
  n   <- length(go_ids)
  mat <- diag(n); dimnames(mat) <- list(go_ids, go_ids)
  if (n < 2) return(mat)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      s <- tryCatch(goSim(go_ids[i], go_ids[j], semData=semData, measure=method),
                   error=function(e) NA_real_)
      s <- if (is.na(s) || !is.finite(s)) 0 else s
      mat[i,j] <- s; mat[j,i] <- s
    }
  }
  mat
}

# Cluster by similar terms
cluster_terms <- function(df, sim_mat, cutoff=CLUSTER_CUTOFF) {
  ok <- df$ID[df$ID %in% rownames(sim_mat)]
  if (length(ok) < 2)
    return(df %>% mutate(cluster=1L, cluster_label=Description[1]))
  sub <- sim_mat[ok, ok, drop=FALSE]
  sub[sub>1] <- 1; sub[sub<0] <- 0
  hc  <- hclust(as.dist(1-sub), method="average")
  cls <- cutree(hc, h=cutoff)
  df  <- df %>%
    left_join(data.frame(ID=ok, cluster=as.integer(cls),
                          stringsAsFactors=FALSE), by="ID") %>%
    mutate(cluster=replace_na(cluster, max(cls, na.rm=TRUE)+1L))
  nms <- df %>%
    group_by(cluster) %>%
    slice_min(p.adjust, n=1, with_ties=FALSE) %>%
    dplyr::select(cluster, cluster_label=Description) %>% ungroup()
  df %>% left_join(nms, by="cluster")
}

cluster_aware_layout <- function(df, sim_mat, scale=NET_SCALE) {
  n  <- nrow(df)
  if (n == 1) return(data.frame(x=0, y=0))
  clusters <- unique(df$cluster)
  n_cl     <- length(clusters)
  cl_cols  <- ceiling(sqrt(n_cl * 1.5))
  cl_cx    <- ((seq_len(n_cl)-1) %% cl_cols) * (scale * 0.55)
  cl_cy    <- -((seq_len(n_cl)-1) %/% cl_cols) * (scale * 0.55)
  all_x    <- numeric(n); all_y <- numeric(n)
  for (ci in seq_along(clusters)) {
    cl     <- clusters[ci]
    idx    <- which(df$cluster == cl)
    n_cl_i <- length(idx)
    if (n_cl_i == 1) { all_x[idx] <- cl_cx[ci]; all_y[idx] <- cl_cy[ci]; next }
    sub_ids    <- df$ID[idx]
    sub_ids_ok <- sub_ids[sub_ids %in% rownames(sim_mat)]
    if (length(sub_ids_ok) >= 2) {
      sub_mat <- sim_mat[sub_ids_ok, sub_ids_ok, drop=FALSE]
      g_sub   <- graph_from_adjacency_matrix(sub_mat, mode="undirected",
                                              weighted=TRUE, diag=FALSE)
      set.seed(42+ci)
      lay_sub <- layout_with_fr(g_sub, niter=500,
                                 weights=pmax(E(g_sub)$weight, 0.01))
    } else {
      ang     <- seq(0, 2*pi, length.out=n_cl_i+1)[-(n_cl_i+1)]
      lay_sub <- cbind(cos(ang), sin(ang))
    }
    rx <- range(lay_sub[,1]); ry <- range(lay_sub[,2])
    cr <- scale * 0.50
    if (diff(rx)>0) lay_sub[,1] <- (lay_sub[,1]-rx[1])/diff(rx)*cr*2-cr
    if (diff(ry)>0) lay_sub[,2] <- (lay_sub[,2]-ry[1])/diff(ry)*cr*2-cr
    for (k in seq_along(idx)) {
      pk <- if (sub_ids[k] %in% sub_ids_ok) which(sub_ids_ok==sub_ids[k]) else k
      all_x[idx[k]] <- cl_cx[ci] + lay_sub[min(pk, nrow(lay_sub)), 1]
      all_y[idx[k]] <- cl_cy[ci] + lay_sub[min(pk, nrow(lay_sub)), 2]
    }
  }
  data.frame(x=all_x, y=all_y)
}

# ---- Network  --------------------------------------------------------------

non_empty <- Filter(function(df) nrow(df) >= MIN_GENES, tpm_enrichment)
n_tis     <- length(non_empty)

all_nodes <- list()
all_edges <- list()

for (i in seq_along(names(non_empty))) {
  tname <- names(non_empty)[i]
  df    <- non_empty[[tname]]
  cat(sprintf("[%d/%d] %-12s\n", i, n_tis, tname))

  df_t <- df %>%
    arrange(p.adjust) %>% slice_head(n=MAX_TERMS_TPM) %>%
    mutate(tissue=tname,
           node_id=paste0(tname,"_",seq_len(n())),
           label_short=str_wrap(Description, 18),
           gene_count=Count,
           sz=rescale(Count, to=c(NODE_SIZE_MIN, NODE_SIZE_MAX), from=range(Count)))

  bgee_ids <- if (!is.null(bgee_enrichment[[tname]]) &&
                    nrow(bgee_enrichment[[tname]])>0)
    bgee_enrichment[[tname]]$ID else character(0)
  df_t <- df_t %>% mutate(in_bgee=ID %in% bgee_ids)

  sim_mat <- tryCatch(get_sem_matrix(df_t$ID, semData, SEM_METHOD),
                      error=function(e) NULL)

  df_t <- cluster_terms(df_t, sim_mat)

  # Delete cluster with one go terms only
  cluster_sizes <- table(as.character(df_t$cluster))
    
  valid_clusters <- names(cluster_sizes[cluster_sizes >= 2])

  df_t <- df_t %>%
    filter(as.character(cluster) %in% valid_clusters)

  if (nrow(df_t) == 0) {
    next
  }

  # Update matrix after filter
  sim_mat <- sim_mat[df_t$ID, df_t$ID, drop = FALSE]


  lay <- cluster_aware_layout(df_t, sim_mat, scale=NET_SCALE)
  rx    <- range(lay$x); ry <- range(lay$y)
  if (diff(rx)>0) lay$x <- (lay$x-rx[1])/diff(rx)*NET_SCALE
  if (diff(ry)>0) lay$y <- (lay$y-ry[1])/diff(ry)*NET_SCALE

  df_t$x <- lay$x; df_t$y <- lay$y

  ids <- df_t$ID; nids <- df_t$node_id
  edges_t <- data.frame(from=character(), to=character(),
                         weight=numeric(), same_cluster=logical(),
                         stringsAsFactors=FALSE)
  if (nrow(df_t) >= 2) {
    for (a in 1:(length(ids)-1)) {
      for (b in (a+1):length(ids)) {
        if (ids[a] %in% rownames(sim_mat) && ids[b] %in% rownames(sim_mat)) {
          w <- sim_mat[ids[a], ids[b]]
          if (!is.na(w) && w >= EDGE_THRESH)
            edges_t <- rbind(edges_t, data.frame(
              from=nids[a], to=nids[b], weight=w,
              same_cluster=df_t$cluster[a]==df_t$cluster[b],
              stringsAsFactors=FALSE))
        }
      }
    }
    if (nrow(edges_t) > MAX_EDGES_PER_TISSUE)
      edges_t <- edges_t %>% arrange(desc(weight)) %>% slice_head(n=MAX_EDGES_PER_TISSUE)
  }

  if (nrow(edges_t) > 0) {
    xy  <- df_t %>% dplyr::select(node_id, x, y)
    edf <- edges_t %>%
      left_join(xy, by=c("from"="node_id")) %>% dplyr::rename(x1=x, y1=y) %>%
      left_join(xy, by=c("to"="node_id"))   %>%  dplyr::rename(x2=x, y2=y) %>%
      mutate(tissue=tname)
    all_edges[[tname]] <- edf
  }
  all_nodes[[tname]] <- df_t
}



nodes_all  <- bind_rows(all_nodes)
edges_all  <- if (length(all_edges)>0) bind_rows(all_edges) else NULL
nodes_halo <- nodes_all %>% filter(in_bgee)
tc_vec     <- tissue_colors[names(tissue_colors) %in% nodes_all$tissue]


cluster_centers <- nodes_all %>%
  group_by(tissue, cluster, cluster_label) %>%
  summarise(
    cx = mean(x),
    cy = mean(y),
    n_nd = n(),
    .groups = "drop"
  ) %>%
  mutate(
    label_short = str_wrap(cluster_label, 16)
  )

tissue_order <- nodes_all %>%
  group_by(tissue) %>% summarise(tot=sum(gene_count), .groups="drop") %>%
  arrange(desc(tot)) %>% pull(tissue)
nodes_all$tissue       <- factor(nodes_all$tissue, levels=tissue_order)
if (!is.null(edges_all)) edges_all$tissue <- factor(edges_all$tissue, levels=tissue_order)
cluster_centers$tissue <- factor(cluster_centers$tissue, levels=tissue_order)

# ---- Plot  --------------------------------------------------------------


p <- ggplot() +

  ggforce::geom_mark_hull(
    data=nodes_all,
    aes(x=x, y=y, group=interaction(tissue, cluster),
        fill=tissue, color=tissue),
    expand=unit(HULL_EXPAND,"lines"),
    radius=unit(HULL_EXPAND*0.8,"lines"),
    alpha=HULL_ALPHA_BG, linewidth=0.55, linetype="solid",
    show.legend=FALSE
  ) +


  { if (!is.null(edges_all) && any(!edges_all$same_cluster))
    geom_segment(data=filter(edges_all, !same_cluster),
                 aes(x1,y1,xend=x2,yend=y2,color=tissue),
                 alpha=0.10, linewidth=0.22, lineend="round",
                 show.legend=FALSE)
  } +

  { if (!is.null(edges_all) && any(edges_all$same_cluster))
    geom_segment(data=filter(edges_all, same_cluster),
                 aes(x1,y1,xend=x2,yend=y2,color=tissue,alpha=weight),
                 linewidth=0.55, lineend="round", show.legend=FALSE)
  } +

  { if (nrow(nodes_halo)>0) list(
    geom_point(data=nodes_halo, aes(x,y,size=sz*HALO_MULT*1.3),
               color=HALO_COLOR, alpha=0.30, show.legend=FALSE),
    geom_point(data=nodes_halo, aes(x,y,size=sz*HALO_MULT),
               shape=21, fill=NA, color=HALO_COLOR,
               stroke=2.0, alpha=1, show.legend=FALSE)
  )} +

  geom_point(data=nodes_all, aes(x,y,color=tissue,size=sz), alpha=NODE_ALPHA) +
  geom_point(data=nodes_all, aes(x,y,size=sz*0.35), color="white", alpha=0.60) +

  { if (SHOW_TERM_LABELS_GLOBAL)
    geom_text_repel(
      data=nodes_all,
      aes(x,y,label=label_short,color=tissue),
      size=2.6, lineheight=0.80, fontface="plain",
      box.padding=0.25, point.padding=0.8,
      force=2, force_pull=0.08,
      max.overlaps=12,
      segment.size=0.15, segment.alpha=0.30,
      segment.curvature=0.1,
      show.legend=FALSE,
      min.segment.length=0.25, seed=42
    )
  } +
  geom_text_repel(
    data            = cluster_centers,
    aes(x=cx, y=cy, label=label_short, color=tissue),
    size            = CLUSTER_LABEL_SIZE,
    fontface        = "bold",
    lineheight      = 0.90,
    box.padding     = 1.0,
    point.padding   = 1.0,
    force           = 6,
    force_pull      = 0.05,
    max.overlaps    = 40,
    segment.size    = 0.35,
    segment.alpha   = 0.55,
    segment.color   = "grey40",
    segment.curvature = 0.15,
    arrow           = arrow(length=unit(0.014,"npc"), type="open"),
    show.legend     = FALSE,
    min.segment.length = 0.1,
    seed            = 123,
    bg.color        = "white",
    bg.r            = 0.18
  ) +

  scale_color_manual(values=tc_vec, guide="none") +  
  scale_fill_manual( values=tc_vec, guide="none") +
  scale_size_identity() +
  scale_alpha_continuous(range=EDGE_ALPHA_RANGE, guide="none") +
  facet_wrap(~tissue, ncol=FACET_NCOL, scales="free") +

  theme_minimal(base_size=14) +
  theme(
    panel.grid       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    axis.title       = element_blank(),
    panel.border     = element_rect(color="grey75", fill=NA, linewidth=0.5),
    panel.spacing    = unit(1.4, "lines"),
    strip.text       = element_text(face="bold", size=STRIP_LABEL_SIZE, color="grey15"),
    strip.background = element_rect(fill="grey93", color=NA),
    strip.placement  = "outside",
    plot.title       = element_text(face="bold", size=20, hjust=0.5,
                                    margin=margin(b=4), family="sans"),
    plot.subtitle    = element_text(size=11.5, hjust=0.5, color="grey30",
                                    margin=margin(b=10)),
    plot.caption     = element_text(size=9.5, color="grey45", hjust=0.5,
                                    margin=margin(t=8)),
    plot.background  = element_rect(fill="#FFFFFF", color=NA),
    panel.background = element_rect(fill="#FFFFFF", color=NA),
    plot.margin      = margin(18, 18, 10, 18)
  ) 

# Legend

cnt_all  <- nodes_all$gene_count
sz_brks  <- c(min(cnt_all), round(quantile(cnt_all, 0.5)), max(cnt_all))
sz_vals  <- rescale(sz_brks, to=c(NODE_SIZE_MIN, NODE_SIZE_MAX),
                    from=range(cnt_all))

p_leg_size <- ggplot(
  data.frame(x=1:3, y=1, sz=sz_vals, lab=as.character(sz_brks)),
  aes(x, y)) +
  geom_point(aes(size=sz), color="grey30") +
  geom_text(aes(label=lab, y=y-0.55), size=4, color="grey20", hjust=0.5) +
  scale_size_identity() +
  ylim(-0.2, 1.8) +
  annotate("text", x=2, y=1.7, label="Gene count",
           fontface="bold", size=4.2, color="grey20", hjust=0.5) +
  theme_void() +
  theme(plot.background=element_rect(fill="white", color="grey80",
                                     linewidth=0.5),
        plot.margin=margin(6,14,6,14))


p_final <- p / p_leg_size +
  plot_layout(heights=c(10, 0.7))

ggsave(paste0(OUTPUT_PATH,"figure_5_goclusterplot.svg"), plot = p_final, width = 25, height = 15)
