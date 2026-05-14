library(rhdf5)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(umap)
library(readr)
library(here)
library(vroom)
library(tidytext)
library(cowplot)
library(tidygraph)
library(igraph)
library(vegan)
require(ggraph)
library(ggrastr)

###############################################################################!
# Functions ####
###############################################################################!
# Convert proband_id to a numeric ID within that gene_disease pair
createProbandLabel <- function(df){
  df %>% 
    arrange(gene, MONDO, proband_id) %>% 
    mutate(label = as.numeric(factor(proband_id)), .by = c(gene, MONDO)) %>% 
    unite(label, gene, MONDO, label, sep = "__")
}


#' geneHPOgraph
#'
#' @param g igraph object
#' @param gene String. Gene to be plotted
#' @param disease String, optional. Additional regex string to apply as an additional filter.
#' Used for genes that are associated with more than one disease
#' @param n_node_labels Double. Number of nodes to label. If NULL, no labels plotted. Labels plotted in order of term count
#' @param greater_than_1 Logical. Whether an HPO term must be present in more than 1 proband to be be included in n_node_labels
#' @param subtitle String, optional. Additional custom title text
#'
#' @return
#' @export
#'
#' @examples
geneHPOgraph <- function(g, gene, disease = NULL, n_node_labels = NULL, greater_than_1 = FALSE, subtitle = NULL, include_title = T, wrap_width = 30){
  # browser()
  # browser()
  # Get counts of HPOs associated with gene +/- disease
  hpo_count <- df_proband_hpo[df_proband_hpo$gene == gene,]
  if (!is.null(disease)){
    hpo_count <- hpo_count[grepl(disease, hpo_count$disease, ignore.case = T),]
  }
  hpo_count <- deframe(count(hpo_count, hpo_id, sort = T) )
  
  # Add HPO counts to nodes
  node_names <- V(g)$name
  node_count <- unname(hpo_count[node_names])
  node_count <- ifelse(is.na(node_count),0,node_count)
  V(g)$count <- node_count
  
  # Add HPO terms to nodes with counts above a threshold
  if (!is.null(n_node_labels)){
    node_term <- hpo_key[node_names]
    node_term <- stringr::str_wrap(node_term, width = wrap_width)
    threshold <- sort(node_count, decreasing = T)[n_node_labels]
    min_threshold <- ifelse(greater_than_1, 2, 1) # Threshold must be at least a count of 1, but at least 2 if specified
    threshold <- ifelse(threshold < min_threshold, min_threshold, threshold) 
    label <- ifelse(node_count >= threshold, node_term, NA)
    V(g)$label <- label
  }
  
  # Generate base plot
  plt <- ggraph(g, layout = "kk")+
    rasterize(geom_edge_link(alpha = 0.5, color = "grey40"), dpi = 500)+
    rasterize(geom_node_point(aes(size = log2(count), fill = log2(count)), shape = 21), dpi = 150)+
    theme_graph() +
    scale_fill_viridis_c() +
    theme(legend.position="bottom")+
    guides(size = "none")
  
  if (include_title){plt <- ggtitle(gene, subtitle = subtitle)}
  
  # Add labels to plot
  if (!is.null(n_node_labels)){    
    # Adjust plot axis range for label space
    x_range <- layer_scales(plt)$x$range$range
    x_max <- max(x_range)
    expanded_x_max <- x_max*2.5
    label_x_max <- expanded_x_max*0.4
    
    plt <- plt +
      expand_limits(x = expanded_x_max) +
      geom_node_label(
        aes(label = label),
        repel = T,
        size = 3,
        box.padding = 0.3,
        point.padding = 0.3,
        max.overlaps = Inf,
        min.segment.length = 0,
        na.rm = TRUE,
        xlim = c(label_x_max, NA),
        force = 10,
        force_pull = 2
      )
  }
  return(plt)
}

#' termHPOgraph
#'
#' @param g igraph object
#' @param gene String. Gene to be plotted
#' @param disease String, optional. Additional regex string to apply as an additional filter.
#' Used for genes that are associated with more than one disease
#' @param n_node_labels Double. Number of nodes to label. If NULL, no labels plotted. Labels plotted in order of term count
#' @param subtitle String, optional. Additional custom title text
#'
#' @return
#' @export
#'
#' @examples
termHPOgraph <- function(g, terms, label_nodes = TRUE, title = NULL){
  # browser()
  
  # Add HPO counts to nodes
  node_names <- V(g)$name
  V(g)$count <- ifelse(V(g)$name %in% terms, 1, NA)
  
  
  # Add HPO terms to nodes with counts above a threshold
  if (label_nodes){
    node_term <- hpo_key[node_names]
    node_term <- stringr::str_wrap(node_term, width = 30)
    label <- ifelse(node_names %in% terms, node_term, NA)
    V(g)$label <- label
  }
  
  # Generate base plot
  plt <- ggraph(g, layout = "kk")+
    rasterize(geom_edge_link(alpha = 0.5, color = "grey40"), dpi = 500)+
    rasterize(geom_node_point(aes(size = count), shape = 21, fill = "darkorange"), dpi = 150)+
    theme_graph() +
    scale_fill_viridis_c() +
    ggtitle(title)+
    theme(legend.position="bottom")+
    guides(size = "none")
  
  # Add labels to plot
  if (label_nodes){    
    # Adjust plot axis range for label space
    x_range <- layer_scales(plt)$x$range$range
    x_max <- max(x_range)
    expanded_x_max <- x_max*2.5
    label_x_max <- expanded_x_max*0.4
    
    plt <- plt +
      expand_limits(x = expanded_x_max) +
      geom_node_label(
        aes(label = label),
        repel = T,
        size = 3,
        box.padding = 0.3,
        point.padding = 0.3,
        max.overlaps = Inf,
        min.segment.length = 0,
        na.rm = TRUE,
        xlim = c(label_x_max, NA),
        force = 10,
        force_pull = 2
      )
  }
  return(plt)
}

###############################################################################!
# Load data ####
###############################################################################!

df_hpo_adjacency <- read_csv(here("data/scid_hpo_adjacency.csv"))

proband_path <- here("data/clingen/api/scid_hpo_funSimAvg.h5")
df_proband_hpo <- data.frame(h5read(proband_path, "/proband_hpo")) 

hpo_key <- read_csv(here("data/hpo_key.csv"))
hpo_key <- deframe(hpo_key)

h5_path_allClingen <- here("data/clingen/scrape/clingen_scrape_hpo.h5")
mtx_proband_allClingen <- h5read(h5_path_allClingen, "/proband_distance")
df_proband_allClingen <- data.frame(h5read(h5_path_allClingen, "/proband_metadata")) %>% 
  unite(label, gene, disease, proband_id, sep = "__", remove = F)
df_proband_allClingen <- createProbandLabel(df_proband_allClingen)
rownames(mtx_proband_allClingen) <- df_proband_allClingen$label
colnames(mtx_proband_allClingen) <- df_proband_allClingen$label

df_cdwg <- read_csv(here("data/cdwg_key.csv"))
df_gcep <- read_csv(here("data/gcep_key.csv")) %>% 
  mutate(
    GCEP = gsub(" GCEP", "", GCEP),
    GCEP = gsub("Hemostasis Thrombosis", "Hemostasis/Thrombosis", GCEP),
    GCEP = gsub("Charcot-Marie-Tooth", "Charcot-Marie-Tooth Disease", GCEP),
    GCEP = gsub("Monogenic Autoinflammatory Disease", "Monogenic Autoinflammatory Diseases", GCEP),
  )

df_key <- full_join(df_gcep, df_cdwg, by = "GCEP") %>% 
  mutate(CDWG = ifelse(is.na(CDWG), "Other", CDWG))

df_key <- df_key %>% 
  select(Gene, MONDO, GCEP) %>% 
  distinct() %>% 
  drop_na()

iuis_files <- list.files(here("data/gene_tables"), 
                         pattern = "IUIS_1", 
                         full.names = TRUE)
table1_genes <- lapply(iuis_files, read_lines)
table1_genes <- unlist(table1_genes)

gene_metadata_path <- here("data/clingen/api/scid_hpo_s20201101_e20251231_funSimAvg.h5")
df_gene_metadata <- data.frame(h5read(gene_metadata_path, "/gene_metadata"))

df_gene <- df_key %>% 
  filter(GCEP == "SCID-CID") %>% 
  mutate(count = n(), .by = Gene) %>% 
  filter(count == 1) %>% 
  left_join(df_gene_metadata, by = c("Gene" = "gene")) 

###############################################################################!
# HPO consistency ####
###############################################################################!

scid_probands <- df_proband_allClingen %>% 
  separate(label, into = c("Gene", "MONDO", "Proband"), sep = "__", remove = F) %>% 
  left_join(filter(df_key, GCEP == "SCID-CID"), by = c("Gene", "MONDO")) %>% 
  pull(GCEP) %>% 
  {. == "SCID-CID"} %>% 
  which()

df_clingen_hpo <- vroom::vroom(here("data/clingen/scrape/clingen_proband_hpo_terms.csv.gz"))

df_scid_dist <- broom::tidy(as.dist(mtx_proband_allClingen[scid_probands, scid_probands]))

df_scid_dist <- df_scid_dist %>% 
  separate(item1, into = c("gene1", "mondo1", "proband1"), sep = "__") %>% 
  separate(item2, into = c("gene2", "mondo2", "proband2"), sep = "__")

df_scid_plot <- df_scid_dist %>% 
  mutate(distance = 1-distance) %>% 
  filter(gene1 == gene2 & mondo1 == mondo2) %>% 
  filter(gene1 %in% table1_genes) %>% 
  unite(label, gene1, mondo1) %>% 
  left_join(
    df_gene %>% 
      unite(label, Gene, MONDO) %>% 
      select(label, classification),
    by = c("label")
  ) %>% 
  mutate(classification = factor(classification, levels = rev(c("Limited", "Moderate", "Strong", "Definitive")))) 

scid_groups <- df_scid_plot %>% 
  unite(label2, gene2, mondo2, sep = "_") %>% 
  select(label, label2) %>% 
  distinct() %>% 
  pivot_longer(everything(), names_to = "label", values_to = "value") %>% 
  select(value) %>% 
  distinct() %>% 
  pull(value)

df_classification_key <- df_scid_plot %>% 
  select(label, classification) %>% 
  distinct()

df_scid_plot_jaccard <- df_clingen_hpo %>% 
  unite(group, Gene, MONDO, sep = "_", remove = F) %>% 
  filter(group %in% scid_groups) %>% 
  unite(id, Gene, MONDO, label, sep = "__") %>% 
  select(id, HPO_ID) %>% 
  distinct() %>% 
  mutate(count = 1) %>% 
  pivot_wider(names_from = "HPO_ID", values_from = "count", values_fill = 0) %>% 
  column_to_rownames("id") %>% 
  vegan::vegdist(method = "jaccard") %>% 
  broom::tidy() %>% 
  separate(item1, into = c("gene1", "mondo1", "proband1"), sep = "__") %>% 
  separate(item2, into = c("gene2", "mondo2", "proband2"), sep = "__") %>% 
  filter(gene1 == gene2 & mondo1 == mondo2) %>% 
  unite(label, gene1, mondo1, sep = "_") %>% 
  distinct() %>% 
  left_join(df_classification_key, by = "label") %>% 
  mutate(distance = 1-distance)

df_scid_plot_count <- df_proband_allClingen %>% 
  separate(label, into = c("gene", "mondo", "proband"), sep = "__") %>% 
  count(gene, mondo) %>% 
  filter(gene %in% table1_genes) %>% 
  unite(label, gene, mondo, sep = "_") %>% 
  filter(label %in% unique(df_scid_plot$label)) %>% 
  left_join(
    df_gene %>% 
      unite(label, Gene, MONDO) %>% 
      select(label, classification),
    by = c("label")
  )

###############################################################################!
# Plots ####
###############################################################################!
g_hpo <- graph_from_data_frame(df_hpo_adjacency, directed = FALSE)

custom_theme <- list(
  theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 7)
  )
)


plt_consistency <- df_scid_plot %>% 
  mutate(metric = "GIC") %>% 
  bind_rows(mutate(df_scid_plot_jaccard, metric = "Jaccard")) %>% 
  mutate(metric = factor(metric, levels = c("GIC", "Jaccard"))) %>% 
  mutate(classification = factor(classification, levels = c("Definitive", "Strong", "Moderate", "Limited"))) %>% 
  ggplot() +
  geom_hline(yintercept = seq(1.5, 100, 1), color = "gray50", linewidth = 0.3, linetype = "dotted") +  # Adjust upper limit as needed
  geom_boxplot(
    aes(y = forcats::fct_reorder(label, distance, .fun = median, .desc = F), x = distance, fill = metric), 
    position = position_dodge(width = 0.7),
    outliers = F, width = 0.6, alpha = 0.5) +
  geom_label(data = mutate(df_scid_plot_count, classification = factor(classification, levels = c("Definitive", "Strong", "Moderate", "Limited"))), 
             aes(y = label, label = n, x = 1.05), 
             size = 3, label.padding = unit(0.2, "lines"), label.size = NA)+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major.y = element_blank(),  # Turn off default gridlines
    panel.grid.minor.y = element_blank(),
    panel.spacing.y = unit(0.5, "lines"),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  labs(y = NULL, x = "Phenotypic consistency", fill = "Metric") +
  facet_grid(classification ~ metric, scales = "free", space = "free")+
  ggforce::facet_col(facets = vars(classification), scales = "free_y", space = "free") +
  scale_fill_brewer(palette = "Dark2") +
  guides(fill = guide_legend(reverse = TRUE))+
  custom_theme



n <- 5
plt_cd3d <- geneHPOgraph(g_hpo, gene = "CD3D", n_node_labels = n, include_title = F, wrap_width = 30 )
plt_cd3e <- geneHPOgraph(g_hpo, gene = "CD3E", n_node_labels = n, include_title = F )


diarrhea_terms <-
  c(
    "HP:0002014",
    "HP:0002028",
    "HP:0025086",
    "HP:0002041",
    "HP:0005208",
    "HP:0004385",
    "HP:0002254",
    "HP:0025085",
    "HP:0033343"
  )
plt_diarrhea <- termHPOgraph(g_hpo, diarrhea_terms)


###############################################################################!
# Assemble ####
###############################################################################!

guide_format <- list(scale_fill_viridis_c(
  guide = guide_colorbar(
    barwidth = 10, 
    barheight = 0.5
  )
))

margin_format <- list(theme(plot.margin = margin(5,25,0,2, "pt")))

plt_consistency_format <- plt_consistency
plt_diarrhea_format <- plt_diarrhea 
plt_cd3d_format <- plt_cd3d +  guide_format
plt_cd3e_format <- plt_cd3e + guide_format

plt_assembled <- plot_grid(
  plt_consistency_format, plot_grid(
    plt_diarrhea_format + theme(plot.margin = margin(20,25,20,2, "pt")), 
    plt_cd3d_format + theme(plot.margin = margin(0,25,0,2, "pt")) , 
    plt_cd3e_format + theme(plot.margin = margin(0,25,0,2, "pt")), 
    ncol = 1, 
    labels = c('B','C','D')
  ),
  nrow = 1, labels = 'A', rel_widths = c(5,6)
)



ggsave(here("figures/4_hpo_consistency.pdf"),
       plt_assembled,
       width = 8, height = 9,
       device = cairo_pdf)
