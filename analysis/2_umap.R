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
library(patchwork)

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

# Plot that facets data, but also includes a background of all data points
# in each facet
facetOverlayPlot <- function(df, facet_var){
  background_data <- df
  background_data[[facet_var]] <- NULL
  
  ggplot(df, aes_string(x = "V1", y = "V2")) +
    ggrastr::geom_point_rast(data = background_data,  color = "grey") +
    ggrastr::geom_point_rast(size = 2, alpha = 1, shape = 21, aes_string(fill = facet_var)) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(alpha = 1, size = 2))) +
    scale_fill_brewer(palette = "Dark2") +
    facet_wrap(sym(facet_var)) +
    labs(x = "UMAP1", y = "UMAP2", color = "CDWG") +
    theme(legend.position = "none")
}

# Plot that takes a list of genes and highlights all probands with genes 
# from that list as one group and all other probands as part of a null group
selectGenePlot <- function(df, genes, label){
  ggplot(df, aes_string(x = "V1", y = "V2")) +
    ggrastr::geom_point_rast(size = 2, color = "grey") +
    ggrastr::geom_point_rast(
      data = filter(df, Gene %in% genes), 
      size = 2, shape = 21, fill = "forestgreen") +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
    labs(x = "UMAP1", y = "UMAP2") +
    ggtitle(label)
  # scale_fill_viridis_d(begin = 0.1, end = 0.95, direction = -1)
}

hpoFeaturePlot <- function(df_umap, hpo_terms){
  hpo_terms <- sprintf("(?i)%s", hpo_terms) # Add case insensitivity
  hpo_pattern <- paste(hpo_terms, collapse = "|")
  df_proband_term_annotation <- df_proband_hpo %>% 
    rename(gene = Gene, proband_id = label) %>%  
    createProbandLabel() %>% 
    mutate(term = stringr::str_extract(
      HPO_term, 
      pattern = hpo_pattern
    )) %>% 
    mutate(term = stringr::str_to_title(term))
  
  df_umap <- df_umap %>% 
    left_join(df_proband_term_annotation, by = "label")
  
  ggplot(drop_na(df_umap), aes_string(x = "V1", y = "V2")) +
    ggrastr::geom_point_rast(data = select(df_umap, -term),  color = "grey") +
    ggrastr::geom_point_rast(size = 2, alpha = 1, shape = 21, fill = "darkorange" ) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(alpha = 1, size = 2))) +
    labs(x = "UMAP1", y = "UMAP2") +
    facet_wrap(~term) +
    theme(legend.position = "none")
}

readUMAP <- function(path, cdwgs){
  df <- readRDS(path)
  df <- filter(df, CDWG %in% cdwgs)
  return(df)
}

###############################################################################!
# Load data ####
###############################################################################!

keep_cdwg <- c("Cardiovascular", "Hemostasis/Thrombosis", "Hereditary Cancer",
               "Immunology", "Inborn Errors of Metabolism",
               "Neurodevelopmental Disorders", "Neurological Disorders",
               "Pulmonary")


h5_path_allClingen <- here("data/clingen/scrape/clingen_scrape_hpo.h5")
mtx_proband_allClingen <- h5read(h5_path_allClingen, "/proband_distance")
df_proband_allClingen <- data.frame(h5read(h5_path_allClingen, "/proband_metadata")) %>% 
  unite(label, gene, disease, proband_id, sep = "__", remove = F)

df_proband_allClingen <- createProbandLabel(df_proband_allClingen)

rownames(mtx_proband_allClingen) <- df_proband_allClingen$label
colnames(mtx_proband_allClingen) <- df_proband_allClingen$label

df_key <- read_csv(here("data/clingen_key.csv")) %>% 
  dplyr::select(Gene, Disease, MONDO, GCEP, CDWG)

df_proband_hpo <- vroom::vroom(here("data/clingen/scrape/clingen_proband_hpo_terms.csv.gz"))

allClingen_path <- here("data/clingen_scape_umap_n400_d0-6.RDS")
df_umap_allClingen <- readUMAP(allClingen_path, keep_cdwg)
immunoCDWG_path <- here("data/clingen_immunoCDWG_umap_n50_d0-2.RDS")
df_umap_immunoCDWG <- readRDS(immunoCDWG_path)
scidGCEP_path <- here("data/clingen_scidGCEP_umap_n25_d0-2.RDS")
df_umap_scidGCEP <- readRDS(scidGCEP_path)


iuis_files <- list.files(here("data/gene_tables"), 
                         pattern = "IUIS_1", 
                         full.names = TRUE)
table1_genes <- lapply(iuis_files, read_lines)
table1_genes <- unlist(table1_genes)


###############################################################################!
# Plots ####
###############################################################################!


plt_all <- df_umap_allClingen %>% 
  mutate(CDWG = stringr::str_wrap(CDWG, width = 30)) %>% 
  ggplot( aes_string(x = "V1", y = "V2", fill = "CDWG")) +
  ggrastr::geom_point_rast(size = 2, alpha = 1, shape = 21, color = "grey10", stroke = 0.2) +
  theme_bw() +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "All-UMAP1", y = "All-UMAP2", fill = NULL) 

plt_all_iuis <- selectGenePlot(df_umap_allClingen, table1_genes, NULL)+
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "All-UMAP1", y = "All-UMAP2", color = "CDWG") 

plt_immuno <- df_umap_immunoCDWG %>% 
  mutate(GCEP = stringr::str_wrap(GCEP, width = 30)) %>% 
  ggplot(aes_string(x = "V1", y = "V2", fill = "GCEP")) +
  ggrastr::geom_point_rast(size = 2, alpha = 1, shape = 21, stroke = 0.3) +
  theme_bw() +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(x = "Immunology-UMAP1", y = "Immunology-UMAP2", fill = NULL) +
  scale_fill_brewer(palette = "Dark2") 

plt_immuno_iuis <- selectGenePlot(df_umap_immunoCDWG, table1_genes, NULL) +
  labs(x = "Immunology-UMAP1", y = "Immunology-UMAP2")

df_hpo_tfidf <- df_proband_hpo %>% 
  left_join(df_key, by = c("Gene", "MONDO"), relationship = "many-to-many") %>%
  rename(gene = Gene, proband_id = label) %>% 
  createProbandLabel() %>% 
  unite(hpo, HPO_ID, HPO_term, sep = "___") %>% 
  count(GCEP, hpo, sort = T) %>% 
  mutate(total = sum(n), .by = GCEP) %>% 
  bind_tf_idf(hpo, GCEP, n)

plt_scid_features <-
  hpoFeaturePlot(
    df_umap_scidGCEP,
    hpo_terms = c(
      "failure to thrive",
      "pneumocystis",
      "chronic mucocutaneous candidiasis"
    )
  ) + labs(x = "SCID-UMAP1", y = "SCID-UMAP2")

###############################################################################!
# Assemble ####
###############################################################################!
plt_assembled <- (plt_all + plt_all_iuis + plt_immuno + plt_immuno_iuis + plot_layout(ncol = 2)) /
  plt_scid_features +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    tag_levels = 'A',
  ) &
  theme(plot.tag = element_text(size = 18, face = "bold"), plot.caption = ggtext::element_textbox_simple(size = 10, hjust = 0))

pdf(here("figures/2_umap.pdf"), width = 8, height = 9)
plt_assembled
dev.off()