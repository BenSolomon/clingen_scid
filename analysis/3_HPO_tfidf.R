library(here)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(rhdf5)
library(tidytext)
library(ggplot2)
library(patchwork)

###############################################################################!
# Functions ####
###############################################################################!

plotHPOheatmap <- function(df, x_var, filtering_var, filter_selections, n_hpo = 10, x_wrap_width = 40, y_wrap_width = 40){
  # browser()
  df <- df %>% 
    mutate(hpo = stringr::str_wrap(hpo, width = y_wrap_width))
  
  topHPO <- df %>% 
    filter(.data[[filtering_var]] %in% filter_selections) %>% 
    slice_max(order_by = tf_idf, n = n_hpo, by = !!sym(x_var)) %>% 
    arrange(!!sym(x_var), desc(tf_idf)) 
  
  orderHPO <- unique(topHPO$hpo)
  orderVar <- unique(topHPO[[x_var]])
  
  
  df <- df[df$hpo %in% orderHPO, ]
  df <- df[df[[filtering_var]] %in% filter_selections, ]
  df <- df[, c(x_var, "hpo", "tf_idf")]
  df <- complete(df, !!sym(x_var), hpo, fill = list(tf_idf = 0))
  df$hpo <- factor(df$hpo, levels = rev(orderHPO))
  df[[x_var]] <- stringr::str_wrap(df[[x_var]], width = x_wrap_width)
  df[[x_var]] <- forcats::fct_relevel(
    df[[x_var]], 
    stringr::str_wrap(orderVar, width = x_wrap_width))
  
  
  ggplot(df, aes_string(x = x_var, y = "hpo", fill = "tf_idf"))+
    geom_tile()+
    scale_fill_viridis_c() +
    theme_minimal()+
    theme(axis.text.x = element_text(angle= 90, hjust= 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL, fill = "TF-IDF")
  
}

###############################################################################!
# Load ####
###############################################################################!
base_dir <- c("/labs/khatrilab/solomonb/gcep")

proband_path <- file.path(base_dir, "data/clingen/api/scid_hpo_funSimAvg.h5")
df_scid_hpo <- data.frame(h5read(proband_path, "/proband_hpo")) 
df_clingen_hpo <- vroom::vroom(file.path(base_dir, "data/clingen_scrape/clingen_proband_hpo_terms.csv.gz"))
df_key <- read_csv(file.path(base_dir, "data/clingen_key.csv"))

###############################################################################!
# Calculate TF-IDF ####
###############################################################################!
df_cdwg_tfidf <- df_clingen_hpo %>%
  select(-MONDO) %>% 
  left_join(
    distinct(df_key[, c("Gene", "Disease", "CDWG")])
    , by = c("Gene", "Disease")) %>% 
  unite(hpo, HPO_ID, HPO_term, sep = " - ") %>% 
  count(CDWG, hpo) %>% 
  mutate(total = sum(n), .by = CDWG) %>% 
  bind_tf_idf(hpo, CDWG, n) 


df_gcep_tfidf <- df_clingen_hpo %>%
  select(-MONDO) %>% 
  left_join(
    distinct(df_key[, c("Gene", "Disease", "GCEP")])
    , by = c("Gene", "Disease")) %>% 
  unite(hpo, HPO_ID, HPO_term, sep = " - ") %>% 
  count(GCEP, hpo) %>% 
  mutate(total = sum(n), .by = GCEP) %>% 
  bind_tf_idf(hpo, GCEP, n) 

df_gdr_tfidf <- df_clingen_hpo %>%
  select(-MONDO) %>% 
  left_join(
    distinct(df_key[, c("Gene", "Disease")])
    , by = c("Gene", "Disease")) %>% 
  unite(gdr, Gene, Disease, sep = " - ") %>% 
  unite(hpo, HPO_ID, HPO_term, sep = " - ") %>% 
  count(gdr, hpo) %>% 
  mutate(total = sum(n), .by = gdr) %>% 
  bind_tf_idf(hpo, gdr, n) 

###############################################################################!
# Plot ####
###############################################################################!

textsize <- 7
titlesize <- textsize*1.4

custom_theme <- list(
  scale_y_discrete(position = "right"),
  theme(
    axis.text = element_text(size = textsize), 
    legend.text = element_text(size = textsize), 
    legend.title = element_text(size = titlesize, vjust = 1),
    legend.direction = "horizontal",
    legend.key.height = unit(0.2, "cm")
    # legend.justification = c(0,0)
  )
)

remove_hpo <-
  c(
    "HP:0001882 - Leukopenia",
    "HP:0002090 - Pneumonia",
    "HP:0002205 - Recurrent respiratory infections",
    "HP:0005403 - T lymphocytopenia",
    "HP:0031690 - Opportunistic infection",
    "HP:0002718 - Recurrent bacterial infections",
    "HP:0001888 - Lymphopenia",
    "HP:0004430 - Severe combined immunodeficiency"
  )

###############################################################################!
## Individual ####

filter_gcep <- c("SCID-CID", "Primary Immune Regulatory Disorders", "Antibody Deficiencies")
plt_heatmap_gcep <- df_gcep_tfidf %>%
  filter(GCEP %in% filter_gcep) %>%
  # filter(!(hpo %in% remove_hpo)) %>% 
  plotHPOheatmap(filtering_var = "GCEP", filter_selections = filter_gcep, x_var = "GCEP", n_hpo = 15, x_wrap_width = 20, y_wrap_width = 45)+
  custom_theme + theme(legend.position = c(1.6, -0.05))


filter_genes_1 <- c("DCLRE1C", "NHEJ1", "PRKDC","RAG1", "RAG2")
plt_heatmap_gdr_1 <- df_gdr_tfidf %>%
  separate(gdr, into = c("Gene", "Disease"), sep =" - ", remove = F) %>% 
  left_join(df_key, by = c("Gene", "Disease")) %>% 
  filter(GCEP == "SCID-CID") %>% 
  filter(!(hpo %in% remove_hpo)) %>% 
  plotHPOheatmap(filtering_var = "Gene", filter_selections = filter_genes_1, x_var = "Gene", n_hpo = 5, x_wrap_width = 50)+
  custom_theme+ theme(legend.position = c(1.7, -0.1))

filter_genes_2 <- c("B2M", "CD8A", "TAP1", "TAP2", "TAPBP")
plt_heatmap_gdr_2 <- df_gdr_tfidf %>%
  separate(gdr, into = c("Gene", "Disease"), sep =" - ", remove = F) %>% 
  left_join(df_key, by = c("Gene", "Disease")) %>% 
  filter(GCEP == "SCID-CID") %>% 
  filter(!(hpo %in% remove_hpo)) %>% 
  plotHPOheatmap(filtering_var = "Gene", filter_selections = filter_genes_2, x_var = "Gene", n_hpo = 5, y_wrap_width = 45) +
  custom_theme+ theme(legend.position = c(1.7, -0.1))

###############################################################################!
## Combined ####
plt_combined <- (plt_heatmap_gcep | (plt_heatmap_gdr_1   / plt_heatmap_gdr_2 + plot_layout(heights = c(18,25)))  ) + 
  plot_layout(widths = c(1,1)) + 
  plot_annotation(
    tag_levels = 'A'
  ) &
  theme(plot.tag = element_text(size = 18, face = "bold"), plot.caption = ggtext::element_textbox_simple(size = 10, hjust = 0))
plt_combined

pdf(here("figures/3_tf_idf.pdf"), width = 8, height = 9.5)
plt_combined
dev.off()
