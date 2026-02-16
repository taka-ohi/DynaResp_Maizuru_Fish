####
#### R script for Ohigashi et al (2025)
#### visualize moving-window UIC
#### 2025.10.30 written by Ohigashi
#### R 4.5.1
#### 


### load packages and functions
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(tidyr); packageVersion("tidyr")
library(stringr); packageVersion("stringr")
library(ggplot2); packageVersion("ggplot2")


### load data
# moving-window UIC result
final_df <- readRDS("07_movingwindow_UIC_out/movingwindow_UIC_resultdf.rds")

# fish ecological/utilization information
tax_sheet <- read.csv("01_dataformatting_out/fishinfo_w_ecology.csv")

# sample information
sample_sheet <- read.table("01_dataformatting_out/metadata.txt", header = T)


### calculation of p.adj and formatting
# calculate adjusted p-values
final_df_adj <- final_df %>%
  filter(tp <= 0) %>% 
  group_by(window_start, tp) %>%
  mutate(p.adj = p.adjust(surr_p, method = "BH")) %>% # adjusted among fish within same window and tp
  ungroup()

# add global adjusted p-values (i.e., considering inter-window adjustment)
final_df_adj <- final_df_adj %>%
  group_by(tp) %>%
  mutate(p.adj_global = p.adjust(surr_p, method = "BH")) %>% # adjusted among fish x window within same tp
  ungroup()


### create a heatmap
# consider p < 0.05 values as significant at any tp
heatmap_df <- final_df_adj %>%
  group_by(window_start, effect_var) %>%
  summarise(sig = any(p.adj_global < 0.05), .groups = "drop")

## annotate Fish species names
# fish species name
heatmap_df <- heatmap_df %>%
  left_join(tax_sheet %>% select(Fish_ID, FishBase_name),
            by = c("effect_var" = "Fish_ID"))

## cluster the species
mat <- heatmap_df %>%
  mutate(sig_num = as.numeric(sig)) %>%
  select(-sig, -FishBase_name) %>%
  pivot_wider(names_from = window_start, values_from = sig_num, values_fill = 0) %>%
  as.data.frame()

mat <- mat %>% column_to_rownames("effect_var")

# clustering
row_clust <- hclust(dist(mat, method = "euclidean"), method = "ward.D2")
row_order <- rownames(mat)[row_clust$order]

# add the clustered levels to the Fish ID
heatmap_df$effect_var <- factor(heatmap_df$effect_var, levels = row_order)

# labels for species
sp_labs <- heatmap_df %>%
  distinct(effect_var, FishBase_name) %>%
  mutate(
    genus  = str_extract(FishBase_name, "^[A-Za-z]+"),
    species = str_replace(FishBase_name, "^[A-Za-z]+_", ""),
    short  = paste0(substr(genus, 1, 1), ". ", species)
  )
make_ital <- function(x) {
  # x = "O. altipennis"
  parts <- str_split(x, " ", n = 2)[[1]]
  genus_abbr <- parts[1]
  species    <- parts[2]
  # return plotmath (italicized)
  paste0("italic('", genus_abbr, "')~italic('", species, "')")
}

sp_labs <- sp_labs %>%
  mutate(label_expr = sapply(short, make_ital))

# fix Platycephalus_sp._2
sp_labs <- sp_labs %>%
  mutate(label_expr = case_when(
    FishBase_name == "Platycephalus_sp._2" ~ "italic('Platycephalus')~'sp. 2'",
    TRUE ~ label_expr
  ))

# plot heatmap
heat <- ggplot(heatmap_df, aes(x = window_start, y = effect_var, fill = sig)) +
  # geom_tile(color = "white") +
  geom_tile(color = "grey90", size = 0.2) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "red"),
                    name = "Including\nadjusted p < 0.05") +
  scale_x_continuous(breaks = c(1, 73, 145, 217, 289, 361),
                     labels = c("2002-Jan", "2005-Jan", "2008-Jan", "2011-Jan",
                                "2014-Jan", "2017-Jan"),
                     expand = c(0, 0)) +
  scale_y_discrete(
    labels = rlang::set_names(
      lapply(sp_labs$label_expr, function(z) parse(text = z)[[1]]),
      sp_labs$effect_var
    )
  ) +
  theme_minimal() +
  labs(x = "Window start time", y = "Fish species",
       title = "UIC-detected windows for each fish (window size = 120)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())
heat

### save data
dir.create("FigCode/FigS_windowUIC_out")
# save heatmap
ggsave("FigCode/FigS_windowUIC_out/UIC_detected_windows_sp.png", plot = heat,
       width = 9, height = 14, bg = "white")
cairo_pdf("FigCode/FigS_windowUIC_out/UIC_detected_windows_sp.pdf", width = 9, height = 14)
print(heat)
dev.off()