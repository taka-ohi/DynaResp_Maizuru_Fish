####
#### R script for Ohigashi et al (2025)
#### plot correlations between variables
#### 2025.06.17 written by Ohigashi
#### R 4.5.0
####


### load packages
library(tidyverse); packageVersion("tidyverse")
library(ggrepel); packageVersion("ggrepel")
library(ggimage); packageVersion("ggimage")
library(cowplot); packageVersion("cowplot")
library(phyloseq); packageVersion("phyloseq")
library(tidyr); packageVersion("tidyr")
source("Function/F2_HelperFunction_vizualize.R")


### load data
# mean IS, latitude, and so on
analysis_df <- readRDS("08_correlations_out/analysis_df.rds")

# fish ecological/utilization information
tax_sheet <- read.csv("01_dataformatting_out/fishinfo_w_ecology.csv")

# correlation result
cor_normal <- read.csv("08_correlations_out/correlation_normal.csv")

# fishcount data (to get number of observed times)
ps_all <- readRDS("01_dataformatting_out/fishdata_ps.rds")


### 0. preparation
# get fishcount table
fishcount <- otu_table(ps_all) %>% as.data.frame()
fishcount <- fishcount %>% rownames_to_column("Sample")

# calculate observed times
fish_obs_count <- fishcount %>%
  summarise(across(
    starts_with("Fish"),
    ~ sum(. > 0, na.rm = TRUE)
  )) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Fish_ID",
    values_to = "n_observed"
  )

# combine with the other data 
analysis_df <- analysis_df %>%
  left_join(fish_obs_count, by = "Fish_ID")

### 1. visualize correlation between distribution center and IS
## latitude x mean IS correlation
# extract correlation result
r.tmp <- cor_normal %>%
  filter(v1 == "Lat_center", v2 == "mean_IS") %>% pull(r)
p.tmp <- cor_normal %>%
  filter(v1 == "Lat_center", v2 == "mean_IS") %>% pull(p)
# make annotation text
anno.tmp <- paste0("italic(r) == ", scaleFUN2(r.tmp), ' * ", " * ', p_to_sig(p.tmp, type = "threshold"))
anno.tmp <- sub("\U2212", "-", anno.tmp)

## extract important species to show in the plot
# cutoff of top 5 species
top5_cutoff <- analysis_df %>%
  arrange(desc(mean_abun)) %>%
  slice_head(n = 5) %>%
  pull(mean_abun) %>%
  min()
# extract importance information from fishbase
important_df <- analysis_df %>%
  left_join(tax_sheet %>% select(Fish_ID, Importance, Order, Family), by = "Fish_ID") %>%
  filter(!is.na(Lat_center), !is.na(mean_IS), n_window_sig >= 1) %>%
  filter(mean_abun >= top5_cutoff | Importance == "highly commercial")

# add italic label for plotting
important_df <- important_df %>%
  mutate(
    Genus = str_extract(FishBase_name, "^[^_]+"),
    Species = str_extract(FishBase_name, "(?<=_).*"),
    genus_initial = paste0(substr(Genus, 1, 1), "."),
    label_expr = paste0("italic('", genus_initial, " ", Species, "')")
  )

# format data frame for plot
plot_df <- analysis_df %>%
  filter(!is.na(Lat_center), !is.na(mean_IS), n_window_sig >= 1) %>%
  left_join(tax_sheet %>% select(Fish_ID, Order, Family), by = "Fish_ID")

## plot
p_lat_IS <- ggplot(plot_df,
       aes(x = Lat_center, y = mean_IS, size = n_observed, color = Order)) +
  geom_point(alpha = 0.7) +
  theme_classic() +
  theme(legend.title = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        axis.text=element_text(size=14, color = "black"),
        axis.title=element_text(size=14, color = "black")) +
  scale_size_continuous(breaks = c(10, 100, 250)) +
  geom_smooth(method = "lm", color = "black", show.legend = F) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  scale_y_continuous(labels = scaleFUN) +
  scale_x_continuous(labels = scaleFUN_int) +
  labs(x = sprintf("Latitudinal center of distribution (%sN)", "\U00B0"),
       y = "Mean dynamic response to temperature",
       size = "Number of obsevations\nin Maizuru Bay") +
  geom_label_repel(data = important_df,
                  aes(label = label_expr, color = Order),
                  parse = TRUE,
                  size = 3,
                  fill = alpha("white", 0.5),
                  # color = "red",
                  segment.alpha = 1, show.legend = F) +
  annotate("text",
           x = Inf, y = Inf, hjust = "right", vjust = 1,
           label = anno.tmp, size = 5,
           parse = TRUE
           )

## plot with fish shapes
important_df <- important_df %>%
  mutate(
    image_path = paste0("Data/Fish_images/", FishBase_name, ".svg")
  )
# edit x coordinates for fish shapes
important_df <- important_df %>%
  mutate(Lat_center_ed = ifelse(FishBase_name == "Sebastes_cheni",
                                Lat_center -1,
                                Lat_center + 1))

p_lat_IS_shape <- ggplot(plot_df,
                   aes(x = Lat_center, y = mean_IS, size = n_observed, color = Order)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_smooth(method = "lm", color = "black", show.legend = F) +
  geom_point(alpha = 0.4) +
  theme_classic() +
  theme(legend.title = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        axis.text=element_text(size=14, color = "black"),
        axis.title=element_text(size=14, color = "black")) +
  scale_size_continuous(breaks = c(10, 100, 250)) +
  scale_y_continuous(labels = scaleFUN) +
  scale_x_continuous(labels = scaleFUN_int) +
  labs(x = sprintf("Latitudinal center of distribution (%sN)", "\U00B0"),
       y = "Mean dynamic response to temperature",
       size = "Number of observations\nin Maizuru Bay") +
  geom_image(data = important_df,
             aes(x = Lat_center_ed, y = mean_IS, image = image_path, color = Order),
             size = 0.06,
             asp = 1.2
             ) + 
  geom_label_repel(data = important_df,
                   aes(label = label_expr, color = Order),
                   parse = TRUE,
                   size = 3,
                   fill = alpha("white", 0.5),
                   # color = "red",
                   segment.alpha = 1, show.legend = F) +
  annotate("text",
           x = Inf, y = Inf, hjust = "right", vjust = 1,
           label = anno.tmp, size = 5,
           parse = TRUE
  )

p_lat_IS_shape

# it is a little bit tricky, but I get legend from the plot w/o fish shape, and join it to the plot w/ fish shape.
lat_IS_leg <- get_legend(p_lat_IS)

p_lat_IS_shape_fin <- plot_grid(p_lat_IS_shape + theme(legend.position = "none"),
                                lat_IS_leg,
                                nrow = 1, rel_widths = c(1, 0.3))

## save plot data
dir.create("FigCode/Fig_correlations_out")
# basic
saveRDS(p_lat_IS, file = "FigCode/Fig_correlations_out/cor_Lat_center_mean_IS.rds")
ggsave("FigCode/Fig_correlations_out/cor_Lat_center_mean_IS.png", plot = p_lat_IS,
       width = 10, height = 7, bg = "white")
# with fish shapes
saveRDS(p_lat_IS_shape_fin, file = "FigCode/Fig_correlations_out/cor_Lat_center_mean_IS_fishapes.rds")
ggsave("FigCode/Fig_correlations_out/cor_Lat_center_mean_IS_fishapes.png", plot = p_lat_IS_shape_fin,
       width = 10, height = 6, bg = "white")


### 2. visualize correlation between IS and range shift velocity
# extract correlation result
r.tmp <- cor_normal %>%
  filter(v1 == "mean_IS", v2 == "slope.w") %>% pull(r)
p.tmp <- cor_normal %>%
  filter(v1 == "mean_IS", v2 == "slope.w") %>% pull(p)
# make annotation text
anno.tmp <- paste0("italic(r) == ", scaleFUN2(r.tmp), ' * ", " * ', p_to_sig(p.tmp, type = "threshold"))
anno.tmp <- sub("\U2212", "-", anno.tmp)

## extract important species to show in the plot
important_shift_df <- plot_df %>%
  filter(!is.na(slope.w)) %>% 
  left_join(tax_sheet %>% select(Fish_ID, Importance), by = "Fish_ID") %>%
  filter(Importance %in% c("highly commercial", "commercial"))

# add italic label for plotting
important_shift_df <- important_shift_df %>%
  mutate(
    Genus = str_extract(FishBase_name, "^[^_]+"),
    Species = str_extract(FishBase_name, "(?<=_).*"),
    genus_initial = paste0(substr(Genus, 1, 1), "."),
    label_expr = paste0("italic('", genus_initial, " ", Species, "')")
  )

scales::show_col(scales::hue_pal()(7))
# Clupeiformes: #CD9600, Perciformes: #00BFC4, Tetradontiformes: #FF61CC

## plot
p_IS_shift <- ggplot(plot_df %>% filter(!is.na(slope.w)),
                     aes(x = mean_IS, y = slope.w*111, color = Order,
                         size = n_observed
                         )) +
  geom_point(alpha = 0.7) +
  theme_classic() +
  theme(legend.title = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        axis.text=element_text(size=14, color = "black"),
        axis.title.x=element_text(size=13, color = "black"),
        axis.title.y=element_text(size=13, color = "black")) +
  geom_smooth(method = "lm", color = "black", show.legend = F) +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  scale_size_continuous(breaks = c(10, 100, 250)) +
  scale_y_continuous(labels = scaleFUN_int) +
  scale_x_continuous(labels = scaleFUN2) +
  scale_color_manual(values = c("Clupeiformes" = "#C49A00", "Perciformes" = "#00C094", "Tetradontiformes" = "#FB61D7"))+
  labs(x = "Mean dynamic response to temperature",
       y = expression("Estimated poleward shift velocity (km year"^{'\U2212'*"1"}*")"),
       size = "Number of observations\nin Maizuru Bay"
       ) +
  geom_label_repel(data = important_shift_df,
                   aes(label = label_expr, color = Order),
                   parse = TRUE,
                   size = 3,
                   fill = alpha("white", 0.5),
                   # color = "red",
                   segment.alpha = 1, show.legend = F) +
  annotate("text",
           x = Inf, y = Inf, hjust = "right", vjust = 1,
           label = anno.tmp, size = 5,
           parse = TRUE
  )
p_IS_shift


## plot with fish shapes
important_shift_df <- important_shift_df %>%
  mutate(
    image_path = paste0("Data/Fish_images/", FishBase_name, ".svg")
  ) 

# edit x coordinates for fish shapes
important_shift_df <- important_shift_df %>%
  mutate(mean_IS_ed = ifelse(FishBase_name %in% c("Acanthopagrus_schlegelii", "Parajulis_poecilepterus"),
                             mean_IS - 0.015,
                             mean_IS + 0.015))


p_IS_shift_shape <- p_IS_shift +
  geom_image(data = important_shift_df,
             aes(x = mean_IS_ed, y = slope.w*111, color = Order,
                 image = image_path),
             size = 0.06,
             # nudge_x = -0.01,
             asp = 1.2, inherit.aes = FALSE
  )
p_IS_shift_shape


# it is a little bit tricky, but I get legend from the plot w/o fish shape, and join it to the plot w/ fish shape.
IS_shift_leg <- get_legend(p_IS_shift)

p_IS_shift_shape_fin <- plot_grid(p_IS_shift_shape + theme(legend.position = "none"),
                                IS_shift_leg,
                                nrow = 1, rel_widths = c(1, 0.3))


## save plot data
# basic
saveRDS(p_IS_shift, file = "FigCode/Fig_correlations_out/cor_mean_IS_shift.rds")
ggsave("FigCode/Fig_correlations_out/cor_mean_IS_shift.png", plot = p_IS_shift,
       width = 8, height = 7, bg = "white")
# with fish shape
saveRDS(p_IS_shift_shape_fin, file = "FigCode/Fig_correlations_out/cor_mean_IS_shift_fishapes.rds")
ggsave("FigCode/Fig_correlations_out/cor_mean_IS_shift_fishapes.png", p_IS_shift_shape_fin,
       width = 10, height = 6, bg = "white")


### 3. visualize other correlations between Lat_center and lm and gam-based coefficients
# create plot list
cor_normal_plotlist <- list()

# make a list of variable pairs
var_pairs <- list(
  # c("Lat_center", "mean_IS"),
  c("Lat_center", "lm_coef"),
  c("Lat_center", "gam_slope_ave"),
  # c("mean_IS", "slope.w"),
  c("lm_coef", "slope.w"),
  c("gam_slope_ave", "slope.w")
)

# create plots
for (pair in var_pairs) {
  v1 <- pair[1]
  v2 <- pair[2]
  
  # remove undetected species by the moving-window UIC in the case of mean IS
  if (v1 == "mean_IS" | v2 == "mean_IS") {
    sub_df0 <- analysis_df %>%
      filter(n_window_sig >= 1) # significant at least one time
  } else {
    sub_df0 <- analysis_df
  }
  
  # select focal data and omit NA
  sub_df <- sub_df0 %>%
    select(FishBase_name, all_of(v1), all_of(v2)) %>%
    na.omit()
  
  # call correlation result
  r.tmp <- cor_normal %>%
    filter(v1 == !!v1, v2 == !!v2) %>% pull(r)
  p.tmp <- cor_normal %>%
    filter(v1 == !!v1, v2 == !!v2) %>% pull(p)
  anno.tmp <- paste0("italic(r) == ", scaleFUN2(r.tmp), ' * ", " * ', p_to_sig(p.tmp, type = "threshold"))
  anno.tmp <- sub("\U2212", "-", anno.tmp)
  
  # set axis titles
  if (v1 == "Lat_center") {
    xtitle <- sprintf("Latitudinal center of distribution (%sN)", "\U00B0")
  } else if (v1 == "lm_coef") {
    xtitle <- "LM-based response to temperature"
  } else if (v1 == "gam_slope_ave") {
    xtitle <- "GAM-based response to temperature"
  } else {
    xtitle <- v1
  }
  if (v2 == "slope.w") {
    ytitle <- expression("Estimated poleward shift velocity (km year"^{'\U2212'*"1"}*")")
  } else if (v2 == "lm_coef") {
    ytitle <- "LM-based response to temperature"
  } else if (v2 == "gam_slope_ave") {
    ytitle <- "GAM-based response to temperature"
  } else {
    ytitle <- v2
  }
  
  sub_df2 <- sub_df %>%
    mutate(y_plot = if (v2 == "slope.w") .data[[v2]] * 111 else .data[[v2]])
  
  # plot
  plo.tmp <- ggplot(sub_df2,
                    aes(x = !!sym(v1), y = y_plot
                    )) +
    geom_point(alpha = 0.7) +
    theme_classic() +
    theme(legend.title = element_text(size = 12, color = "black"),
          legend.text = element_text(size = 12, color = "black"),
          axis.text=element_text(size=13, color = "black"),
          axis.title=element_text(size=13, color = "black")) +
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
    scale_y_continuous(labels = scaleFUN) +
    scale_x_continuous(labels = scaleFUN2) +
    labs(x = xtitle,
         y = ytitle) +
    annotate("text",
             x = Inf, y = Inf, hjust = "right", vjust = 1,
             label = anno.tmp, size = 5,
             parse = TRUE
    )
  
  plo.tmp.fin <- plo.tmp
  if (v1 %in% c("lm_coef", "gam_slope_ave")) {
    plo.tmp.fin <- plo.tmp.fin +
      geom_vline(xintercept = 0, color = "grey", linetype = "dashed")
  }
  if (p.tmp < 0.05) {
    plo.tmp.fin <- plo.tmp.fin +
      geom_smooth(method = "lm", se = TRUE, color = "black")
  }
  
  # store to the list
  plotname <- paste(v1, v2, sep = "_")
  cor_normal_plotlist[[plotname]] <- plo.tmp.fin
}


# save plots as RDS
for (name in names(cor_normal_plotlist)) {
  saveRDS(cor_normal_plotlist[[name]], file = sprintf("FigCode/Fig_correlations_out/cor_%s.rds", name))
}



