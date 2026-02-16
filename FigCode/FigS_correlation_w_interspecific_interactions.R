####
#### R script for Ohigashi et al (2025)
#### plot correlations between variables
#### 2025.11.25 written by Ohigashi
#### R 4.5.2
####


### load packages
library(tidyverse); packageVersion("tidyverse")
library(ggrepel); packageVersion("ggrepel")
library(ggimage); packageVersion("ggimage")
library(cowplot); packageVersion("cowplot")
source("Function/F2_HelperFunction_vizualize.R")


### load data
# data frame including shift velocity
analysis_df <- readRDS("08_correlations_out/analysis_df.rds")

# mean IS accounting for interspecific interactions
IS_w.int <- readRDS("09_interspecific_interaction_out/meanIS_w.sp.interaction_df.rds")

# number of significant windows in moving-window UIC
sig_counts <- readRDS("07_movingwindow_UIC_out/n_sig_mwUIC_df.rds")

# fish ecological/utilization information
tax_sheet <- read.csv("01_dataformatting_out/fishinfo_w_ecology.csv")


### format data
analysis_df_ed <- analysis_df %>%
  select(Fish_ID, FishBase_name, Lat_center, mean_IS, slope.w, mean_abun) %>%
  left_join(IS_w.int %>% select(Fish_ID, FishBase_name, mean_IS_w.int = mean_IS),
            by = c("Fish_ID", "FishBase_name")) %>% 
  left_join(sig_counts, by = c("Fish_ID" = "effect_var"))

analysis_df_ed <- analysis_df_ed %>%
  filter(n_window_sig >= 1)


### calculate correlation between mean IS accounting for interspecific interactions and shift velocity
# correlation test
test_result <- cor.test(analysis_df_ed[["mean_IS_w.int"]], analysis_df_ed[["slope.w"]], use = "complete.obs")
# store the result
r <- test_result$estimate
p <- test_result$p.value


### visualize
r.tmp <- r
p.tmp <- p
# make annotation text
anno.tmp <- paste0("italic(r) == ", scaleFUN2(r.tmp), ' * ", " * ', p_to_sig(p.tmp, type = "threshold"))
anno.tmp <- sub("\U2212", "-", anno.tmp)

# format data frame for plot
plot_df <- analysis_df_ed %>%
  filter(!is.na(Lat_center), !is.na(mean_IS_w.int)) %>%
  left_join(tax_sheet %>% select(Fish_ID, Order, Family), by = "Fish_ID")

# extract important species
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

## plot
p_IS_shift <- ggplot(plot_df %>% filter(!is.na(slope.w)),
                     aes(x = mean_IS_w.int, y = slope.w*111, color = Order
                     )) +
  geom_point(alpha = 0.7) +
  theme_classic() +
  theme(legend.title = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        axis.text=element_text(size=14, color = "black"),
        axis.title=element_text(size=13, color = "black")) +
  geom_smooth(method = "lm", color = "black", show.legend = F) +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  scale_y_continuous(labels = scaleFUN_int) +
  scale_x_continuous(labels = scaleFUN2) +
  scale_color_manual(values = c("Clupeiformes" = "#C49A00", "Perciformes" = "#00C094", "Tetradontiformes" = "#FB61D7"))+
  labs(x = "Mean dynamic response to temperature\n(excluding interspecific interactions)",
       y = expression("Estimated poleward shift velocity (km year"^{'\U2212'*"1"}*")")) +
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
  mutate(mean_IS_w.int_ed = ifelse(FishBase_name %in% c("Acanthopagrus_schlegelii", "Parajulis_poecilepterus"),
                             mean_IS_w.int - 0.015,
                             mean_IS_w.int + 0.015))


p_IS_shift_shape <- p_IS_shift +
  geom_image(data = important_shift_df,
             aes(x = mean_IS_w.int_ed, y = slope.w*111, color = Order,
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


### interaction effect of SST change and mean IS (accounting for interspecific interactions)
# load data frame
eo_SST_IS_shift <- readRDS("08_correlations_out/SSTchange_shift_IS_EAsiaOceania_all.sp.rds")

# add IS accounting for interspecific interaction
IS_w.int <- readRDS("09_interspecific_interaction_out/meanIS_w.sp.interaction_df.rds")

eo_SST_IS_shift_ed <- eo_SST_IS_shift %>%
  mutate(SST_slope_NS_ave = rowMeans(
    dplyr::select(., SST_slope_Northern, SST_slope_Southern),
    na.rm = TRUE
  ))

eo_SST_IS_shift_ed <- eo_SST_IS_shift_ed %>%
  left_join(IS_w.int %>% select(FishBase_name, mean_IS_w.int = mean_IS), by = "FishBase_name") %>%
  filter(!is.na(mean_IS_w.int) & !is.na(slope.w))

eo_SST_IS_shift_ed <- eo_SST_IS_shift_ed %>%
  select(-n_window_sig) %>% # remove old column
  left_join(sig_counts, by = c("Fish_ID" = "effect_var")) %>%
  filter(n_window_sig >= 1)

# create a model
model <- lm(slope.w ~ SST_slope_NS_ave * mean_IS_w.int, data = eo_SST_IS_shift_ed)
model_coef <- coef(model)

### create a plot with continuously predicted area by the model
# create a data frame
pred_df <- expand.grid(SST_slope_NS_ave = seq(min(eo_SST_IS_shift_ed$SST_slope_NS_ave), max(eo_SST_IS_shift_ed$SST_slope_NS_ave), length.out = 100),
                       mean_IS_w.int = seq(min(eo_SST_IS_shift_ed$mean_IS_w.int), max(eo_SST_IS_shift_ed$mean_IS_w.int), 0.001))

# prediction by the model
pred_df$predicted_slope.w <- predict(model, newdata = pred_df, type = "response")

# plot
p_cont <- ggplot(pred_df, aes(SST_slope_NS_ave, predicted_slope.w*111, color = mean_IS_w.int,
                              group = mean_IS_w.int)) +
  geom_line(alpha = 0.4) +
  geom_point(data = eo_SST_IS_shift_ed,
             aes(SST_slope_NS_ave, slope.w*111, fill = mean_IS_w.int), size = 3, shape = 21, color = "white", show.legend = F) +
  labs(x = expression("SST change rate (\U00B0"*"C year"^{'\U2212'*"1"}*")"),
       y = expression("Estimated poleward shift velocity (km year"^{'\U2212'*"1"}*")"),
       color = "Mean dynamic response\nto temperature\n(excluding interspecific interactions)") +
  theme_classic() +
  scale_color_viridis_c(labels = scaleFUN2) +
  scale_fill_viridis_c() +
  scale_y_continuous(labels = scaleFUN) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  coord_cartesian(ylim = range((eo_SST_IS_shift_ed$slope.w)*111), xlim = range(eo_SST_IS_shift_ed$SST_slope_NS_ave)) +
  theme(legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 10, color = "black"),
        axis.text=element_text(size=14, color = "black"),
        axis.title=element_text(size=13, color = "black"))

### create a plot with categorically fitted lines
# set break points for mean IS
breaks <- c(min(eo_SST_IS_shift_ed$mean_IS_w.int), 0, max(eo_SST_IS_shift_ed$mean_IS_w.int))

# set labels
labels_custom <- c(
  sprintf("Negative (< %.0f, N = 5)", breaks[2]),
  sprintf("Positive (> %.0f, N = 8)", breaks[2])
)
labels_custom <- gsub("-", "\U2212", labels_custom)

# 
eo_SST_IS_shift_ed <- eo_SST_IS_shift_ed %>%
  mutate(IS_category_tertile = cut(mean_IS_w.int,
                                   breaks = breaks,
                                   include.lowest = TRUE,
                                   labels = labels_custom))

viridis_colors <- viridis::viridis(2)
names(viridis_colors) <- levels(eo_SST_IS_shift_ed$IS_category_tertile)

# plot
p_cate <- ggplot(data = eo_SST_IS_shift_ed,
                 aes(SST_slope_NS_ave, slope.w*111, fill = mean_IS_w.int)) +
  geom_point(size = 3, shape = 21, color = "black") +
  labs(x = expression("SST change rate (\U00B0"*"C year"^{'\U2212'*"1"}*")"),
       y = expression("Estimated poleward shift velocity (km year"^{'\U2212'*"1"}*")"),
       fill ="Mean dynamic response\nto temperature\n(excluding interspecific interactions)",
       color = "Response category") +
  geom_smooth(method = "lm", aes(group = IS_category_tertile, color = IS_category_tertile), se = FALSE) +
  scale_fill_viridis_c(labels = scaleFUN2) +
  scale_color_manual(values = viridis_colors) +
  scale_y_continuous(labels = scaleFUN) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme(legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 10, color = "black"),
        axis.text=element_text(size=14, color = "black"),
        axis.title=element_text(size=12, color = "black"))


### save result
dir.create("FigCode/FigS_correlation_w_interspecific_interactions_out")

# IS (accounting for interspecific interaction) and shift velocity with fish shape
saveRDS(p_IS_shift_shape_fin, file = "FigCode/FigS_correlation_w_interspecific_interactions_out/cor_mean_IS_shift_fishapes_w.sp.int.rds")
ggsave("FigCode/FigS_correlation_w_interspecific_interactions_out/cor_mean_IS_shift_fishapes_w.sp.int.png", p_IS_shift_shape_fin,
       width = 10, height = 6, bg = "white")

## SST change x IS on shift velocity
# continuous
ggsave("FigCode/FigS_correlation_w_interspecific_interactions_out/SST_shift_IS_continuous_w.sp.int.png", plot = p_cont,
       width = 8, height = 5, bg = "white")
# categorical
ggsave("FigCode/FigS_correlation_w_interspecific_interactions_out/SST_shift_IS_categorical_w.sp.int.png", plot = p_cate,
       width = 8, height = 5, bg = "white")
