####
#### R script for Ohigashi et al (2025)
#### plot relationships among IS, SST change, and shift
#### 2025.06.17 written by Ohigashi
#### R 4.5.0
####


### load packages
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
source("Function/F2_HelperFunction_vizualize.R")


### load data
# data frame
eo_SST_IS_shift <- readRDS("08_correlations_out/SSTchange_shift_IS_EAsiaOceania.rds")

# fishcount data (to get number of observed times)
ps_all <- readRDS("01_dataformatting_out/fishdata_ps.rds")

# fish ecological/utilization information
tax_sheet <- read.csv("01_dataformatting_out/fishinfo_w_ecology.csv")


### preparation
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
eo_SST_IS_shift <- eo_SST_IS_shift %>%
  left_join(tax_sheet %>% select(Fish_ID, FishBase_name), by = "FishBase_name")
eo_SST_IS_shift <- eo_SST_IS_shift %>%
  left_join(fish_obs_count, by = "Fish_ID")


# create a model
model <- lm(slope.w ~ SST_slope_NS_ave * mean_IS, data = eo_SST_IS_shift)
model_coef <- coef(model)


### create a plot with continuously predicted area by the model
# create a data frame
pred_df <- expand.grid(SST_slope_NS_ave = seq(min(eo_SST_IS_shift$SST_slope_NS_ave), max(eo_SST_IS_shift$SST_slope_NS_ave), length.out = 100),
                       mean_IS = seq(min(eo_SST_IS_shift$mean_IS), max(eo_SST_IS_shift$mean_IS), 0.001))

# prediction by the model
pred_df$predicted_slope.w <- predict(model, newdata = pred_df, type = "response")

# plot
p_cont <- ggplot(pred_df, aes(SST_slope_NS_ave, predicted_slope.w*111, color = mean_IS,
                    group = mean_IS)) +
  geom_line(alpha = 0.4) +
  geom_point(data = eo_SST_IS_shift,
             aes(SST_slope_NS_ave, slope.w*111, fill = mean_IS), size = 3, shape = 21, color = "white", show.legend = F) +
  labs(x = expression("SST change rate at species’ mean occurrence latitude (\U00B0"*"C year"^{'\U2212'*"1"}*")"),
       y = expression("Estimated poleward shift velocity (km year"^{'\U2212'*"1"}*")"),
       color = "Mean dynamic\nresponse\nto temperature") +
  theme_classic() +
  scale_color_viridis_c(labels = scaleFUN2) +
  scale_fill_viridis_c() +
  scale_y_continuous(labels = scaleFUN_int) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  coord_cartesian(ylim = range((eo_SST_IS_shift$slope.w)*111), xlim = range(eo_SST_IS_shift$SST_slope_NS_ave)) +
  theme(legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 10, color = "black"),
        axis.text=element_text(size=14, color = "black"),
        axis.title=element_text(size=13, color = "black"))


### create a plot with categorically fitted lines
# set break points for mean IS
breaks <- quantile(eo_SST_IS_shift$mean_IS, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

# set labels
labels_custom <- c(
  sprintf("Negative (< %.2f)", breaks[2]),
  sprintf("Near neutral (%.2f to %.2f)", breaks[2], breaks[3]),
  sprintf("Positive (> %.2f)", breaks[3])
)
labels_custom <- gsub("-", "\U2212", labels_custom)

# 
eo_SST_IS_shift <- eo_SST_IS_shift %>%
  mutate(IS_category_tertile = cut(mean_IS,
                                   breaks = breaks,
                                   include.lowest = TRUE,
                                   labels = labels_custom))

viridis_colors <- viridis::viridis(3)
names(viridis_colors) <- levels(eo_SST_IS_shift$IS_category_tertile)

# plot
p_cate <- ggplot(data = eo_SST_IS_shift,
                 aes(SST_slope_NS_ave, slope.w*111, fill = mean_IS)) +
  geom_point(size = 3, shape = 21, color = "black") +
  labs(x = expression("SST change rate at species’ mean occurrence latitude (\U00B0"*"C year"^{'\U2212'*"1"}*")"),
       y = expression(" Estimated poleward shift velocity (km year"^{'\U2212'*"1"}*")"),
       fill = "Mean dynamic\nresponse\nto temperature",
       color = "Response category") +
  geom_smooth(method = "lm", aes(group = IS_category_tertile, color = IS_category_tertile), se = FALSE) +
  scale_fill_viridis_c(labels = scaleFUN2) +
  scale_color_manual(values = viridis_colors) +
  scale_y_continuous(labels = scaleFUN_int) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme(legend.title = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 10, color = "black"),
        axis.text=element_text(size=14, color = "black"),
        axis.title=element_text(size=13, color = "black"))
p_cate


### save data
# continuous area
saveRDS(p_cont, "FigCode/Fig_correlations_out/SST_shift_IS_continuous.rds")
ggsave("FigCode/Fig_correlations_out/SST_shift_IS_continuous.png", plot = p_cont,
       width = 8, height = 5, bg = "white")

# categorical lines
saveRDS(p_cate, "FigCode/Fig_correlations_out/SST_shift_IS_categorical.rds")
ggsave("FigCode/Fig_correlations_out/SST_shift_IS_categorical.png", plot = p_cate,
       width = 8, height = 5, bg = "white")
