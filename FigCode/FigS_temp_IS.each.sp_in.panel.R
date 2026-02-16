####
#### R script for Ohigashi et al (2026)
#### visualize interaction strength calculated in MDR S-map
#### 2026.01.15 written by Ohigashi
#### R 4.5.2
####


### load packages and functions
library(tidyverse); packageVersion("tidyverse")
library(tidyr); packageVersion("tidyr")
library(ggtext); packageVersion("ggtext")
source("Function/F2_HelperFunction_vizualize.R")


### load data
# MDR S-map result
smap_summary_list <- readRDS("03_MDR_Smap_out/MDR_Smap_compiled.rds")

# fish ecological/utilization information
tax_sheet <- read.csv("01_dataformatting_out/fishinfo_w_ecology.csv")

# mean IS, latitude, and so on
analysis_df <- readRDS("08_correlations_out/analysis_df.rds")


### format data
# extract bottom temp => fishcount result
smap_df <- smap_summary_list[["mdr_surr_bt.to.fc"]] %>%
  filter(grepl("Water_temp", cause_var))

# combine smap table with fish information
smap_df_ed <- smap_df %>%
  select(Sample, effect_var, cause_var, IS, cause_var_val) %>%
  left_join(analysis_df %>% select(Fish_ID, FishBase_name, Lat_center, n_window_sig),
            by = c("effect_var" = "Fish_ID"))

# use only species with N of significant window >= 1, and arrange species by latitudinal center
plot_df <- smap_df_ed %>%
  filter(n_window_sig >= 1) %>%
  # 1. Replace underscores with spaces
  mutate(FishBase_name_clean = str_replace_all(FishBase_name, "_", " ")) %>%
  # 2. Wrap the cleaned name in italics and add latitude
  mutate(FishBase_name_ed = sprintf("*%s* (%.1f)", FishBase_name_clean, Lat_center)) %>%
  # 3. Reorder the factor levels by latitude
  mutate(FishBase_name_ed = fct_reorder(FishBase_name_ed, Lat_center))


### create plot
facet_limits <- plot_df %>%
  group_by(FishBase_name_ed) %>%
  summarise(ymax = max(abs(IS), na.rm = TRUE)*1.1) %>% 
  ungroup() %>% 
  mutate(ymin = -ymax,
         cause_var_val = 0)

# plot
plo <- ggplot(plot_df, aes(x = cause_var_val, y = IS)) +
  geom_point(alpha = 0.5, color = "grey") +
  geom_smooth(method = "gam",
              formula = y ~ s(x),
              se = TRUE, color = "black", linetype = "solid") + # loess
  facet_wrap(vars(FishBase_name_ed), scales = "free_y") + # FishBase_name 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_classic() +
  geom_blank(data = facet_limits, aes(y = ymin)) +
  geom_blank(data = facet_limits, aes(y = ymax)) +
  scale_y_continuous(labels = scaleFUN2) +
  xlim(5, 30) +
  labs(x = "Water temperature (\U2103)", y = "Dynamic response to temperature") +
  theme(strip.text = element_markdown(size = 7),
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 10)
  )


# <-------> #
# save
# <-------> #

ggsave("FigCode/FigS_temp_IS.each.sp_in.panel.png", plot = plo,
       width = 20, height = 10, bg = "white")
