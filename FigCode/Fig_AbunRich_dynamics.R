####
#### R script for Ohigashi et al (2025)
#### Plot abundance dynamics of fish
#### 2025.05.23 written by Ohigashi
#### R 4.5.0
####


### load packages
library(phyloseq); packageVersion("phyloseq")
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(ggplot2); packageVersion("ggplot2")
library(cowplot); packageVersion("cowplot")
library(ggbreak); packageVersion("ggbreak")
library(purrr); packageVersion("purrr")
library(tidyverse); packageVersion("tidyverse")


### load data
ps_all <- readRDS("01_dataformatting_out/fishdata_ps.rds")


### 1. dynamics of top 5 species and others
# convert to melted data frame
ps_df <- ps_all %>%
  psmelt()

# extract top 5 fish in abundance
top5_taxa <- ps_df %>%
  group_by(FishBase_name) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 5) %>%
  pull(FishBase_name)

# make other fish "Others
ps_df <- ps_df %>%
  mutate(Taxon_label = if_else(FishBase_name %in% top5_taxa, FishBase_name, "Others"))

# create a data frame for plotting
ps_plot_df <- ps_df %>%
  group_by(Sample, Taxon_label) %>%
  summarise(abundance = sum(Abundance), .groups = "drop") %>%
  mutate(log_abundance = log(abundance + 0.5))

# add continuous time data (continuous)
ps_plot_df <- ps_plot_df %>%
  mutate(time_cont = as.numeric(substr(Sample, 7, 9)))

# set levels of species
ps_plot_df$Taxon_label <- factor(ps_plot_df$Taxon_label, levels = c(top5_taxa, "Others"))


# convert taxa labels to italic
taxa_levels <- levels(ps_plot_df$Taxon_label)
taxa_labels_expr <- map(taxa_levels, function(taxon) {
  if (taxon == "Others") {
    expression("Others")
  } else {
    parts <- strsplit(taxon, "_")[[1]]
    genus_initial <- paste0(substr(parts[1], 1, 1), ".")  # "S."
    bquote(italic(.(genus_initial) ~ .(parts[2])))
    # bquote(italic(.(parts[1]) ~ .(parts[2])))
  }
})

# set color vector
color_vec <- c(scales::hue_pal()(length(taxa_levels) - 1), "grey")


# plot
p_abun_top5 <- ggplot(ps_plot_df, aes(x = time_cont, y = abundance, color = Taxon_label)) +
  geom_point(size = 1, alpha = 0.3) +
  geom_line(aes(group = Taxon_label)) +
  theme_classic() +
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 10, color = "black"),
        title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 10, color = "black"),
        legend.position = "bottom",
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank()
  ) +
  labs(x = NULL, y = "Abundance", color = NULL) +
  geom_hline(yintercept = 2800, color = "grey", linetype = "dotdash") +
  scale_x_continuous(breaks = c(1, 73, 145, 217, 289, 361, 433, 505),
                     labels = c("2002-Jan", "2005-Jan", "2008-Jan", "2011-Jan",
                                "2014-Jan", "2017-Jan", "2020-Jan", "2023-Jan")) +
  scale_y_break(c(1800, 2800), scales = 0.1, ticklabels = c(0, 500, 1000, 1500, 3000, 4000)) +
  scale_color_manual(
    values = color_vec,
    labels = taxa_labels_expr 
  ) +
  guides(color = guide_legend(nrow = 1))

p_abun_top5


### 2. dynamics with latitudinal category (high and low)
# format to plot by latitudinal category
sample_sheet <- sample_data(ps_all) %>%
  as(., "data.frame") %>% 
  rownames_to_column("Sample")

sample_long <- sample_sheet %>%
  pivot_longer(
    cols = c(lat.high.total, lat.low.total, lat.high.rich, lat.low.rich),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    latcat = case_when(
      str_detect(variable, "high") ~ "High latitude sp.",
      str_detect(variable, "low") ~ "Low latitude sp."
    ),
    metric = case_when(
      str_detect(variable, "total") ~ "abundance",
      str_detect(variable, "rich") ~ "richness"
    )
  ) %>%
  select(Sample, Water_temp_bottom, latcat, metric, value)

abun_sheet <- sample_long %>%
  filter(metric == "abundance") %>%
  mutate(time_cont = as.numeric(substr(Sample, 7, 9)))

rich_sheet <- sample_long %>%
  filter(metric == "richness") %>%
  mutate(time_cont = as.numeric(substr(Sample, 7, 9)))

## 2-1. abundance x time
p_abun.ts.latcat <- ggplot(abun_sheet %>% filter(value > 0), aes(x = time_cont, y = log(value), color = latcat, group = latcat)) +
  geom_point(alpha = 0.1) +
  geom_smooth(aes(group = latcat, color = latcat, fill = latcat),
              method = "loess", se = TRUE, alpha = 0.2) +
  theme_classic() +
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 12, color = "black")
  ) +
  labs(x = NULL, y = "Log(total abundance + 0.5)", color = NULL, fill = NULL,
       title = "Abundance") +
  scale_x_continuous(breaks = c(49, 193, 337, 481),
                     labels = c("2004-Jan", "2010-Jan", "2016-Jan", "2022-Jan")) +
  scale_color_manual(values = c("High latitude sp."="blue", "Low latitude sp."="red")) +
  scale_fill_manual(values = c("High latitude sp." = "blue", "Low latitude sp." = "red"))

## 2-2. richness x time
p_rich.ts.latcat <- ggplot(rich_sheet, aes(x = time_cont, y = value, color = latcat, group = latcat)) +
  geom_point(alpha = 0.1) +
  geom_smooth(aes(group = latcat, color = latcat, fill = latcat),
              method = "loess", se = TRUE, alpha = 0.2) +
  theme_classic() +
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 12, color = "black")
  ) +
  labs(x = NULL, y = "Number of species", color = NULL, fill = NULL,
       title = "Species richness") +
  scale_x_continuous(breaks = c(49, 193, 337, 481),
                     labels = c("2004-Jan", "2010-Jan", "2016-Jan", "2022-Jan")) +
  scale_color_manual(values = c("High latitude sp."="blue", "Low latitude sp."="red")) +
  scale_fill_manual(values = c("High latitude sp." = "blue", "Low latitude sp." = "red"))

## 2-3. abundance x temperature
p_abun.temp.latcat <- ggplot(abun_sheet %>% filter(value > 0), aes(x = Water_temp_bottom, y = log(value), color = latcat, group = latcat)) +
  geom_point(alpha = 0.1) +
  geom_smooth(aes(group = latcat, color = latcat, fill = latcat),
              method = "loess", se = TRUE, alpha = 0.2) +
  theme_classic() +
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 12, color = "black")
  ) +
  labs(x = sprintf("Water temperature (%sC)", "\U00B0"), y = "Log(total abundance + 0.5)", color = NULL, fill = NULL,
       title = "Abundance") +
  scale_color_manual(values = c("High latitude sp."="blue", "Low latitude sp."="red")) +
  scale_fill_manual(values = c("High latitude sp." = "blue", "Low latitude sp." = "red"))

## 2-4. richness x temperature
p_rich.temp.latcat <- ggplot(rich_sheet, aes(x = Water_temp_bottom, y = value, color = latcat, group = latcat)) +
  geom_point(alpha = 0.1) +
  geom_smooth(aes(group = latcat, color = latcat, fill = latcat),
              method = "loess", se = TRUE, alpha = 0.2) +
  theme_classic() +
  theme(axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 12, color = "black")
  ) +
  labs(x = sprintf("Water temperature (%sC)", "\U00B0"), y = "Number of species", color = NULL, fill = NULL,
       title = "Species richness") +
  scale_color_manual(values = c("High latitude sp."="blue", "Low latitude sp."="red")) +
  scale_fill_manual(values = c("High latitude sp." = "blue", "Low latitude sp." = "red"))



### save plot
dir.create("FigCode/Fig_AbunRich_dynamics_out")
# top 5 fish
saveRDS(p_abun_top5, "FigCode/Fig_AbunRich_dynamics_out/abun_top5_others.rds")
ggsave("FigCode/Fig_AbunRich_dynamics_out/abun_top5_others.png", plot = p_abun_top5,
       width = 10, height = 5, bg = "white")

# latitude category
# time series
saveRDS(p_abun.ts.latcat, "FigCode/Fig_AbunRich_dynamics_out/abun_ts_latcat.rds")
ggsave("FigCode/Fig_AbunRich_dynamics_out/abun_ts_latcat.png", plot = p_abun.ts.latcat,
       width = 8, height = 7, bg = "white")
saveRDS(p_rich.ts.latcat, "FigCode/Fig_AbunRich_dynamics_out/rich_ts_latcat.rds")
ggsave("FigCode/Fig_AbunRich_dynamics_out/rich_ts_latcat.png", plot = p_rich.ts.latcat,
       width = 8, height = 7, bg = "white")
# temperature
saveRDS(p_abun.temp.latcat, "FigCode/Fig_AbunRich_dynamics_out/abun_temp_latcat.rds")
ggsave("FigCode/Fig_AbunRich_dynamics_out/abun_temp_latcat.png", plot = p_abun.temp.latcat,
       width = 8, height = 7, bg = "white")
saveRDS(p_rich.temp.latcat, "FigCode/Fig_AbunRich_dynamics_out/rich_temp_latcat.rds")
ggsave("FigCode/Fig_AbunRich_dynamics_out/rich_temp_latcat.png", plot = p_rich.temp.latcat,
       width = 8, height = 7, bg = "white")


