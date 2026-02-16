####
#### R script for Ohigashi et al (2025)
#### Validate species shift velocity (based on citizen science data) using Maizuru data
#### 2025.10.10 written by Ohigashi
#### R 4.5.1
####


### load packages
library(ggplot2); packageVersion("ggplot2")
library(ggrepel); packageVersion("ggrepel")
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(stringr); packageVersion("stringr")
library(phyloseq); packageVersion("phyloseq")
library(ggpp); packageVersion("ggpp")
library(grid); packageVersion("grid")
source("Function/F2_HelperFunction_vizualize.R")


### loda data
# species poleward shift velocity
eo_SST_IS_shift <- readRDS("08_correlations_out/SSTchange_shift_IS_EAsiaOceania.rds")

# species abundance data
ps_all <- readRDS("01_dataformatting_out/fishdata_ps.rds")

# species info
tax_sheet <- tax_table(ps_all) %>%
  as.data.frame() %>%
  rownames_to_column("Fish_ID")


### format data
## fish abundance dynamics
# convert to melted data frame
ps_df <- ps_all %>%
  psmelt()

# pick up species with velocity data
ps_df_w.velo <- ps_df %>%
  filter(FishBase_name %in% eo_SST_IS_shift$FishBase_name)

# add continuous time data (continuous)
ps_df_w.velo <- ps_df_w.velo %>%
  mutate(time_cont = as.numeric(substr(Sample, 7, 9)))


# create plot list for the species
plot_list <- list()

# modify species names
ps_df_w.velo <- ps_df_w.velo %>%
  mutate(
    Genus = str_extract(FishBase_name, "^[^_]+"),
    Species = str_extract(FishBase_name, "(?<=_).*"),
    sci_abr = paste0(substr(Genus, 1, 1), ". ", Species),
    sci_abr_italic = paste0("italic('", substr(Genus, 1, 1), ". ", Species, "')")
  )

# create plots for each fish's abundance trend
for (sp in unique(ps_df_w.velo$FishBase_name)) {
  p <- ggplot(ps_df_w.velo %>% filter(FishBase_name == sp), aes(x = time_cont, y = log10(Abundance+0.5))) +
    geom_smooth(method = "loess", color = "black") +
    facet_wrap(~ sci_abr_italic, labeller = label_parsed) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA)) #+
  plot_list[[sp]] <- p
}


# plot for distribution center and shift velocity
vel_lat_df <- eo_SST_IS_shift %>%
  left_join(tax_sheet %>% select(FishBase_name, Lat_center, latcat),
            by = "FishBase_name")
vel_lat_df$Lat_center <- as.numeric(vel_lat_df$Lat_center)

p_lat_velo <- ggplot(vel_lat_df, aes(x = Lat_center, y = slope.w*111)) +
  geom_blank() +
  theme_classic() +
  geom_vline(xintercept = 35.47, color = "grey", linetype = "dashed") +
  labs(x = sprintf("Latitudinal center of distribution (%sN)", "\U00B0"),
       y = expression("Estimated poleward shift velocity (km year"^{'\U2212'*"1"}*")"),
       title = "Public biodiversity data-based shift velocity and abundance trends in Maizuru Bay") +
  annotate("text", x = 35.47 + 0.5, y = min(vel_lat_df$slope.w * 111),
           label = "Maizuru", hjust = 0, size = 6) + 
  theme(axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 13, color = "black"),
        title = element_text(size = 14))
p_lat_velo

# base points
pos_df <- vel_lat_df %>%
  transmute(
    FishBase_name,
    cx = Lat_center,
    sv = slope.w * 111
  ) %>%
  distinct()

# dummy plot to avoid overlapping
p_tmp <- ggplot(pos_df, aes(cx, sv, label = paste0(FishBase_name, "\nXXXXXXXX\nXXXXX"))) +
  geom_text_repel(
    box.padding   = 0.6,
    point.padding = 0.4,
    force         = 6,
    force_pull     = 0.08,
    max.iter      = 20000,
    max.time      = 5,
    max.overlaps  = Inf,
    direction      = "y",
    min.segment.length = 0,
    size = 5
  )+
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))


built <- ggplot_build(p_tmp)
repelled_xy <- built$data[[1]][, c("x", "y")]
pos_df$rx <- repelled_xy$x
pos_df$ry <- repelled_xy$y

pos_df$grob <- plot_list[pos_df$FishBase_name]

# check the plot
p_out <- p_lat_velo +
  geom_point(data = pos_df, aes(x = cx, y = sv), size = 2) +
  geom_segment(data = pos_df,
               aes(x = cx, y = sv, xend = rx, yend = ry),
               linewidth = 0.3, alpha = 0.7) +
  # add species trend plot
  ggpp::geom_plot(
    data = pos_df,
    aes(x = rx, y = ry, label = grob),
    vp.width  = unit(30, "mm"),
    vp.height = unit(20, "mm")
  ) +
  coord_cartesian(clip = "off")+
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

p_out

# adjust coordinate
pos_df <- pos_df %>%
  mutate(
    rx_jit = case_when(
    FishBase_name == "Paramonacanthus_oblongus" ~ rx - 3,
    FishBase_name == "Scorpaenodes_evides" ~ rx + 0.5,
    FishBase_name == "Takifugu_poecilonotus" ~ rx + 2.5,
    FishBase_name == "Takifugu_niphobles" ~ rx + 2.5,
    FishBase_name == "Istigobius_campbelli" ~ rx - 2.5,
    TRUE ~ rx
  ),
  ry_jit = case_when(
    FishBase_name == "Sphyraena_pinguis" ~ ry - 22.5,
    FishBase_name == "Scorpaenodes_evides" ~ ry - 7.5,
    FishBase_name == "Siganus_fuscescens" ~ ry + 5,
    FishBase_name == "Chromis_notata" ~ ry - 20,
    FishBase_name == "Parajulis_poecilepterus" ~ ry + 15,
    FishBase_name == "Sebastiscus_marmoratus" ~ ry + 30,
    FishBase_name == "Tridentiger_trigonocephalus" ~ ry - 5,
    TRUE ~ ry
  )
  )

p_out_ed <- p_lat_velo +
  geom_point(data = pos_df, aes(x = cx, y = sv), size = 2) +
  geom_segment(data = pos_df,
               aes(x = cx, y = sv, xend = rx, yend = ry),
               linewidth = 0.3, alpha = 0.7) +
  # add species trend plot
  ggpp::geom_plot(
    data = pos_df,
    aes(x = rx_jit, y = ry_jit, label = grob),
    vp.width  = unit(33, "mm"),
    vp.height = unit(20, "mm")
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

p_out_ed



### save result
ggsave("FigCode/Fig_shiftvelocity_validation_out/maizuru.abun.dynamics_shift_ed.png", plot = p_out_ed,
       width = 11, height = 8, bg = "white")

cairo_pdf("FigCode/Fig_shiftvelocity_validation_out/maizuru.abun.dynamics_shift_ed.pdf",
          width = 11, height = 8)
print(p_out_ed)
dev.off()
