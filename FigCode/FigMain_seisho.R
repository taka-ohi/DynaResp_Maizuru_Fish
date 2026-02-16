####
#### R script for Ohigashi et al (2025)
#### Combine the plots and create figures to show in the main article
#### 2025.06.17 written by Ohigashi
#### R 4.5.0
####


### load packages and functions
library(ggplot2); packageVersion("ggplot2")
library(cowplot); packageVersion("cowplot")
library(ggbreak); packageVersion("ggbreak")


### create a directory to save the panels
dir.create("FigCode/FigMain_out")


### load data and create panels (combine figures)

##### 1. Abundance dynamics, NMDS, and northern/southern species dynamics #####
## import data
top5_abun <- readRDS("FigCode/Fig_AbunRich_dynamics_out/abun_top5_others.rds")
nmds_all <- readRDS("FigCode/Fig_NMDS_out/Fig_NMDS_tempyear.rds")
abun_latcat <- readRDS("FigCode/Fig_AbunRich_dynamics_out/abun_ts_latcat.rds")
rich_latcat <- readRDS("FigCode/Fig_AbunRich_dynamics_out/rich_ts_latcat.rds")


top5_abun_ed2 <- top5_abun +
  coord_cartesian(ylim = c(0, 1800)) + # break the plot
  theme(legend.position = c(0.65, 0.75),
        legend.text = element_text(size = 11, color = "black")#,
        # legend.spacing = unit(0.8, "cm")
        ) # show legend inside the plot


# extract the legend for latitude category
latcat_leg <- get_legend(abun_latcat)

# combine them
figM1 <- plot_grid(top5_abun_ed2 + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   plot_grid(nmds_all + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                             abun_latcat + theme(legend.position = "none",
                                                 plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                                                 axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)),
                             rich_latcat + theme(legend.position = "none",
                                                 plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                                                 
                                                 axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)),
                             latcat_leg,
                             nrow = 1,
                             rel_widths = c(1.5, 0.9, 0.9, 0.4),
                             labels = c("b", "c", "d", NA),
                             label_size = 20
                             ),
                   nrow = 2,
                   rel_heights = c(1, 1),
                   labels = c("a", NA),
                   label_size = 20
                   )

# save PDF
cairo_pdf("FigCode/FigMain_out/FigM_abun_nmds_latcat.pdf", width = 16, height = 10)
print(figM1)
dev.off()
# save png
ggsave(filename = "FigCode/FigMain_out/FigM_abun_nmds_latcat.png",
       plot = figM1, width = 16.5, height = 10, bg = "white")


##### 2. correlation between latitudinal distribution center and mean IS #####
## import data
cor_lat_IS <- readRDS("FigCode/Fig_correlations_out/cor_Lat_center_mean_IS_fishapes.rds")

# save PDF
cairo_pdf("FigCode/FigMain_out/FigM_corrlation_latcenter_IS.pdf", width = 8, height = 6)
print(cor_lat_IS)
dev.off()


##### 3. correlation between mean IS and shift velocity, and SST change #####
## import data
cor_IS_shift <- readRDS("FigCode/Fig_correlations_out/cor_mean_IS_shift_fishapes.rds")
sst_shift_IS <- readRDS("FigCode/Fig_correlations_out/SST_shift_IS_continuous.rds")

# combine plots
figM2 <- plot_grid(cor_IS_shift + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   sst_shift_IS + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                   nrow = 2,
                   rel_heights = c(1, 1),
                   # ncol = 2,
                   # rel_widths = c(1, 1),
                   labels = c("a", "b"),
                   label_size = 20)
# figM2

# save PDF
cairo_pdf("FigCode/FigMain_out/FigM_corrlation_IS_shift_SST.pdf", width = 8, height = 10)
print(figM2)
dev.off()
# save png
ggsave(filename = "FigCode/FigMain_out/FigM_corrlation_IS_shift_SST.png",
       plot = figM2, width = 8, height = 10, bg = "white")



