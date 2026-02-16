####
#### R script for Ohigashi et al (2025)
#### Combine the plots and create figures to show in the supplementary materials
#### 2025.11.03 written by Ohigashi
#### R 4.5.2
####


### load packages and functions
library(ggplot2); packageVersion("ggplot2")
library(cowplot); packageVersion("cowplot")
source("Function/F2_HelperFunction_vizualize.R")


### create a directory to save the panels
dir.create("FigCode/FigSupple_out")


### load data and create panels (combine figures)

##### 1. latitudinal center x dynamic response (PIC version) #####
## import data
pic_latIS <- readRDS("FigCode/Fig_correlations_out/PIC_Lat_center_mean_IS.rds")

## small modification
# convert x-axis labels appropriately (digits)
pic_latIS <- pic_latIS + scale_x_continuous(labels = scaleFUN_int)

# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_pic_latIS.pdf", width = 8, height = 6)
print(pic_latIS)
dev.off()


##### 2. latitudinal center x lm or GAM based response (normal and PIC) #####
## import data
latlmres <- readRDS("FigCode/Fig_correlations_out/cor_Lat_center_lm_coef.rds")
latgamres <- readRDS("FigCode/Fig_correlations_out/cor_Lat_center_gam_slope_ave.rds")
pic_latlmres <- readRDS("FigCode/Fig_correlations_out/PIC_Lat_center_lm_coef.rds")
pic_latgamres <- readRDS("FigCode/Fig_correlations_out/PIC_Lat_center_gam_slope_ave.rds")

## small modification
# convert x-axis labels appropriately (digits)
latlmres <- latlmres + scale_x_continuous(labels = scaleFUN_int)
latgamres <- latgamres + scale_x_continuous(labels = scaleFUN_int)
pic_latlmres <- pic_latlmres + scale_x_continuous(labels = scaleFUN_int)
pic_latgamres <- pic_latgamres + scale_x_continuous(labels = scaleFUN_int)

# combine plots
figS_lat_res <- plot_grid(latlmres + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                          latgamres + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                          pic_latlmres + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                          pic_latgamres + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                          nrow = 2,
                          rel_heights = c(1, 1),
                          ncol = 2,
                          rel_widths = c(1, 1),
                          labels = c("a", "b", "c", "d"),
                          label_size = 20)
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_lat_res.pdf", width = 10, height = 10)
print(figS_lat_res)
dev.off()


##### 3. dynamic response x shift velocity (PIC version) #####
## import data
pic_ISshift <- readRDS("FigCode/Fig_correlations_out/PIC_mean_IS_slope.w.rds")

## small modification
# convert x-axis labels appropriately (digits)
pic_ISshift <- pic_ISshift + scale_x_continuous(labels = scaleFUN)

# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_pic_ISshift.pdf", width = 8, height = 6)
print(pic_ISshift)
dev.off()


##### 4. lm or GAM based response x shift velocity (normal and PIC) #####
## import data
lmresshift <- readRDS("FigCode/Fig_correlations_out/cor_lm_coef_slope.w.rds")
gamresshift <- readRDS("FigCode/Fig_correlations_out/cor_gam_slope_ave_slope.w.rds")
pic_lmresshift <- readRDS("FigCode/Fig_correlations_out/PIC_lm_coef_slope.w.rds")
pic_gamresshift <- readRDS("FigCode/Fig_correlations_out/PIC_gam_slope_ave_slope.w.rds")

## small modification
# convert x-axis labels appropriately (digits)
lmresshift <- lmresshift +
  scale_x_continuous(labels = scaleFUN)+
  scale_y_continuous(labels = scaleFUN_int)
gamresshift <- gamresshift +
  scale_x_continuous(labels = scaleFUN)+
  scale_y_continuous(labels = scaleFUN_int)
pic_lmresshift <- pic_lmresshift + scale_x_continuous(labels = scaleFUN)
pic_gamresshift <- pic_gamresshift  + scale_x_continuous(labels = scaleFUN)

# combine plots
figS_res_shift <- plot_grid(lmresshift + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                            gamresshift + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                            pic_lmresshift + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                            pic_gamresshift + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
                            nrow = 2,
                            rel_heights = c(1, 1),
                            ncol = 2,
                            rel_widths = c(1, 1),
                            labels = c("a", "b", "c", "d"),
                            label_size = 20)
# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_res_shift.pdf", width = 10, height = 10)
print(figS_res_shift)
dev.off()


##### 5. Temperature x richness or abundance #####
## import data
temp_abun <- readRDS("FigCode/Fig_AbunRich_dynamics_out/abun_temp_latcat.rds")
temp_rich <- readRDS("FigCode/Fig_AbunRich_dynamics_out/rich_temp_latcat.rds")

# get legend
leg_latcat <- get_legend(temp_abun)

# combine plots
figS_temp_abunrich <- plot_grid(temp_abun + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                                                 legend.position = "none"),
                                temp_rich + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                                                  legend.position = "none"),
                                leg_latcat,
                                ncol = 3,
                                rel_widths = c(1, 1, 0.4),
                                labels = c("a", "b", NA),
                                label_size = 20)

# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_temp_abunrich.pdf", width = 10, height = 5)
print(figS_temp_abunrich)
dev.off()


##### 6. prediction by IS considering interspecific interactions #####
## import data
IS.int_shift <- readRDS("FigCode/FigS_correlation_w_interspecific_interactions_out/cor_mean_IS_shift_fishapes_w.sp.int.rds")
vali_IS.int_IS <- readRDS("09_interspecific_interaction_out/validation_IS_and_IS.w.species.interaction.rds")

# combine plots
figS_IS.int <- plot_grid(IS.int_shift + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                                                  legend.position = "none"),
                         vali_IS.int_IS + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                                                  legend.position = "none"),
                                nrow = 1,
                                rel_widths = c(1.5, 1),
                                labels = c("a", "b"),
                                label_size = 20)

# save PDF
cairo_pdf("FigCode/FigSupple_out/FigS_IS.int.pdf", width = 11.5, height = 5)
print(figS_IS.int)
dev.off()


