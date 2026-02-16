####
#### R script for Ohigashi et al (2025)
#### plot PIC-based correlations between variables
#### 2025.06.17 written by Ohigashi
#### R 4.5.0
####


### load packages
library(tidyverse); packageVersion("tidyverse")
library(ape); packageVersion("ape")
source("Function/F2_HelperFunction_vizualize.R")


### load data
# mean IS, latitude, and so on
analysis_df <- readRDS("08_correlations_out/analysis_df.rds")

# correlation result
cor_pic <- read.csv("08_correlations_out/correlation_PIC.csv")

# phylogenetic tree of Maizuru fish
fishtree <- read.tree("05_phylogeny_out/RAxML_tree_w.B.nwk")


### visualize correlations between variables
# create plot list
cor_pic_plotlist <- list()

# make a list of variable pairs
var_pairs <- list(
  c("Lat_center", "mean_IS"),
  c("Lat_center", "lm_coef"),
  c("Lat_center", "gam_slope_ave"),
  c("mean_IS", "slope.w"),
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
  r.tmp <- cor_pic %>%
    filter(v1 == !!v1, v2 == !!v2) %>% pull(PIC_r)
  p.tmp <- cor_pic %>%
    filter(v1 == !!v1, v2 == !!v2) %>% pull(PIC_p)
  anno.tmp <- paste0("italic(r) == ", scaleFUN2(r.tmp), ' * ", " * ', p_to_sig(p.tmp, type = "threshold"))
  anno.tmp <- sub("\U2212", "-", anno.tmp)
  
  # set axis titles
  if (v1 == "Lat_center") {
    xtitle <- "Latitudinal center of distribution (PIC)"
  } else if (v1 == "lm_coef") {
    xtitle <- "LM-based response to temperature (PIC)"
  } else if (v1 == "gam_slope_ave") {
    xtitle <- "GAM-based response to temperature (PIC)"
  } else if (v1 == "mean_IS") {
    xtitle <- "Mean dynamic response to temperature (PIC)"
  } else {
    xtitle <- v1
  }
  if (v2 == "slope.w") {
    ytitle <- "Estimated poleward shift velocity (PIC)"
  } else if (v2 == "lm_coef") {
    ytitle <- "LM-based response to temperature (PIC)"
  } else if (v2 == "gam_slope_ave") {
    ytitle <- "GAM-based response to temperature (PIC)"
  } else if (v2 == "mean_IS") {
    ytitle <- "Mean dynamic response to temperature (PIC)"
  } else {
    ytitle <- v2
  }
  
  # convert data to PIC values
  pruned_tree <- drop.tip(fishtree, setdiff(fishtree$tip.label, sub_df$FishBase_name))
  # order data frame as the tree is
  sub_df_reo <- sub_df[match(pruned_tree$tip.label,sub_df$FishBase_name),]
  
  # Calculate PICs
  pic_v1 <- pic(sub_df_reo[[v1]], pruned_tree)
  pic_v2 <- pic(sub_df_reo[[v2]], pruned_tree)
  pic_df <- data.frame(pic_v1 = pic_v1, pic_v2 = pic_v2)
  
  # plot
  plo.tmp <- ggplot(pic_df,
                    aes(x = pic_v1, y = pic_v2
                    )) +
    geom_point() +
    theme_classic() +
    theme(legend.title = element_text(size = 12, color = "black"),
          legend.text = element_text(size = 12, color = "black"),
          axis.text=element_text(size=14, color = "black"),
          axis.title=element_text(size=14, color = "black")) +
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
  if (p.tmp < 0.05) {
    plo.tmp.fin <- plo.tmp.fin +
      geom_smooth(method = "lm", se = TRUE, color = "black")
  }
  
  # store to the list
  plotname <- paste(v1, v2, sep = "_")
  cor_pic_plotlist[[plotname]] <- plo.tmp.fin
}


# save plots as RDS
for (name in names(cor_pic_plotlist)) {
  saveRDS(cor_pic_plotlist[[name]], file = sprintf("FigCode/Fig_correlations_out/PIC_%s.rds", name))
}
