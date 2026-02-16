####
#### R script for Ohigashi et al (2025)
#### NMDS plot for fish with latitudinal category
#### 2025.05.23 written by Ohigashi
#### R 4.5.0
####


### load packages
library(phyloseq); packageVersion("phyloseq")
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(ggplot2); packageVersion("ggplot2")
library(ggsci); packageVersion("ggsci")
library(ggrepel); packageVersion("ggrepel")
library(vegan); packageVersion("vegan")
library(cowplot); packageVersion("cowplot")
source("Function/F2_HelperFunction_vizualize.R")


### load data
# fish community data all time
ps_all <- readRDS("01_dataformatting_out/fishdata_ps.rds")

### plot NMDS
set.seed(123) # set rng.seed

## all seasons at once
ps_all.bray <- ordinate(ps_all, "NMDS", "bray")

# fit abundance and richness information to the ordination
# extract temp and year vars
env_all <- data.frame(Sample = sample_names(ps_all),
                      Temp = sample_data(ps_all)$Water_temp_bottom,
                      Year = sample_data(ps_all)$Year,
                      lat.high.total = sample_data(ps_all)$lat.high.total,
                      lat.low.total = sample_data(ps_all)$lat.low.total,
                      lat.high.rich = sample_data(ps_all)$lat.high.rich,
                      lat.low.rich = sample_data(ps_all)$lat.low.rich
)
# add mean temperature for each year
env_all <- env_all %>%
  group_by(Year) %>%
  mutate(Temp.yr.ave = mean(Temp)) %>%
  ungroup()


ef <- envfit(ps_all.bray, env_all[,-1], permu = 99999)
ef_pvals <- as.data.frame(ef$vectors$pvals) # extract p-values
ef_r <- as.data.frame(ef$vectors$r) # extract r2
ef_arrows <- as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r)) # add weight by r
ef_table <- cbind(ef_arrows, ef_pvals); colnames(ef_table)[3] <- "pvals"
ef_table <- ef_table %>% rownames_to_column(var = "factor")
ef_table <- ef_table %>%
  mutate(factor = case_when(
    factor == "lat.high.total" ~ "High.lat.total",
    factor == "lat.low.total" ~ "Low.lat.total",
    factor == "lat.high.rich" ~ "High.lat.rich",
    factor == "lat.low.rich" ~ "Low.lat.rich",
    TRUE ~ factor
  ),
  factor_sig = case_when(
    pvals < 0.05 ~ paste0(factor, "*"),
    TRUE ~ factor
  ),
  year_div = "2002~2006" # dummy data
  )

# plot
pnmds_all <- plot_ordination(ps_all, ps_all.bray, color = "Season", shape = "year_div") +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill=Season, group = Season)) +
  geom_point(size = 2) +
  scale_color_manual(values = c(
    "Spring" = "green", 
    "Summer" = "orange", 
    "Fall" = "brown", 
    "Winter" = "blue"
  )) +
  scale_fill_manual(values = c(
    "Spring" = "green", 
    "Summer" = "orange", 
    "Fall" = "brown", 
    "Winter" = "blue"
  )) +
  scale_shape_manual(values = c(16, 17, 15, 1, 2)) +
  scale_y_continuous(labels = scaleFUN) +
  scale_x_continuous(labels = scaleFUN) +
  theme_bw() +
  theme(panel.border = element_rect(size = 1.2, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        plot.title = element_text(size=14, color = "black"),
        axis.text=element_text(size=12, color = "black"),
        axis.title=element_text(size=14, color = "black"))+
  labs(title = "Fish community",
       fill = "Season", color = "Season", shape = "Year")+
  geom_segment(data = ef_table, aes(x=0, y=0, xend=NMDS1*0.9, yend=NMDS2*0.9),
               size=0.5, arrow=arrow(type = "open", length = unit(0.08, "inches")), color="black", show.legend = F) +
  geom_label_repel(data = ef_table, aes(x = NMDS1, y = NMDS2, label = factor_sig),
                   size = 3, color="black", show.legend = F)

# plot only with Year, Temp, and Temp.yr.ave
pnmds_sub <- plot_ordination(ps_all, ps_all.bray, color = "Season", shape = "year_div") +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill=Season, group = Season)) +
  geom_point(size = 2) +
  scale_color_manual(values = c(
    "Spring" = "green", 
    "Summer" = "orange", 
    "Fall" = "brown", 
    "Winter" = "blue"
  )) +
  scale_fill_manual(values = c(
    "Spring" = "green", 
    "Summer" = "orange", 
    "Fall" = "brown", 
    "Winter" = "blue"
  )) +
  scale_shape_manual(values = c(16, 17, 15, 1, 2)) +
  scale_y_continuous(labels = scaleFUN) +
  scale_x_continuous(labels = scaleFUN) +
  theme_bw() +
  theme(panel.border = element_rect(size = 1.2, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        plot.title = element_text(size=14, color = "black"),
        axis.text=element_text(size=12, color = "black"),
        axis.title=element_text(size=14, color = "black"))+
  labs(title = "Fish community",
       fill = "Season", color = "Season", shape = "Year")+
  geom_segment(data = ef_table %>% filter(factor %in% c("Year", "Temp", "Temp.yr.ave")), aes(x=0, y=0, xend=NMDS1*0.9, yend=NMDS2*0.9),
               size=0.5, arrow=arrow(type = "open", length = unit(0.08, "inches")), color="black", show.legend = F) +
  geom_label_repel(data = ef_table %>% filter(factor %in% c("Year", "Temp", "Temp.yr.ave")), aes(x = NMDS1, y = NMDS2, label = factor_sig),
                   size = 3, color="black", show.legend = F)


### save data
dir.create("FigCode/Fig_NMDS_out")
saveRDS(pnmds_all, "FigCode/Fig_NMDS_out/Fig_NMDS_allvars.rds")
saveRDS(pnmds_sub, "FigCode/Fig_NMDS_out/Fig_NMDS_tempyear.rds")
ggsave("FigCode/Fig_NMDS_out/Fig_NMDS_allvars.png", plot = pnmds_all,
       width = 8, height = 7, bg = "white")
ggsave("FigCode/Fig_NMDS_out/Fig_NMDS_tempyear.png", plot = pnmds_sub,
       width = 8, height = 7, bg = "white")



