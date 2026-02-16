####
#### R script for Ohigashi et al (2025)
#### Plot state space of Trachurus japonicus for a concept diagram (state-dependent, dynamic, etc.)
#### 2025.05.23 written by Ohigashi
#### R 4.5.0
####


### load packages
library(plotly); packageVersion("plotly")
library(phyloseq); packageVersion("phyloseq")
library(tidyverse); packageVersion("tidyverse")

library(reticulate)
path_to_python <- "XXXX"
use_miniconda(path_to_python)


### load data
ps_all <- readRDS("01_dataformatting_out/fishdata_ps.rds")


### format data
# get target taxa 's ID
target_fish_id <- tax_table(ps_all) %>%
  as.data.frame() %>%
  rownames_to_column("Fish_ID") %>%
  filter(FishBase_name == "Pseudolabrus_sieboldi") %>%
  pull(Fish_ID)

# extract time series of T. japonicus
abun_ts <- otu_table(ps_all)[, target_fish_id] %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  rename(Abundance = !!target_fish_id) %>% # change column name
  arrange(Sample) %>%
  mutate(Abundance_lag1 = lag(Abundance, 1)) # add one-timepoint delayed abundance

# extract time series of bottom water temperature
temp_ts <- sample_data(ps_all) %>%
  as("data.frame") %>%
  rownames_to_column("Sample") %>%
  select(Sample, Water_temp_bottom) %>%
  mutate(Water_temp_bottom_lag1 = lag(Water_temp_bottom, 1)) # add one-timepoint delayed abundance

# combine them
ts_df <- left_join(abun_ts, temp_ts, by = "Sample")


### linear model to show on the 3D figure
fit <- lm(Abundance ~ Water_temp_bottom, data = ts_df[97:193,])

x_seq <- seq(min(ts_df[97:193,]$Water_temp_bottom), max(ts_df[97:193,]$Water_temp_bottom), length.out = 100)
z_pred <- predict(fit, newdata = data.frame(Water_temp_bottom = x_seq))
y_fixed <- min(ts_df[97:193,]$Abundance_lag1)

line_df <- data.frame(x = x_seq, y = y_fixed, z = z_pred)


### plot 3D figure
df <- ts_df[97:193, ]

df <- df %>%
  mutate(
    u = lead(Water_temp_bottom) - Water_temp_bottom,
    v = lead(Abundance_lag1) - Abundance_lag1,
    w = lead(Abundance) - Abundance
  ) %>%
  na.omit()


p2 <- plot_ly() %>%
  add_trace(
    type = "cone",
    x = df$Water_temp_bottom,
    y = df$Abundance_lag1,
    z = df$Abundance,
    u = df$u,
    v = df$v,
    w = df$w,
    sizemode = "absolute",
    sizeref = 70,
    colorscale = list(c(0, 1), c("orange", "orange")),
    showscale = FALSE
  ) %>% 
  add_trace(
    type = "scatter3d",
    mode = "lines+marker",
    x = df$Water_temp_bottom,
    y = df$Abundance_lag1,
    z = df$Abundance,
    line = list(color = 'grey', width = 2),
    showlegend = FALSE
  ) %>%
  add_trace(data = line_df, x = ~x, y = ~y, z = ~z, 
            type = "scatter3d", mode = "lines", line = list(color = 'red', width = 4, linetype = "dashed"),
            name = "Regression line", showlegend = FALSE) %>% 
  layout(
    scene = list(
      xaxis = list(title = "Water temperature [t]", titlefont = list(size = 25),
                   showline= T, linewidth=2, linecolor="black"),
      yaxis = list(title = "Abundance [t-1]", titlefont = list(size = 25),
                   showline= T, linewidth=2, linecolor="black"),
      zaxis = list(title = "Abundance [t]", titlefont = list(size = 25),
                   showline= T, linewidth=2, linecolor="black")
    )
  )

# set view point of the 3d plot
p_camera <- p2 %>%
  layout(scene = list(
    camera = list(
      eye = list(x = 1.44, y = -1.92, z = 1.14)
    )
  ))
p_camera

# save plot
dir.create("FigCode/Fig_concept_out")
plotly::save_image(p_camera, file = "FigCode/Fig_concept_out/concept_3d.png",
                   width = 800, height = 700, scale = 2)


### plot 2D figure
p_2d <- ggplot(df, aes(x = Water_temp_bottom, y = Abundance)) +
  geom_point(color = "orange", size = 3, alpha = 0.3) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(x = "Water temperature", y = "Abundance") +
  theme_classic() +
  theme(axis.text=element_text(size=16, color = "black"),
        axis.title=element_text(size=16, color = "black")) 
p_2d

# save plot
ggsave("FigCode/Fig_concept_out/concept_2d.png",
       width = 4, height = 3, bg = "white")
