####
#### R script for Ohigashi et al (2025)
#### Analysis of causal effect from temperature to fish using UIC, considering interspecific interactions
#### 2025.11.26 written by Ohigashi
#### R 4.5.2
#### 


### load packages
library(phyloseq); packageVersion("phyloseq")
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(rUIC); packageVersion("rUIC") # 0.9.12
library(rEDM); packageVersion("rEDM") # 0.7.5
library(macamts); packageVersion("macamts") # 0.1.4
library(pbapply); packageVersion("pbapply")
source("Function/F1_HelperFunction_stats.R")


### load data
# fish abundance dynamics
ps_all <- readRDS("01_dataformatting_out/fishdata_ps.rds")

# # fish taxonomy
fishinfo <- read.csv("01_dataformatting_out/fishinfo_w_ecology.csv")


### check observed frequency of species
# convert to melted data frame
ps_df <- ps_all %>%
  psmelt()

# rank of species frequency observed
top_sp_obs <- ps_df %>%
  group_by(FishBase_name) %>%
  summarise(obs_count = sum(Abundance >= 1)) %>%
  arrange(desc(obs_count))

# extract species that were observed at least 5% of the sampling
total_events <- n_distinct(ps_df$Sample)
focal_fish <- top_sp_obs %>%
  filter(obs_count >= total_events * 0.05) # 46 fish extracted


##### multi-species UIC conditioned by water temperature #####
### format data
# fish count data
fishcount <- otu_table(ps_all) %>% as.data.frame()
fishcount <- fishcount %>% rownames_to_column("Sample")

# extract focal fish counts
focal_fish_id <- focal_fish %>%
  left_join(fishinfo, by = "FishBase_name") %>%
  pull(Fish_ID)
focal_fishcount <- fishcount %>%
  select(Sample, any_of(focal_fish_id))

# sample data
sampledata <- data.frame(sample_data(ps_all))
sampledata <- sampledata %>% rownames_to_column("Sample")

# select focal variables
sample_df <- sampledata %>%
  # use only bottom temperature and remove unused variables
  select(Sample, Water_temp_bottom)

# combine fish and sample data
sample_fishcount_df <- merge(sample_df, focal_fishcount, by = "Sample", sort = F)

df_std <- as.data.frame(apply(sample_fishcount_df[, -1], 2, function(x) as.numeric(scale(x))))

## loop to calculate how other species affected one species
# prepare output object
# get Fish IDs
sp_id <- colnames(df_std)[-1]

# loop
uic_sp_int <- pblapply(
  sp_id,
  cl = 32,
  FUN = function(id) {
    
    effect_var <- id
    
    # perform UIC for other species => effect var (influenced species)
    set.seed(123) 
    uic_res <- uic_across_TOed(
      df_std,
      effect_var,
      E_range = 0:24,
      tp_range = -12:0,
      silent = TRUE,
      cond_var = "Water_temp_bottom"
    )
    
    return(uic_res)
  }
)

names(uic_sp_int) <- sp_id

### save result
dir.create("09_interspecific_interaction_out")

saveRDS(uic_sp_int, "09_interspecific_interaction_out/intersp_effect_list_condv.temp.rds")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/09_1_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))
