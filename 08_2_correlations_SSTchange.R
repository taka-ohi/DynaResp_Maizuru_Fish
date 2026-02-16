####
#### R script for Ohigashi et al (2025)
#### relationships among IS, SST change, and shift
#### 2025.06.17 written by Ohigashi
#### R 4.5.0
####


### load packages
library(raster); packageVersion("raster")
library(dplyr); packageVersion("dplyr")
library(purrr); packageVersion("purrr")
library(ggplot2); packageVersion("ggplot2")
library(akima); packageVersion("akima")


### load data
# SST data
# downloaded from here: https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html 
sst <- brick("Data/HadISST_sst.nc")

# sp. shift velocity, IS, and moving-UIC
analysis_df <- readRDS("08_correlations_out/analysis_df.rds")


### calculate SST change at the mean latitude where the species observed
# subset SST data pf 1960~2024
layer_names <- names(sst)
years <- as.numeric(substr(layer_names, 2, 5))
target_idx <- which(years >= 1960 & years <= 2024)
sst_target <- sst[[target_idx]]
years_target <- years[target_idx]


# function to calculate slop of SST
get_sst_slope <- function(lat, sst_brick, years_vec, long_min, long_max) {
  if (is.na(lat)) return(NA_real_)  # return NA if the species not observed
  
  # crop by the latitude where the species located
  sst_subset <- crop(sst_brick, extent(long_min, long_max, lat - 1, lat + 1))
  
  # spatial average
  sst_mean_ts <- cellStats(sst_subset, stat = "mean")
  
  # average by year
  sst_annual <- tapply(sst_mean_ts, years_vec, mean)
  year_numeric <- as.numeric(names(sst_annual))
  
  # calculate slope
  if (length(year_numeric) < 10) return(NA_real_)
  slope <- coef(lm(sst_annual ~ year_numeric))[2]
  return(slope)
}


# add annual SST change (degree C / yr) at the mean latitude where the species observed (East Asia & Oceania)
ana_df_SST <- analysis_df %>%
  mutate(
    SST_slope_Northern = map_dbl(mean_lat_hemiNorthern, ~ get_sst_slope(.x, sst_target, years_target, 105, 180)),
    SST_slope_Southern = map_dbl(mean_lat_hemiSouthern, ~ get_sst_slope(.x, sst_target, years_target, 105, 180))
  )


### linear regression considering interaction 
# East Asia & Oceania (Average SST change of the sp. distribution in North and South hemispheres)
# average the SST change if the sp. in both hemispheres, otherwise put one of the values
eo_SST_IS_shift <- ana_df_SST %>%
  mutate(SST_slope_NS_ave = rowMeans(
    dplyr::select(., SST_slope_Northern, SST_slope_Southern),
    na.rm = TRUE
  ))

eo_SST_IS_shift <- eo_SST_IS_shift %>%
  dplyr::select(FishBase_name, SST_slope_NS_ave, mean_IS, slope.w, n_window_sig) %>%
  filter(n_window_sig >= 1) %>% # extract the species with at least one time moving-window significance (adjusted p < 0.05)
  na.omit()

mod_eo_ave <- lm(slope.w ~ SST_slope_NS_ave * mean_IS, data = eo_SST_IS_shift)
mod_eo_ave_summary <- summary(mod_eo_ave) # significant interaction between sst change and mean IS
mod_eo_ave_summary <- coef(mod_eo_ave_summary)
mod_eo_ave_coef <- coef(mod_eo_ave)




### save data
# mean IS and shift velocity
saveRDS(eo_SST_IS_shift, "08_correlations_out/SSTchange_shift_IS_EAsiaOceania.rds")
write.csv(mod_eo_ave_summary, "08_correlations_out/lm_shift_SSTchange_IS_EAsiaOceania.csv", quote = F)

# all shift velocity and SST change
saveRDS(ana_df_SST, "08_correlations_out/SSTchange_shift_IS_EAsiaOceania_all.sp.rds")

### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/08_2_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))




