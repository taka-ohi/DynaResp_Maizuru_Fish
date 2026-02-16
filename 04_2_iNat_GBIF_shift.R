####
#### R script for Ohigashi et al (2025)
#### Calculation of shift velocity based on iNaturalist and GBIF data
#### 2025.05.21 written by Ohigashi
#### R 4.5.0
####


### load packages
library(dplyr); packageVersion("dplyr")
library(tidyr); packageVersion("tidyr")


### load data
# iNat + GBIF data
combined_data_uni_cl <- readRDS("04_iNat_GBIF_out/combined_data_1960to2024_cleaned.rds")
combined_data_uni_cl <- combined_data_uni_cl %>%
  select(FishBase_name, latitude = decimalLatitude, longitude = decimalLongitude, year, source_group)


###### calculate shift in various areas ######
# set areas
areapatterns <- c("EAsia_Oceania", "outof_EAsia_Oceania", "All")

# create a list to store the result
shift_list <- list()

# loop to calculate the shift velocity for each species
for (area in areapatterns) {
  ## 1. limit area and exclude outliers in the area
  # add absolute latitude column
  filtered_tmp <- combined_data_uni_cl %>%
    mutate(abs_lat = abs(latitude))  %>%
    filter(year >= 0)
  # limit the area by the pattern
  if (area == "EAsia_Oceania") {
    filtered_tmp <- filtered_tmp %>%
      filter(longitude >= 105, longitude <= 180)
  } else if (area == "outof_EAsia_Oceania") {
    filtered_tmp <- filtered_tmp %>%
      filter(longitude > -180, longitude <= 105)
  } else {
    filtered_tmp <- filtered_tmp # All
  }
  
  # exclude outliers within the area
  # poleward (absolute latitude based)
  filtered_data_pw <- filtered_tmp %>%
    filter(!is.na(abs_lat)) %>%
    group_by(FishBase_name, year) %>%
    filter(n() >= 10) %>%
    # filter(n() >= 20) %>% # trial
    filter(between(abs_lat, quantile(abs_lat, 0.05), quantile(abs_lat, 0.95))) %>%
    ungroup() %>%
    group_by(FishBase_name) %>%
    filter(n_distinct(year) >= 5) %>%
    ungroup()
  
  ## 2. calculate effort to be used as weight
  # poleward (absolute latitude based)
  # count number of observation
  latitude_band_effort_pw <- filtered_data_pw %>%
    mutate(lat_band = floor(abs_lat)) %>%
    group_by(year) %>% 
    count(lat_band, name = "band_effort_yr")
  combined_weighted_pw <- filtered_data_pw %>%
    mutate(lat_band = floor(abs_lat)) %>%
    left_join(latitude_band_effort_pw, by = c("lat_band", "year")) %>%
    mutate(weight = 1 / band_effort_yr,
           weight_log = 1 / log1p(band_effort_yr),
           weight_sq = 1/ sqrt(band_effort_yr))
  
  ## 3. calculate range shift velocity
  # poleward shift
  shift_lm_pw <- combined_weighted_pw %>%
    group_by(FishBase_name) %>%
    summarise(
      # normal (unweighted)
      slope = tryCatch(coef(lm(abs_lat ~ year))[["year"]], error = function(e) NA),
      pval = tryCatch(summary(lm(abs_lat ~ year))$coefficients["year", "Pr(>|t|)"], error = function(e) NA),
      r2 = tryCatch(summary(lm(abs_lat ~ year))$r.squared, error = function(e) NA),
      # weighted by 1 / band total observation
      slope.w = tryCatch(coef(lm(abs_lat ~ year, weights = weight))[["year"]], error = function(e) NA),
      pval.w = tryCatch(summary(lm(abs_lat ~ year, weights = weight))$coefficients["year", "Pr(>|t|)"], error = function(e) NA),
      r2.w = tryCatch(summary(lm(abs_lat ~ year, weights = weight))$r.squared, error = function(e) NA),
      # weighted by 1 / log(band total observation)
      slope.wlog = tryCatch(coef(lm(abs_lat ~ year, weights = weight_log))[["year"]], error = function(e) NA),
      pval.wlog = tryCatch(summary(lm(abs_lat ~ year, weights = weight_log))$coefficients["year", "Pr(>|t|)"], error = function(e) NA),
      r2.wlog = tryCatch(summary(lm(abs_lat ~ year, weights = weight_log))$r.squared, error = function(e) NA),
      # weighted by 1 / square-root(band total observation)
      slope.wsq = tryCatch(coef(lm(abs_lat ~ year, weights = weight_sq))[["year"]], error = function(e) NA),
      pval.wsq = tryCatch(summary(lm(abs_lat ~ year, weights = weight_sq))$coefficients["year", "Pr(>|t|)"], error = function(e) NA),
      r2.wsq = tryCatch(summary(lm(abs_lat ~ year, weights = weight_sq))$r.squared, error = function(e) NA)
    )
  
  
  ## 4. Mean latitude where the species observed
  sp_mean_abslat <- filtered_data_pw %>%
    group_by(FishBase_name) %>%
    summarise(mean_abslat = mean(abs_lat))
  sp_mean_lat_hemi <- filtered_data_pw %>%
    mutate(hemisphere = if_else(latitude >= 0, "Northern", "Southern")) %>%
    group_by(FishBase_name, hemisphere) %>%
    summarise(mean_lat_hemi = mean(latitude, na.rm = TRUE), .groups = "drop")
  sp_mean_lat_hemi_wide <- sp_mean_lat_hemi %>%
    pivot_wider(
      names_from = hemisphere,
      values_from = mean_lat_hemi,
      names_prefix = "mean_lat_hemi"
    )
  
  
  shift_lm_pw <- shift_lm_pw %>% 
    left_join(sp_mean_abslat, by = "FishBase_name") %>%
    left_join(sp_mean_lat_hemi_wide, by = "FishBase_name")
  
  
  shift_list[[area]] <- shift_lm_pw
  
}


### save data
saveRDS(shift_list, "04_iNat_GBIF_out/shift_velocity.rds")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/04_2_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))


