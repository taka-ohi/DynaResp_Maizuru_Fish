####
#### R script for Ohigashi et al (2025)
#### Obtain data from iNaturalist and GBIF
#### 2025.04.14 written by Ohigashi
#### R 4.5.0
####


### load packages
library(rinat); packageVersion("rinat")
library(rgbif); packageVersion("rgbif")
library(CoordinateCleaner); packageVersion("CoordinateCleaner")
library(countrycode); packageVersion("countrycode")
library(dplyr); packageVersion("dplyr")
library(purrr); packageVersion("purrr")
library(stringr); packageVersion("stringr")


### load data
# fish information
tax_sheet <- read.csv("01_dataformatting_out/fishinfo_w_ecology.csv", header = T)


### 1. Preparation
# make a vector of Maizuru species
fish_mz <- tax_sheet$FishBase_name

# set words that the rows will be excluded
exclude_keywords <- c(
  # English
  "aquarium", "zoo", "market",
  # Japanese
  "水族館", "動物園", "市場",
  # Spanish
  "acuario", "zoológico", "mercado",
  # Chinese
  "水族馆", "动物园", "市场",
  # Korean
  "수족관", "동물원", "시장",
  # French
  "aquarium", "zoo", "marché"
)


### 2. search Maizuru fishes on iNaturalist
# create a function to obtain iNaturalist data
safe_get_inat <- function(species_name, years = 1960:2024) {
  all_obs <- list()
  
  for (yr in years) {
    Sys.sleep(runif(1, 1, 2))  # polite delay
    
    obs.tmp <- tryCatch({
      get_inat_obs(query = gsub("_", " ", species_name),
                   year = yr,
                   maxresults = 10000)
    }, error = function(e) {
      message("Error for ", species_name, " in ", yr, ": ", e$message)
      return(NULL)
    })
    
    if (!is.null(obs.tmp) && nrow(obs.tmp) > 0) {
      obs.tmp$year <- as.numeric(substr(obs.tmp$observed_on, 1, 4))
      
      obs.tmp <- obs.tmp %>%
        select(id, observed_on, year, latitude, longitude, place_guess, description, species_guess, iconic_taxon_name) %>%
        mutate(across(c(place_guess, description, species_guess, iconic_taxon_name), as.character))
      
      
      # exclude trash data
      obs.tmp <- obs.tmp %>%
        filter(
          !is.na(latitude), !is.na(longitude), !is.na(year),
          iconic_taxon_name == "Actinopterygii",
          !str_detect(tolower(coalesce(place_guess, "")), str_c(exclude_keywords, collapse = "|")),
          !str_detect(tolower(coalesce(description, "")), str_c(exclude_keywords, collapse = "|"))
        )
      
      all_obs[[as.character(yr)]] <- obs.tmp
    }
  }
  
  # combine all years
  if (length(all_obs) > 0) {
    bind_rows(all_obs)
  } else {
    NULL
  }
}

# make a list for species data
iNat_list <- list()

# create a directory to save data
dir.create("04_iNat_GBIF_out")

# # (if you resume from the middle part)
# iNat_list <- readRDS("04_iNat_GBIF_out/iNat_partial_results.rds")
# fish_mz_done <- names(iNat_list)
# fish_mz_undone <- setdiff(fish_mz, fish_mz_done)

for (fish in fish_mz) { # change fish_mz to fish_mz_undone if you resume from the middle part #
  message("Processing: ", fish)
  result <- safe_get_inat(species_name = fish, years = 1960:2024)
  iNat_list[[fish]] <- result
  # save the temporal list
  saveRDS(iNat_list, file = "04_iNat_GBIF_out/iNat_partial_results.rds")
}
# save
saveRDS(iNat_list, file = "04_iNat_GBIF_out/iNat_raw_list_1960to2024.rds")


### 3. search Maizuru fishes on GBIF
# create a function to obtain GBIF data
safe_get_gbif <- function(species_name, years = 1960:2024) {
  all_obs <- list()
  for (yr in years) {
    Sys.sleep(runif(1, 1, 2))  # wait for 1 or 2 sec
    obs.tmp <- tryCatch({
      dat <- occ_search(scientificName = gsub("_", " ", species_name),
                        limit = 10000,
                        year = yr,
                        hasCoordinate = TRUE)
      dat$data
    }, error = function(e) {
      message("Error for ", species_name, ": ", e$message)
      return(NULL)
    })
    
    if (!is.null(obs.tmp) && nrow(obs.tmp) > 0) {
      if ("occurrenceRemarks" %in% colnames(obs.tmp)) {
        obs.tmp <- obs.tmp %>%
          filter(!str_detect(tolower(coalesce(occurrenceRemarks, "")), str_c(exclude_keywords, collapse = "|")),
                 basisOfRecord %in% c("HUMAN_OBSERVATION", "MACHINE_OBSERVATION", "OCCURRENCE", "OBSERVATION")
          )
      } else {
        obs.tmp <- obs.tmp %>%
          filter(basisOfRecord %in% c("HUMAN_OBSERVATION", "MACHINE_OBSERVATION", "OCCURRENCE", "OBSERVATION"))
      }
      
      all_obs[[as.character(yr)]] <- obs.tmp
    }
  }
  # combine all years
  if (length(all_obs) > 0) {
    bind_rows(all_obs)
  } else {
    NULL
  }
}

# make a list for species data
GBIF_list <- list()
# # if you resume from the middle part
# GBIF_list <- readRDS("04_iNat_GBIF_out/GBIF_partial_results.rds")
# fish_mz_undone <- setdiff(fish_mz, names(GBIF_list))

for (fish in fish_mz_undone) { # change fish_mz to fish_mz_undone if you resume from the middle part #
  message("Processing: ", fish)
  result <- safe_get_gbif(species_name = fish, years = 1960:2024)
  GBIF_list[[fish]] <- result
  # save the temporal list
  saveRDS(GBIF_list, file = "04_iNat_GBIF_out/GBIF_partial_results.rds")
}
# save
saveRDS(GBIF_list, file = "04_iNat_GBIF_out/GBIF_raw_list_1960to2024.rds")


### 4. combine iNat and GBIF data, and clean it
# read lists
iNat_list <- readRDS("04_iNat_GBIF_out/iNat_raw_list_1960to2024.rds")
GBIF_list <- readRDS("04_iNat_GBIF_out/GBIF_raw_list_1960to2024.rds")

## convert list to data frame
# iNaturalist
inat_data <- iNat_list %>%
  bind_rows(.id = "FishBase_name") %>%
  select(FishBase_name, decimalLatitude=latitude, decimalLongitude=longitude,
         year, IDofSource=id) %>%
  mutate(source = "iNaturalist", IDofSource=as.character(IDofSource)) # add source column


# GBIF
gbif_data <- GBIF_list %>%
  bind_rows(.id = "FishBase_name") %>%
  select(FishBase_name, decimalLatitude,
         decimalLongitude, year,
         IDofSource=gbifID#,
         # countryCode
  ) %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  mutate(source = "GBIF", IDofSource=as.character(IDofSource)) # add source column

## combine iNat and GBIF data
combined_data <- bind_rows(inat_data, gbif_data)

combined_data_clean <- combined_data %>%
  mutate(
    decimalLatitude = round(decimalLatitude, 5),
    decimalLongitude = round(decimalLongitude, 5),
    year = as.integer(year),
    source = trimws(source)
  )
combined_data_unique <- combined_data_clean %>%
  group_by(FishBase_name, decimalLatitude, decimalLongitude, year) %>%
  summarise(source_group = paste(sort(unique(source)), collapse = "+"), .groups = "drop") %>%
  arrange(FishBase_name, decimalLatitude, decimalLongitude)

flags <- clean_coordinates(x = combined_data_unique,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           # countries = "countryCode",
                           species = "FishBase_name",
                           tests = c("capitals", "centroids", "duplicates",
                                     "equal", "zeros", "institutions", "gbif"
                           ))
# Exclude problematic records
combined_data_uni_cl <- combined_data_unique[flags$.summary,]

# save the clean result
saveRDS(combined_data_uni_cl, "04_iNat_GBIF_out/combined_data_1960to2024_cleaned.rds")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/04_1_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))



