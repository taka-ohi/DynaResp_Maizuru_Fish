####
#### R script for Ohigashi et al (2025)
#### calculate correlations between variables
#### 2025.06.17 written by Ohigashi
#### R 4.5.0
####


### load packages
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(purrr); packageVersion("purrr")
library(ape); packageVersion("ape")
library(ggplot2); packageVersion("ggplot2")
library(cowplot); packageVersion("cowplot")
source("Function/F2_HelperFunction_vizualize.R")


### load data
# IS and distribution center with moving UIC results
meanIS_taxinfo <- readRDS("07_movingwindow_UIC_out/meanIS_taxinfo_mwUIC_df.rds")

# shift velocity
shift_velocity_list <- readRDS("04_iNat_GBIF_out/shift_velocity.rds")
shift_v_EAsiaOceania <- shift_velocity_list$EAsia_Oceania

# conventional response traits
conv_res <- readRDS("06_conventional_responsetrait_out/lm_gam_coef.rds")

# phylogenetic tree of Maizuru fish
fishtree <- read.tree("05_phylogeny_out/RAxML_tree_w.B.nwk")

# fishdata
tax_sheet <- read.csv("01_dataformatting_out/fishinfo_w_ecology.csv")


### format data
# combine data frames
full_df <- tax_sheet %>%
  select(Fish_ID, FishBase_name, Lat_center) 
full_df <- full_df %>%
  left_join(conv_res, by = "Fish_ID") %>%
  left_join(meanIS_taxinfo %>% select(-Lat_center), by = c("Fish_ID", "FishBase_name")) %>%
  left_join(shift_v_EAsiaOceania, by = "FishBase_name")

# make a list of variable pairs
var_pairs <- list(
  c("Lat_center", "mean_IS"),
  c("Lat_center", "lm_coef"),
  c("Lat_center", "gam_slope_ave"),
  c("mean_IS", "slope.w"),
  c("lm_coef", "slope.w"),
  c("gam_slope_ave", "slope.w")
)


### calculate normal correlations
cor_results <- map_dfr(var_pairs, function(pair) {
  v1 <- pair[1]
  v2 <- pair[2]
  
  # remove undetected species by the moving-window UIC in the case of mean IS
  if (v1 == "mean_IS" | v2 == "mean_IS") {
    sub_df <- full_df %>%
      filter(n_window_sig >= 1) # significant at least one time
  } else {
    sub_df <- full_df
  }
  
  # omit NA and calculate correlations
  test_result <- cor.test(sub_df[[v1]], sub_df[[v2]], use = "complete.obs")
  
  tibble(
    v1 = v1,
    v2 = v2,
    r = test_result$estimate,
    p = test_result$p.value
  )
})


### calculate PIC-based correlations
pic_results <- map_dfr(var_pairs, function(pair) {
  v1 <- pair[1]
  v2 <- pair[2]
  
  # remove undetected species by the moving-window UIC in the case of mean IS
  if (v1 == "mean_IS" | v2 == "mean_IS") {
    sub_df0 <- full_df %>%
      filter(n_window_sig >= 1) # significant at least one time
  } else {
    sub_df0 <- full_df
  }
  
  # select focal data and omit NA
  sub_df <- sub_df0 %>%
    select(FishBase_name, all_of(v1), all_of(v2)) %>%
    na.omit()
  
  pruned_tree <- drop.tip(fishtree, setdiff(fishtree$tip.label, sub_df$FishBase_name))
  # order data frame as the tree is
  sub_df_reo <- sub_df[match(pruned_tree$tip.label,sub_df$FishBase_name),]
  
  # Calculate PICs
  pic_v1 <- pic(sub_df_reo[[v1]], pruned_tree)
  pic_v2 <- pic(sub_df_reo[[v2]], pruned_tree)
  
  
  test_result <- cor.test(pic_v1, pic_v2, use = "complete.obs")
  
  tibble(
    v1 = v1,
    v2 = v2,
    PIC_r = test_result$estimate,
    PIC_p = test_result$p.value
  )
})


### save data
dir.create("08_correlations_out")
saveRDS(full_df, "08_correlations_out/analysis_df.rds")
write.csv(cor_results, "08_correlations_out/correlation_normal.csv", quote = F, row.names = F)
write.csv(pic_results, "08_correlations_out/correlation_PIC.csv", quote = F, row.names = F)


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/08_1_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))




