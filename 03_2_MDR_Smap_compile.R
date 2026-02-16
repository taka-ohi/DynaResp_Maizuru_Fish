####
#### R script for Ohigashi et al (2025)
#### compile interaction strength calculated in MDR S-map
#### 2025.05.21 written by Ohigashi
#### R 4.5.0
####


### load packages and functions
source("Function/F1_HelperFunction_stats.R")
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")


### load data
# MDR S-map results
load("03_MDR_Smap_out/MDR_Smap.RData")

# fish data
ps_all <- readRDS("01_dataformatting_out/fishdata_ps.rds")


### preparation
# set names of used environmental variables
env_vars <- c("Water_temp_bottom", "Water_temp_bottom.diff")

# extract not "NULL" data from the result lists
bt_mdr_res_surr <- bt_mdr_res_surr_all[!sapply(bt_mdr_res_surr_all, is.null)]
bt_mdr_res_con <- bt_mdr_res_con_all[!sapply(bt_mdr_res_con_all, is.null)]
bt_mdr_res_surr_macro <- bt_mdr_res_surr_macro[!sapply(bt_mdr_res_surr_macro, is.null)]
bt_mdr_res_con_macro <- bt_mdr_res_con_macro[!sapply(bt_mdr_res_con_macro, is.null)]


# make a list of MDR S-map result
mdrres_list <- list("mdr_surr_bt.to.fc" = bt_mdr_res_surr,
                    "mdr_con_btdiff.to.fcdiff" = bt_mdr_res_con,
                    "mdr_surr_bt.to.macro" = bt_mdr_res_surr_macro,
                    "mdr_con_btdiff.to.macrodiff" = bt_mdr_res_con_macro
                    )

# <---------------------------------------------> #
# Summarize S-map coefficients to data.frame
# <---------------------------------------------> #
# Assign interaction strength to each element
# ref: https://github.com/ong8181/eDNA-BosoFish-network/blob/main/06_CompileSmapCoef.R 

# looping by the list
smap_summary_list <- list()

for (l in names(mdrres_list)) {
  lis <- mdrres_list[[l]]
  
  # modify sample data considering number of data (because "diff" does not have the first sample)
  if (grepl("diff", l)) {
    safc_safcdiff_tmp <- safc_safcdiff %>% filter(Sample != "Sample001")
  } else {
    safc_safcdiff_tmp <- safc_safcdiff
  }
  samdata_tmp <- safc_safcdiff_tmp |>
      select(all_of(names(safc_safcdiff_tmp)[!grepl("Fish", names(safc_safcdiff_tmp))]))
  
  for (i in 1:length(lis)) {
    ## Extract an effect variable
    tax_i <- names(lis)[i]
    ## Extract embedding information
    block_vars_i <- lis[[tax_i]]$used_vars %>%
      str_split(pattern = "_tp") %>% sapply(`[`, 1)
    # Remove time-delayed self-interactions
    tax_valid_idx <- c(1, which(block_vars_i != tax_i))
    # Extract tp information
    tp_vars_i <- lis[[tax_i]]$used_vars %>%
      str_split(pattern = "_tp") %>% sapply(`[`, 2) %>% 
      as.numeric
    ## Extract S-map coefficients column names
    smapc_colnames <- paste0(sprintf("c_%s", 1:length(block_vars_i)))
    ## Extract S-map coefficients
    smapc_df <- lis[[tax_i]]$mdr_res$smap_coefficients %>%
      select(all_of(smapc_colnames)) 
    
    ## Make delayed block to check whether a focal species pair is
    block_delay_sps <- block_vars_i[tax_valid_idx][!(block_vars_i[tax_valid_idx] %in% env_vars)] %>% 
      safc_safcdiff_tmp[.]
    block_delay_env <- block_vars_i[tax_valid_idx][(block_vars_i[tax_valid_idx] %in% env_vars)] %>% 
      safc_safcdiff_tmp[.]
    block_delay <- cbind(block_delay_sps, block_delay_env)
    for (block_i in 1:ncol(block_delay)) {
      tp_i <- tp_vars_i[tax_valid_idx][block_i]
      block_delay[block_i] <- dplyr::lag(block_delay[block_i], n = abs(tp_i))
    }
    
    
    ## Create one tidy data.frame that includes all information
    smapc_df_tmp <- smapc_df[,tax_valid_idx]
    colnames(smapc_df_tmp) <- sprintf("effect_from_%s", block_vars_i[tax_valid_idx])
    smapc_df_tmp <- cbind(smapc_df_tmp, samdata_tmp)
    smapc_df_tmp$sample_id <- rownames(safc_safcdiff_tmp)
    smapc_df_tmp$effect_var <- tax_i
    smapc_df_tmp$effect_var_val <- block_delay[,tax_i]
    smapc_df_tmp_long <- pivot_longer(smapc_df_tmp,
                                      cols = -c(c(colnames(samdata_tmp), "sample_id", "effect_var", "effect_var_val")),
                                      names_to = "cause_var", values_to = "IS")
    # Assign causal var's value.
    smapc_df_tmp_long$cause_var_val <- NA
    for (j in 1:nrow(smapc_df_tmp_long)) {
      cause_var_i <- smapc_df_tmp_long$cause_var[j] %>% 
        str_split(pattern = "effect_from_") %>% .[[1]] %>% .[2]
      sample_id_i <- smapc_df_tmp_long$sample_id[j]
      # Assign values to the data.frame
      smapc_df_tmp_long$cause_var_val[j] <- block_delay[sample_id_i, cause_var_i]
    }
    valid_smapc_cond1 <- smapc_df_tmp_long$cause_var_val > 0 & !is.na(smapc_df_tmp_long$cause_var_val)
    valid_smapc_cond2 <- smapc_df_tmp_long$effect_var_val > 0 & !is.na(smapc_df_tmp_long$effect_var_val)
    valid_smapc_cond3 <- !is.na(smapc_df_tmp_long$IS)
    
    if (i == 1) {
      smapc_df_long <- smapc_df_tmp_long[valid_smapc_cond1 & valid_smapc_cond2 & valid_smapc_cond3,] %>%
        as.data.frame
    } else {
      smapc_df_tmp_long <- smapc_df_tmp_long[valid_smapc_cond1 & valid_smapc_cond2 & valid_smapc_cond3,] %>% 
        as.data.frame
      row.names(smapc_df_tmp_long) <- as.character((nrow(smapc_df_long)+1):(nrow(smapc_df_long)+nrow(smapc_df_tmp_long)))
      smapc_df_long <- rbind(smapc_df_long, smapc_df_tmp_long)
    }
  }
  smap_summary_list[[l]] <- smapc_df_long
  
}


### extract mean IS data (bottom T -> Fish count) for later use
summary_df <- smap_summary_list[["mdr_surr_bt.to.fc"]] %>%
  filter(grepl("Water_temp", cause_var)) %>%
  group_by(effect_var) %>%
  summarise(
    mean_IS = mean(IS, na.rm = TRUE),
    sd_IS = sd(IS, na.rm = TRUE)
  )
# remove outlier
Q1 <- quantile(summary_df$mean_IS, 0.25)
Q3 <- quantile(summary_df$mean_IS, 0.75)
IQR <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR
summary_df_wo_out <- summary_df %>%
  filter(mean_IS >= lower_bound & mean_IS <= upper_bound)

# combine with taxonomy data
tax_sheet <- tax_table(ps_all) %>%
  as.data.frame() %>%
  rownames_to_column("Fish_ID")
tax_sheet <- tax_sheet %>%
  mutate(across(
    everything(),
    ~ if (all(grepl("^\\s*-?\\d*\\.?\\d*\\s*$", .x)) & !anyNA(suppressWarnings(as.numeric(.x)))) {
      as.numeric(.x)
    } else {
      .x
    }
  ))

fishIS_data <- summary_df_wo_out %>%
  left_join(tax_sheet %>% select(Fish_ID, FishBase_name, Lat_center),
            by = c("effect_var" = "Fish_ID")) %>% 
  select(Fish_ID = effect_var, FishBase_name, mean_IS, sd_IS, Lat_center)

## add mean abundance data
# calculate mean abundance
fishcount_tbl <- otu_table(ps_all)
if (!taxa_are_rows(ps_all)) {
  fishcount_tbl <- t(fishcount_tbl)
}
fish_aveabun <- as.data.frame(fishcount_tbl) %>%
  rownames_to_column("Fish_ID") %>%
  mutate(mean_abun = rowMeans(select(., -Fish_ID))) %>%
  select(Fish_ID, mean_abun)
# add mean abundance
fishIS_data <- merge(fishIS_data, fish_aveabun, by = "Fish_ID", all.x = TRUE)


# <---------------------------------------------> #
# Save data
# <---------------------------------------------> #

# save S-map coefficient data
saveRDS(smap_summary_list, "03_MDR_Smap_out/MDR_Smap_compiled.rds")

# save mean IS with abundance data
saveRDS(fishIS_data, "03_MDR_Smap_out/meanIS_taxinfo_df.rds")

# save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/03_2_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))


