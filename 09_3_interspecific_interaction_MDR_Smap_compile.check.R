####
#### R script for Ohigashi et al (2025)
#### complile MDR S-map data that consider interspecific interactions
#### 2025.11.25 written by Ohigashi
#### R 4.5.2
####


### load packages
library(tidyverse); packageVersion("tidyverse") # 2.0.0
library(macamts); packageVersion("macamts") # 0.1.4
library(phyloseq); packageVersion("phyloseq")
source("Function/F2_HelperFunction_vizualize.R")


### load data
# MDR S-map result
all_mdr_res <- readRDS("09_interspecific_interaction_out/MDRSmap_results_w.sp.interaction.rds")

# time series data
# fish phyloseq object
ps_all <- readRDS("01_dataformatting_out/fishdata_ps.rds")


### Preparation
# remove effect vars that could not be done S-map
all_mdr_res_ed <- na.omit(all_mdr_res)

## preparation for fish count data
fishcount <- otu_table(ps_all) %>% as.data.frame()
fishcount <- fishcount %>% rownames_to_column("Sample")
# sample data
sampledata <- data.frame(sample_data(ps_all))
sampledata <- sampledata %>% rownames_to_column("Sample")

# combine fish and sample data
sa_fc <- merge(sampledata, fishcount, by = "Sample", sort = F)

# remove unused data
sa_fc <- sa_fc %>% select(-Year, -Timing, -Water_temp_surface, -Month, -Season, -year_div)

# rename to fit the following loop
full_df <- sa_fc
samdata <- sampledata


##### summarize S-map coefficients to data frame ######
# Assign interaction strength to each element
for (i in 1:nrow(all_mdr_res_ed)) {
  ## Extract an effect variable
  tax_i <- all_mdr_res_ed$effect_var[i]
  ## Extract embedding information
  block_vars_i <- all_mdr_res_ed$block_mvd[i][[1]] %>%
    colnames %>% str_split(pattern = "_tp") %>% sapply(`[`, 1)
  # Remove time-delayed self-interactions
  tax_valid_idx <- c(1, which(block_vars_i != tax_i))
  # Extract tp information
  tp_vars_i <- all_mdr_res_ed$block_mvd[i][[1]] %>%
    colnames %>% str_split(pattern = "_tp") %>% sapply(`[`, 2) %>% 
    as.numeric
  ## Extract S-map coefficients column names
  smapc_colnames <- paste0(sprintf("c_%s", 1:length(block_vars_i)))
  ## Extract S-map coefficients
  smapc_df <- all_mdr_res_ed$smap_coefficients[i][[1]] %>%
    select(all_of(smapc_colnames)) #%>% .[valid_idx,] 
  
  ## Make delayed block to check whether a focal species pair is
  block_delay_sps <- block_vars_i[tax_valid_idx] %>% 
    full_df[.]
  block_delay <- block_delay_sps
  for (block_i in 1:ncol(block_delay)) {
    tp_i <- tp_vars_i[tax_valid_idx][block_i]
    block_delay[block_i] <- dplyr::lag(block_delay[block_i], n = abs(tp_i))
  }
  
  ## Create one tidy data.frame that includes all information
  smapc_df_tmp <- smapc_df[,tax_valid_idx]
  colnames(smapc_df_tmp) <- sprintf("effect_from_%s", block_vars_i[tax_valid_idx])
  smapc_df_tmp <- cbind(smapc_df_tmp, samdata)
  smapc_df_tmp$effect_var <- tax_i
  smapc_df_tmp$effect_var_val <- block_delay[, tax_i, drop = TRUE]
  smapc_df_tmp_long <- pivot_longer(smapc_df_tmp,
                                    cols = -c(c(colnames(samdata), "Sample", "effect_var", "effect_var_val")),
                                    names_to = "cause_var", values_to = "IS")
  # Assign causal varriable value.
  smapc_df_tmp_long$cause_var_val <- NA
  for (j in 1:nrow(smapc_df_tmp_long)) {
    cause_var_i <- smapc_df_tmp_long$cause_var[j] %>% 
      str_split(pattern = "effect_from_") %>% .[[1]] %>% .[2]
    idx <- as.numeric(sub("^Sample", "", smapc_df_tmp_long$Sample[j]))
    # Assign values to the data.frame
    smapc_df_tmp_long$cause_var_val[j] <- as.numeric(block_delay[idx, cause_var_i])
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


### extract mean IS data (bottom T -> Fish count) for later use
summary_df <- smapc_df_long %>%
  filter(grepl("Water_temp", cause_var)) %>%
  group_by(effect_var) %>%
  summarise(
    mean_IS = mean(IS, na.rm = TRUE),
    sd_IS = sd(IS, na.rm = TRUE)
  )

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

# IS_w.sp.int <- summary_df_wo_out %>%
IS_w.sp.int <- summary_df %>%
  left_join(tax_sheet %>% select(Fish_ID, FishBase_name, Lat_center),
            by = c("effect_var" = "Fish_ID")) %>% 
  select(Fish_ID = effect_var, FishBase_name, mean_IS, sd_IS, Lat_center)


### compare with IS data that does not account for intersecific interaction
# load IS data
IS <- readRDS("03_MDR_Smap_out/meanIS_taxinfo_df.rds")

# get the data of number of significant windows
sig_counts <- readRDS("07_movingwindow_UIC_out/n_sig_mwUIC_df.rds")

# get the data of optimal E
optE_df <- read.csv("03_MDR_Smap_out/OptE_table_fishcount.csv")

# combine IS data accounting for inter specific interactions and IS data not accounting for them
comb_df <- IS_w.sp.int %>% 
  select(Fish_ID, FishBase_name, mean_IS_w.int = mean_IS) %>%
  left_join(IS %>% select(Fish_ID, FishBase_name, mean_IS),
            by = c("Fish_ID", "FishBase_name")) %>%
  left_join(sig_counts, by = c("Fish_ID" = "effect_var")) %>% # validated by moving-window UIC
  left_join(optE_df, by = c("Fish_ID" = "effect_var"))

# filter by significance
comb_df <- comb_df %>%
  filter(n_window_sig >= 1)

vali_plo <- ggplot(comb_df, aes(x = mean_IS, y = mean_IS_w.int)) +
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  scale_x_continuous(labels = scaleFUN2) +
  scale_y_continuous(labels = scaleFUN2) +
  labs(
    x = "Mean dynamic response to temperature\n(including indirect effects of interspecific interactions)",
    y = "Mean dynamic response to temperature\n(excluding interspecific interactions)"
  )

# relationship between optimal E and difference (mean IS - mean IS excluding interspecific interaction)
comb_df_ed <- comb_df %>%
  mutate(IS_abs_diff = abs(mean_IS - mean_IS_w.int))

diff_E_plo <- ggplot(comb_df_ed, aes(x = OptE, y = IS_abs_diff)) +
  geom_point(size = 1) +
  theme_bw() +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  scale_y_continuous(labels = scaleFUN2) +
  labs(
    x = "Optimal E",
    y = "diff. mean dynamic responses to temperature\n(including \U2212 excluding interspecific interactions)"
  )

### save results
# compiled data frame
saveRDS(smapc_df_long, "09_interspecific_interaction_out/IS_w.sp.interaction_df.rds")

# mean temp-to-sp IS data frame
saveRDS(IS_w.sp.int, "09_interspecific_interaction_out/meanIS_w.sp.interaction_df.rds")

# plot
ggsave("09_interspecific_interaction_out/validation_IS_and_IS.w.species.interaction.png", plot = vali_plo,
       width = 7, height = 7, bg = "white")
saveRDS(vali_plo, "09_interspecific_interaction_out/validation_IS_and_IS.w.species.interaction.rds")

ggsave("09_interspecific_interaction_out/diff_IS_and_IS.w.species.interaction_optE.png", plot = diff_E_plo,
       width = 7, height = 7, bg = "white")
saveRDS(diff_E_plo, "09_interspecific_interaction_out/diff_IS_and_IS.w.species.interaction_optE.rds")


