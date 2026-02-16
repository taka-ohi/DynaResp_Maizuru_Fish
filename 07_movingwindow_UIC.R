####
#### R script for Ohigashi et al (2025)
#### Analysis of causal effect from temperature to fish using UIC (using moving window)
#### 2025.05.30 written by Ohigashi
#### R 4.5.0
#### eliminate (mitigate) seasonality by using surrogate data in UIC


### load packages and functions
source("Function/F1_HelperFunction_stats.R")
library(rUIC); packageVersion("rUIC") # 0.9.12
library(rEDM); packageVersion("rEDM") # 0.7.5
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(phyloseq); packageVersion("phyloseq")
library(pbapply); packageVersion("pbapply")
library(purrr); packageVersion("purrr")
library(ggplot2); packageVersion("ggplot2")


### load data
load("02_UIC_out/UIC_temp.to.vars_surrogate.RData")


### format data
# exclude macroecology data from looped_vars (use only fish population data)
looped_var_sc.fish <- looped_var_sc %>%
  as.data.frame() %>%
  select(contains("Fish"))
looped_varnames.fish <- colnames(looped_var_sc.fish)


# check the number of dimensions
length(looped_varnames.fish); dim(looped_var_sc.fish); dim(surr_bottomT_sc)
all(looped_varnames.fish == colnames(looped_var_sc.fish))


####### loop for moving windows ######
## whole loop (for 120 time-point windows = 5 years)
# to store the all result
results_all <- list()
# step size
step_size <- 1
# fisrt-start and last-start time points
first_start <- 1
last_start <- 540 - 120 + 1

for (win_start in seq(first_start, last_start, by = step_size)) { # trial with 1 window => 5.5 m (true value is 421 => about 38.5h)
  win_end <- win_start + 120 - 1
  
  # slice the data by window
  temp_obs <- bottomT_sc[win_start:win_end]
  temp_surr <- surr_bottomT_sc[win_start:win_end, ]  # nrow: 120 x 1000
  fish_mat  <- looped_var_sc.fish[win_start:win_end, ]  # nrow: 120 x 113
  
  # for each fish
  window_result <- pblapply(seq_len(ncol(fish_mat)), function(f_idx) {
    fish_id <- colnames(fish_mat)[f_idx]
    fish_ts_win <- fish_mat[, f_idx]
    
    # UIC by observed bottom temperature
    set.seed(123)
    uic_obs <- uic.optimal(
      data.frame(fish_ts_win, temp_obs),
      lib_var = 1, tar_var = 2, E = 0:24, tau = 1, tp = -12:2, num_surr = 1
    )
    uic_obs_fin <- data.frame(effect_var = fish_id,
                              cause_var = "Water_temp_bottom",
                              uic_obs)
    
    # UIC by surrogate bottom temperature
    surr_list <- lapply(seq_len(ncol(temp_surr)), function(i) {
      set.seed(123)
      uic_surr <- uic.optimal(
        data.frame(fish_ts_win, temp_surr[, i]),
        lib_var = 1, tar_var = 2, E = 0:24, tau = 1, tp = -12:2, num_surr = 1
      )
      uic_res.df <- data.frame(effect_var = fish_id,
                               cause_var = colnames(temp_surr)[i],
                               uic_surr)
      uic_res.df
    })
    # convert the list to data frame
    surr_df <- do.call(rbind, surr_list)
    
    # get a combination of parameters (effect varnames x tp)
    params_list <- expand.grid(i = fish_id, j = -12:2)
    
    # get surrogate names
    surr_names_bt <- unique(surr_df$cause_var)
    
    # calculate p value
    results_df <- pmap_dfr(params_list, function(i, j) {
      res <- calculate_surr_p_ci(
        params = list(i, j),
        uic_obs = uic_obs_fin,
        uic_surr = surr_df,
        surr_names = surr_names_bt,
        test_column = "te"
      )
      as_tibble(res)
    })
    
    # store the result
    data.frame(
      window_start = win_start,
      window_end = win_end,
      results_df
    )
  }, cl = 32)
  
  results_all[[as.character(win_start)]] <- do.call(rbind, window_result)
  
  message(sprintf("process %d / %d is done", win_start, last_start))
}


dir.create("07_movingwindow_UIC_out")
saveRDS(results_all, "07_movingwindow_UIC_out/movingwindow_UIC_resultlist.rds")

# bind all result
final_df <- {
  res <- do.call(rbind, results_all)
  rownames(res) <- NULL
  res
}

saveRDS(final_df, "07_movingwindow_UIC_out/movingwindow_UIC_resultdf.rds")


## resume from here
final_df <- readRDS("07_movingwindow_UIC_out/movingwindow_UIC_resultdf.rds")

# calculate adjusted p-values
final_df_adj <- final_df %>%
  filter(tp <= 0) %>% 
  group_by(window_start, tp) %>%
  mutate(p.adj = p.adjust(surr_p, method = "BH")) %>% # adjusted among fish within same window and tp
  ungroup()

# add global adjusted p-values (i.e., considering inter-window adjustment)
final_df_adj <- final_df_adj %>%
  group_by(tp) %>%
  mutate(p.adj_global = p.adjust(surr_p, method = "BH")) %>% # adjusted among fish x window within same tp
  ungroup()

## create heatmap
# consider p < 0.05 values as significant at any tp
heatmap_df <- final_df_adj %>%
  group_by(window_start, effect_var) %>%
  summarise(sig = any(p.adj_global < 0.05), .groups = "drop")

heat <- ggplot(heatmap_df, aes(x = window_start, y = effect_var, fill = sig)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "red"),
                    name = "Including\np.adj < 0.05") +
  theme_minimal() +
  labs(x = "Window start", y = "Fish species",
       title = "UIC-detected windows for each fish (window size = 120)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# save heatmap
ggsave("07_movingwindow_UIC_out/UIC_detected_windows_sp.png", plot = heat,
       width = 10, height = 14, bg = "white")

## summary 
# count the windows where the species had p < 0.05
sig_counts <- heatmap_df %>%
  group_by(effect_var) %>%
  summarise(n_window_sig = sum(sig, na.rm = TRUE)) %>%
  arrange(desc(n_window_sig))

# check IS-quantified species
meanIS_taxinfo <- readRDS("03_MDR_Smap_out/meanIS_taxinfo_df.rds")

uicsig_counts_w_IS <- meanIS_taxinfo %>%
  left_join(sig_counts, by = c("Fish_ID" = "effect_var"))

# save data frame which contains both number of windows with uic-significance and IS
saveRDS(uicsig_counts_w_IS, "07_movingwindow_UIC_out/meanIS_taxinfo_mwUIC_df.rds")

# save data frame of number of windows wiih uic-significance (all species)
saveRDS(sig_counts, "07_movingwindow_UIC_out/n_sig_mwUIC_df.rds")

