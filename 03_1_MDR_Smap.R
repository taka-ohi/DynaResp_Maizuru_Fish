####
#### R script for Ohigashi et al (2025)
#### MDR S-map for the UIC-based causal effects focusing on all species
#### 2025.05.21 written by Ohigashi; 2025.12.23 edited by Ohigashi; 2026.01.21 edited by Ohigashi
#### R 4.5.0
#### 


### load packages and functions
source("Function/F1_HelperFunction_stats.R")
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(macamts); packageVersion("macamts") # 0.2.2
library(rEDM); packageVersion("rEDM")
library(rUIC); packageVersion("rUIC")
library(pbapply); packageVersion("pbapply")


### load data
# UIC result by surrogate tests (Fish count)
uic_surr_bt.to.fc <- readRDS("02_UIC_out/uic_bottomT_w_season.surrogate_fc.rds")
uic_surr_bt.to.macro <- readRDS("02_UIC_out/uic_bottomT_w_season.surrogate_macro.rds")

# UIC result by conditioning (Fish count)
uic_con_bt.diff.to.fc.diff <- readRDS("02_UIC_out/uic_bottomT.diff_to_fc.diff_cond.temp.rds")
uic_con_bt.diff.to.macro.diff <- readRDS("02_UIC_out/uic_bottomT.diff_to_macro.diff_cond.temp.rds")

# time series data
# fish phyloseq object
ps_all <- readRDS("01_dataformatting_out/fishdata_ps.rds")


### 1. preparation
## 1-1. preparation for fish count data
fishcount <- otu_table(ps_all) %>% as.data.frame()
fishcount <- fishcount %>% rownames_to_column("Sample")
# sample data
sampledata <- data.frame(sample_data(ps_all))
sampledata <- sampledata %>% rownames_to_column("Sample")

# combine fish and sample data
sa_fc <- merge(sampledata, fishcount, by = "Sample", sort = F)

# remove unused data
sa_fc <- sa_fc %>% select(-Year, -Timing, -Water_temp_surface, -Month, -Season, -year_div)

# scale
sa_fc_std <- as.data.frame(apply(sa_fc[,-1], 2, function(x) as.numeric(scale(x))))

# divide fishcount and sample data (both should have water temp)
sa_std <- sa_fc_std %>% select(Water_temp_bottom, !contains("Fish"))
fc_std <- sa_fc_std %>% select(Water_temp_bottom, contains("Fish"))

## make scaled difference data
# take difference
sa_fc_diff <- sa_fc |>
  mutate(across(2:ncol(sa_fc),
                .fns = list(diff = ~ . -lag(.)),
                .names = "{col}.diff")
  ) %>%
  select(Sample, contains(".diff")) %>% # leave only difference data
  na.omit() # remove the first time point

# scale
sa_fc_diff_std <- as.data.frame(apply(sa_fc_diff[,-1], 2, function(x) as.numeric(scale(x))))

# divide fishcount diff and sample diff data (both should have water temp)
sa_diff_std <- sa_fc_diff_std %>% select(Water_temp_bottom.diff, !contains("Fish"))
fc_diff_std <- sa_fc_diff_std %>% select(Water_temp_bottom.diff, contains("Fish"))


## 1-2. preparation for UIC result data
# extract data
## surrogate result
# fish count
uic_surr_bt_fin <- uic_surr_bt.to.fc  %>%
  select(effect_var, cause_var, tp, rmse_surr = rmse, te = te, p_surr = surr_p, pval = surr_p.adj) %>% # "te" and "pval" names are needed
  filter(tp <= 0)
# macro index
uic_surr_bt.to.macro_fin <- uic_surr_bt.to.macro %>%
  select(effect_var, cause_var, tp, rmse_surr = rmse, te = te, pval = surr_p)

## conditioned result 
# fish count difference
uic_con_bt_fin <- uic_con_bt.diff.to.fc.diff |>
  select(effect_var, cause_var, tp, rmse_con = rmse, te = te, p_con = pval, pval = p.adj) %>%  # "te" and "pval" names are needed
  filter(tp <= 0)
# macro index
uic_con_bt.to.macro_fin <- uic_con_bt.diff.to.macro.diff %>%
  select(effect_var, cause_var, tp, rmse_con = rmse, te = te, val = pval)


### 2. MDR S-map
# set rng.seed for reproduction
set.seed(123)

# ------------------------------------------- #
# MDR S-map (bottom T -> species count) 
# ------------------------------------------- #
### calculate optimal E for each fish
# set looped varnames (= effect var)
effect_varnames <- unique(uic_surr_bt_fin$effect_var)

# get optimal E table for each fish
optE_res <- pblapply(effect_varnames, function(varname) {
  simp_res <- rUIC::simplex(fc_std, lib_var = varname, E = 0:24, tp = 1, alpha = 0.05)
  optE <- simp_res[which.min(simp_res$rmse), "E"]
  return(data.frame(effect_var = varname, OptE = optE))
},
cl = 32)
optE_tab_fc <- do.call(rbind, optE_res)

### loop for making a list of MDR S-map results
bt_mdr_res_surr_all <- pblapply(effect_varnames, function(efvar) {
  # initialize result
  result <- NULL
  
  ## 0. get targeted dataframe
  # extract time series
  sub_ts <- fc_std |> select(Water_temp_bottom, all_of(efvar)) # check data frame & cause var name
  # extract UIC result
  sub_uicres <- uic_surr_bt_fin |> filter(effect_var == efvar) # check data frame
  # set optimal E
  obs_optE <- optE_tab_fc$OptE[optE_tab_fc$effect_var == efvar]
  
  tryCatch({
    ## 1. create block
    # updated version for macamts package
    block_mvd <- make_block_mvd_TOed(block = sub_ts,
                                     uic_res = sub_uicres,
                                     effect_var = efvar,
                                     max_lag = obs_optE,
                                     cause_var_colname = "cause_var",
                                     include_var = "strongest_only",
                                     p_threshold = 1 # check
    )
    
    ## 2. calculate multiview distance
    # determine number of dimension used in MDR
    low_d_E <- ifelse(obs_optE >= 3, 3, obs_optE)
    multiview_dist <- compute_mvd(block_mvd,
                                  efvar,
                                  E = low_d_E,
                                  tp = 1,
                                  distance_only = FALSE)
    
    ## 3. MDR S-map
    theta_range <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8)
    lambda_range <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8)
    least_rmse <- Inf
    
    for (th in theta_range) {
      for (lamb in lambda_range) {
        # temporal MDR S-map
        mdr_res.tmp <- s_map_mdr(block_mvd,
                                 dist_w = multiview_dist$multiview_dist,
                                 theta = th,
                                 tp = 1,
                                 regularized = TRUE,
                                 alpha = 0,
                                 lambda = lamb,
                                 save_smap_coefficients = TRUE)
        
        # temporal RMSE
        rmse.tmp <- mdr_res.tmp$stats$rmse
        
        # exploration of optimal model
        if (least_rmse > rmse.tmp) {
          least_rmse <- rmse.tmp
          mdr_res_final <- mdr_res.tmp
          smap_theta <- th
          smap_lambda <- lamb
        } 
      }
    }
    
    ## 4. save result
    used_vars <- colnames(block_mvd)
    result <- list(obs_optE = obs_optE, used_vars = used_vars,
                   smap_theta = smap_theta, smap_lambda = smap_lambda,
                   mdr_res = mdr_res_final, multiview_dist_res = multiview_dist)
    
  }, error = function(e) {
    # if error
    # message(sprintf("pvals over threshold for effect variable '%s': %s", efvar, e$message))
    print(sprintf("Error in effect variable '%s': %s", efvar, e$message))
    result <- list(obs_optE = NA, used_vars = NA, smap_theta = NA, 
                   smap_lambda = NA, mdr_res = NA, multiview_dist_res = NA)
  })
  
  # return result
  return(result)
}, cl = 32)

# summarize the result as a list
names(bt_mdr_res_surr_all) <- effect_varnames

# save result
dir.create("03_MDR_Smap_out")
saveRDS(bt_mdr_res_surr_all, "03_MDR_Smap_out/MDR_Smap_res_bottomT.to.fc.rds")
write.csv(optE_tab_fc, "03_MDR_Smap_out/OptE_table_fishcount.csv", quote = F, row.names = F)


# ---------------------------------------------- #
# MDR S-map (bottom T diff -> species count diff) 
# ---------------------------------------------- #
### calculate optimal E for each fish count difference
# set looped varnames (= effect var)
effect_varnames <- unique(uic_con_bt_fin$effect_var)

# get optimal E table for each fish
optE_res <- pblapply(effect_varnames, function(varname) {
  simp_res <- rUIC::simplex(fc_diff_std, lib_var = varname, E = 0:24, tp = 1, alpha = 0.05)
  optE <- simp_res[which.min(simp_res$rmse), "E"]
  return(data.frame(effect_var = varname, OptE = optE))
},
cl = 32)
optE_tab_fc.diff <- do.call(rbind, optE_res)

### loop for making a list of MDR S-map results
bt_mdr_res_con_all <- pblapply(effect_varnames, function(efvar) {
  # initialize result
  result <- NULL
  
  ## 0. get targeted dataframe
  # extract time series
  sub_ts <- fc_diff_std |> select(Water_temp_bottom.diff, all_of(efvar)) # check data frame & cause var name
  # extract UIC result
  sub_uicres <- uic_con_bt_fin |> filter(effect_var == efvar) # check data frame
  # set optimal E
  obs_optE <- optE_tab_fc.diff$OptE[optE_tab_fc.diff$effect_var == efvar]
  
  tryCatch({
    ## 1. create block
    # updated version for macamts package
    block_mvd <- make_block_mvd_TOed(block = sub_ts,
                                     uic_res = sub_uicres,
                                     effect_var = efvar,
                                     max_lag = obs_optE,
                                     cause_var_colname = "cause_var",
                                     include_var = "strongest_only",
                                     p_threshold = 1 # check
    )
    
    ## 2. calculate multiview distance
    # determine number of dimension used in MDR
    low_d_E <- ifelse(obs_optE >= 3, 3, obs_optE)
    multiview_dist <- compute_mvd(block_mvd,
                                  efvar,
                                  E = low_d_E,
                                  tp = 1,
                                  distance_only = FALSE)
    ## 3. MDR S-map
    theta_range <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8)
    lambda_range <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8)
    least_rmse <- Inf
    
    for (th in theta_range) {
      for (lamb in lambda_range) {
        # temporal MDR S-map
        mdr_res.tmp <- s_map_mdr(block_mvd,
                                 dist_w = multiview_dist$multiview_dist,
                                 theta = th,
                                 tp = 1,
                                 regularized = TRUE,
                                 alpha = 0,
                                 lambda = lamb,
                                 save_smap_coefficients = TRUE)
        
        # temporal RMSE
        rmse.tmp <- mdr_res.tmp$stats$rmse
        
        # exploration of optimal model
        if (least_rmse > rmse.tmp) {
          least_rmse <- rmse.tmp
          mdr_res_final <- mdr_res.tmp
          smap_theta <- th
          smap_lambda <- lamb
        } 
      }
    }
    
    ## 4. save result
    used_vars <- colnames(block_mvd)
    result <- list(obs_optE = obs_optE, used_vars = used_vars,
                   smap_theta = smap_theta, smap_lambda = smap_lambda,
                   mdr_res = mdr_res_final, multiview_dist_res = multiview_dist)
    
  }, error = function(e) {
    # if error
    # message(sprintf("pvals over threshold for effect variable '%s': %s", efvar, e$message))
    print(sprintf("Error in effect variable '%s': %s", efvar, e$message))
    result <- list(obs_optE = NA, used_vars = NA, smap_theta = NA, 
                   smap_lambda = NA, mdr_res = NA, multiview_dist_res = NA)
  })
  
  # return result
  return(result)
}, cl = 32)

# summarize the result as a list
names(bt_mdr_res_con_all) <- effect_varnames

# save result
saveRDS(bt_mdr_res_con_all, "03_MDR_Smap_out/MDR_Smap_res_bottomTdiff.to.fcdiff.rds")
write.csv(optE_tab_fc.diff, "03_MDR_Smap_out/OptE_table_fishcount.diff.csv", quote = F, row.names = F)



# ------------------------------------------- #
# MDR S-map (bottom T -> macro) 
# ------------------------------------------- #
### calculate optimal E for each macro index
# set looped varnames (= effect var)
effect_varnames <- unique(uic_surr_bt.to.macro_fin$effect_var)

# get optimal E table for each fish
optE_res <- pblapply(effect_varnames, function(varname) {
  simp_res <- rUIC::simplex(sa_std, lib_var = varname, E = 0:24, tp = 1, alpha = 0.05)
  optE <- simp_res[which.min(simp_res$rmse), "E"]
  return(data.frame(effect_var = varname, OptE = optE))
},
cl = 32)
optE_tab_sa <- do.call(rbind, optE_res)

### loop for making a list of MDR S-map results
bt_mdr_res_surr_macro <- pblapply(effect_varnames, function(efvar) {
  # initialize result
  result <- NULL
  
  ## 0. get targeted dataframe
  # extract time series
  sub_ts <- sa_std |> select(Water_temp_bottom, all_of(efvar)) # check data frame & cause var name
  # extract UIC result
  sub_uicres <- uic_surr_bt.to.macro_fin |> filter(effect_var == efvar) # check data frame
  # set optimal E
  obs_optE <- optE_tab_sa$OptE[optE_tab_sa$effect_var == efvar]
  
  tryCatch({
    ## 1. create block
    # updated version for macamts package
    block_mvd <- make_block_mvd_TOed(block = sub_ts,
                                     uic_res = sub_uicres,
                                     effect_var = efvar,
                                     max_lag = obs_optE,
                                     cause_var_colname = "cause_var",
                                     include_var = "strongest_only",
                                     p_threshold = 1 # check
    )
    
    ## 2. calculate multiview distance
    # determine number of dimension used in MDR
    low_d_E <- ifelse(obs_optE >= 3, 3, obs_optE)
    multiview_dist <- compute_mvd(block_mvd,
                                  efvar,
                                  E = low_d_E,
                                  tp = 1,
                                  distance_only = FALSE)
    
    ## 3. MDR S-map
    theta_range <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8)
    lambda_range <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8)
    least_rmse <- Inf
    
    for (th in theta_range) {
      for (lamb in lambda_range) {
        # temporal MDR S-map
        mdr_res.tmp <- s_map_mdr(block_mvd,
                                 dist_w = multiview_dist$multiview_dist,
                                 theta = th,
                                 tp = 1,
                                 regularized = TRUE,
                                 alpha = 0,
                                 lambda = lamb,
                                 save_smap_coefficients = TRUE)
        
        # temporal RMSE
        rmse.tmp <- mdr_res.tmp$stats$rmse
        
        # exploration of optimal model
        if (least_rmse > rmse.tmp) {
          least_rmse <- rmse.tmp
          mdr_res_final <- mdr_res.tmp
          smap_theta <- th
          smap_lambda <- lamb
        } 
      }
    }
    
    ## 4. save result
    used_vars <- colnames(block_mvd)
    result <- list(obs_optE = obs_optE, used_vars = used_vars,
                   smap_theta = smap_theta, smap_lambda = smap_lambda,
                   mdr_res = mdr_res_final, multiview_dist_res = multiview_dist)
    
  }, error = function(e) {
    # if error
    # message(sprintf("pvals over threshold for effect variable '%s': %s", efvar, e$message))
    print(sprintf("Error in effect variable '%s': %s", efvar, e$message))
    result <- list(obs_optE = NA, used_vars = NA, smap_theta = NA, 
                   smap_lambda = NA, mdr_res = NA, multiview_dist_res = NA)
  })
  
  # return result
  return(result)
}, cl = 32)

# summarize the result as a list
names(bt_mdr_res_surr_macro) <- effect_varnames

# save result
saveRDS(bt_mdr_res_surr_macro, "03_MDR_Smap_out/MDR_Smap_res_bottomT.to.macro.rds")
write.csv(optE_tab_sa, "03_MDR_Smap_out/OptE_table_macro.csv", quote = F, row.names = F)



# ------------------------------------------- #
# MDR S-map (bottom T diff -> macro diff) 
# ------------------------------------------- #
### calculate optimal E for each macro index difference
# set looped varnames (= effect var)
effect_varnames <- unique(uic_con_bt.diff.to.macro.diff$effect_var)

# get optimal E table for each fish
optE_res <- pblapply(effect_varnames, function(varname) {
  simp_res <- rUIC::simplex(sa_diff_std, lib_var = varname, E = 0:24, tp = 1, alpha = 0.05)
  optE <- simp_res[which.min(simp_res$rmse), "E"]
  return(data.frame(effect_var = varname, OptE = optE))
},
cl = 32)
optE_tab_sa.diff <- do.call(rbind, optE_res)

### loop for making a list of MDR S-map results
bt_mdr_res_con_macro <- pblapply(effect_varnames, function(efvar) {
  # initialize result
  result <- NULL
  
  ## 0. get targeted dataframe
  # extract time series
  sub_ts <- sa_diff_std |> select(Water_temp_bottom.diff, all_of(efvar)) # check data frame & cause var name
  # extract UIC result
  sub_uicres <- uic_con_bt.diff.to.macro.diff |> filter(effect_var == efvar) # check data frame
  # set optimal E
  obs_optE <- optE_tab_sa.diff$OptE[optE_tab_sa.diff$effect_var == efvar]
  
  tryCatch({
    ## 1. create block
    # updated version for macamts package
    block_mvd <- make_block_mvd_TOed(block = sub_ts,
                                     uic_res = sub_uicres,
                                     effect_var = efvar,
                                     max_lag = obs_optE,
                                     cause_var_colname = "cause_var",
                                     include_var = "strongest_only",
                                     p_threshold = 1 # check
    )
    
    ## 2. calculate multiview distance
    # determine number of dimension used in MDR
    low_d_E <- ifelse(obs_optE >= 3, 3, obs_optE)
    multiview_dist <- compute_mvd(block_mvd,
                                  efvar,
                                  E = low_d_E,
                                  tp = 1,
                                  distance_only = FALSE)
    
    ## 3. MDR S-map
    theta_range <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8)
    lambda_range <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8)
    least_rmse <- Inf
    
    for (th in theta_range) {
      for (lamb in lambda_range) {
        # temporal MDR S-map
        mdr_res.tmp <- s_map_mdr(block_mvd,
                                 dist_w = multiview_dist$multiview_dist,
                                 theta = th,
                                 tp = 1,
                                 regularized = TRUE,
                                 alpha = 0,
                                 lambda = lamb,
                                 save_smap_coefficients = TRUE)
        
        # temporal RMSE
        rmse.tmp <- mdr_res.tmp$stats$rmse
        
        # exploration of optimal model
        if (least_rmse > rmse.tmp) {
          least_rmse <- rmse.tmp
          mdr_res_final <- mdr_res.tmp
          smap_theta <- th
          smap_lambda <- lamb
        } 
      }
    }
    
    ## 4. save result
    used_vars <- colnames(block_mvd)
    result <- list(obs_optE = obs_optE, used_vars = used_vars,
                   smap_theta = smap_theta, smap_lambda = smap_lambda,
                   mdr_res = mdr_res_final, multiview_dist_res = multiview_dist)
    
  }, error = function(e) {
    # if error
    # message(sprintf("pvals over threshold for effect variable '%s': %s", efvar, e$message))
    print(sprintf("Error in effect variable '%s': %s", efvar, e$message))
    result <- list(obs_optE = NA, used_vars = NA, smap_theta = NA, 
                   smap_lambda = NA, mdr_res = NA, multiview_dist_res = NA)
  })
  
  # return result
  return(result)
}, cl = 32)

# summarize the result as a list
names(bt_mdr_res_con_macro) <- effect_varnames

# save result
saveRDS(bt_mdr_res_con_macro, "03_MDR_Smap_out/MDR_Smap_res_bottomTdiff.to.macrodiff.rds")
write.csv(optE_tab_sa.diff, "03_MDR_Smap_out/OptE_table_macro.diff.csv", quote = F, row.names = F)

## save main results at once
safc_safcdiff <- sa_fc %>%
  left_join(sa_fc_diff, by = "Sample", keep = FALSE)

save(bt_mdr_res_surr_all, optE_tab_fc, bt_mdr_res_con_all, optE_tab_fc.diff,
     bt_mdr_res_surr_macro, optE_tab_sa, bt_mdr_res_con_macro, optE_tab_sa.diff,
     fc_std, fc_diff_std, sa_std, sa_diff_std, safc_safcdiff,
     uic_surr_bt_fin, uic_surr_bt.to.macro_fin, 
     uic_con_bt_fin, uic_con_bt.to.macro_fin,
     file = "03_MDR_Smap_out/MDR_Smap.RData")

### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/03_1_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))
