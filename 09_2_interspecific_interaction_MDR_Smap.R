####
#### R script for Ohigashi et al (2025)
#### MDR S-map considering interspecific interactions
#### 2025.11.25 written by Ohigashi; 2026.01.21 edited by Ohigashi
#### R 4.5.2
####


### load paackeges
source("Function/F1_HelperFunction_stats.R")
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(macamts); packageVersion("macamts")
library(rEDM); packageVersion("rEDM")
library(rUIC); packageVersion("rUIC")
library(pbapply); packageVersion("pbapply")


### load data
# UIC result of interaction effects
uic_res_list <- readRDS("09_interspecific_interaction_out/intersp_effect_list_condv.temp.rds")

# UIC result of temperature to species interaction
uic_surr_bt.to.fc <- readRDS("02_UIC_out/uic_bottomT_w_season.surrogate_fc.rds")

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
# sa_std <- sa_fc_std %>% select(Water_temp_bottom, !contains("Fish"))
fc_std <- sa_fc_std %>% select(Water_temp_bottom, contains("Fish"))

# preparation of UIC result (temp to fish)
uic_surr_bt_fin <- uic_surr_bt.to.fc  %>%
  select(effect_var, cause_var, tp, rmse_surr = rmse, te = te, p_surr = surr_p, pval = surr_p.adj) %>% # "te" and "pval" names are needed
  filter(tp <= 0)


### 2. perform MDR S-map interspecific interaction + temperature effect ###
all_mdr_res_list <- pblapply(
  names(uic_res_list),
  cl = 32,
  FUN = function(ef_var) {
    
    set.seed(123)
    # get data frame of UIC result (species interaction)
    uic_res <- uic_res_list[[ef_var]] %>%
      dplyr::filter(cause_var != "Water_temp_bottom") # remove effect from temp.
    
    # get data frame of UIC result (temp to species)
    uic_res_temp.to.fc <- uic_surr_bt_fin %>% filter(effect_var == ef_var)
    # check max-TE tp
    temp_max.te_tp <- uic_res_temp.to.fc %>%
      arrange(desc(te)) %>%
      slice(1) %>%
      pull(tp)

    # create a lagged series of temperature
    block_temp <- dplyr::lag(fc_std[,"Water_temp_bottom"], n = abs(temp_max.te_tp)) %>%
      as.data.frame() %>%
      setNames(paste0("Water_temp_bottom_tp", temp_max.te_tp))
    
    # Estimate optimal embeding dimension of fish
    simp_x <- rUIC::simplex(fc_std, lib_var = ef_var, E = 0:24, tp = 1)
    (Ex <- simp_x[which.min(simp_x$rmse),"E"])
    
    # try MDR S-map
    res <- tryCatch({
      # Make block to calculate multiview distance
      block_mvd <- tryCatch(
        {
          # when there are causal variables in other species (interspecific interaction)
          block_mvd0 <- make_block_mvd(
            fc_std, uic_res, ef_var,
            max_lag = Ex,
            include_var  = "strongest_only",
            p_threshold  = 0.05
          )
          # combine with lagged temperature (of max TE)
          cbind(block_mvd0, block_temp)
        },
        error = function(e) {
          message(e$message, " Thus only temperature is used for embedding."
                 )
          block_mvd_alt <- make_block_mvd_TOed(
            fc_std, uic_res_temp.to.fc, ef_var,
            max_lag = Ex,
            include_var  = "strongest_only",
            p_threshold = 1
          )
          block_mvd_alt
        }
      )

      # Compute multiview distance
      # determine number of dimension used in MDR
      low_d_E <- ifelse(Ex >= 3, 3, Ex)
      multiview_dist <- compute_mvd(block_mvd,
                                    ef_var,
                                    E = low_d_E,
                                    tp = 1,
                                    distance_only = FALSE)

      # Do MDR S-map
      theta_range <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8)
      lambda_range <- c(0, 0.001, 0.01, 0.1, 0.5, 1, 2, 4, 8)
      least_rmse <- Inf
      smap_theta <- NA_real_
      smap_lambda <- NA_real_
      mdr_res_final <- NULL

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

      # Calculate p-values (if the S-map succeeds)
      if (!is.null(mdr_res_final$model_output)) {
        mdr_res_final$stats$pval <- cor.test(
          mdr_res_final$model_output$obs,
          mdr_res_final$model_output$pred,
          use = "complete.obs"
        )$p.value
      } else {
        mdr_res_final$stats$pval <- NA_real_
      }

      list(ok = TRUE,
           E = Ex,
           block_mvd = block_mvd,
           multiview_dist = multiview_dist,
           mdr = mdr_res_final,
           theta = smap_theta,
           lambda = smap_lambda)

    }, error = function(e) {
      message(sprintf("Error in effect_var '%s': %s", ef_var, e$message))
      list(ok = FALSE,
           E = Ex)
    })
  return(res)
}
)

# name the list as it is
names(all_mdr_res_list) <- names(uic_res_list)



# convert list to data frame
nm <- names(all_mdr_res_list)

all_mdr_res <- do.call(rbind, lapply(nm, function(ef_var_name) {
  res <- all_mdr_res_list[[ef_var_name]]
  
  if (!res$ok) {
    return(data.frame(
      effect_var = ef_var_name, E = res$E,
      N = NA, rho = NA, mae = NA, rmse = NA, pval = NA,
      theta = NA, lambda = NA,
      block_mvd = I(list(NA)),
      mvd_parms = I(list(NA)),
      top_embeddings = I(list(NA)),
      smap_coefficients = I(list(NA))
    ))
  }
  
  df <- data.frame(
    effect_var = ef_var_name,
    E = res$E,
    N = res$mdr$stats$N,
    rho = res$mdr$stats$rho,
    mae = res$mdr$stats$mae,
    rmse = res$mdr$stats$rmse,
    pval = res$mdr$stats$pval,
    theta = res$theta,
    lambda = res$lambda
  )
  
  df$block_mvd         <- I(list(res$block_mvd))
  df$mvd_parms         <- I(list(res$multiview_dist$parms))
  df$top_embeddings    <- I(list(res$multiview_dist$top_embeddings))
  df$smap_coefficients <- I(list(res$mdr$smap_coefficients))
  
  df
}))



### save result
saveRDS(all_mdr_res, "09_interspecific_interaction_out/MDRSmap_results_w.sp.interaction.rds")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/09_2_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))


