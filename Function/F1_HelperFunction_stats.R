### required packages
library(dplyr); packageVersion("dplyr")
library(vegan); packageVersion("vegan")

### function to calculate p values using surrogate data in parallel
calculate_surr_p_ci <- function(params, uic_obs, uic_surr, surr_names, test_column = "te") {
  i <- params[[1]]
  j <- params[[2]]
  
  # get observed TE (or effective TE)
  obs.te <- uic_obs |>
    filter(effect_var == i & tp == j) |>
    pull(all_of(test_column)) # changed to use TE following rUIC
  
  # get a vector for surrogate TEs (or effective TEs)
  surr.te_vec <- c()
  for (k in surr_names) {
    surr.te <- uic_surr |>
      filter(effect_var == i & tp == j & cause_var == k) |>
      pull(all_of(test_column)) # changed to use TE following rUIC
    surr.te_vec <- c(surr.te_vec, surr.te)
  }
  
  # calculate p values
  surr_p <- mean(surr.te_vec >= obs.te)
  
  # calculate lower 95% confidence interval
  ci_upper <- quantile(surr.te_vec, probs = 0.95)
  
  # return the result
  list(effect_var = i, tp = j, surr_p = surr_p, ci_upper = ci_upper)
}


### function to fix the bug in make_block_mvd() in macamts package
# updated on 2026.01.21 due to macamts has updated
make_block_mvd_TOed <- function (block,
                                 uic_res,
                                 effect_var,
                                 # E_effect_var,
                                 max_lag = 1,
                                 cause_var_colname = "cause_var",
                                 include_var = "strongest_only",
                                 p_threshold = 0.050,
                                 sort_tp = TRUE,
                                 silent = FALSE) {
  # Retrieve colnames
  x_names <- colnames(block)
  
  # Check colnames to be used (and critical) in the analysis
  if (is.numeric(effect_var)) effect_var <- x_names[effect_var]
  if (!(effect_var %in% x_names)) stop("No 'effect_var' in the block!")
  if (!(cause_var_colname %in% colnames(uic_res))) stop("No 'cause_var_colname' in `uic_res`! Please specify the correct colname for causal variables")
  if (!("tp" %in% colnames(uic_res))) stop("'tp' column is required for `uic_res`.")
  if (!("te" %in% colnames(uic_res))) stop("'te' column is required for `uic_res`.")
  if (!("pval" %in% colnames(uic_res))) stop("'pval' column is required for `uic_res`.")
  
  # Check arguments
  if (!all(unique(x_names) == x_names)) stop("\"block\" should have unique column names.")
  if (!is.data.frame(block)) stop("\"block\" should be data.frame.")
  if (!is.numeric(p_threshold)) stop("\"p_threshold\" should be numeric.")
  if (!include_var %in% c("all_significant", "strongest_only", "tp0_only")) stop("\"include_var\" should be \"all_significant\", \"strongest_only\", or \"tp0_only\".")
  # if (E_effect_var < 1) stop("\"max_delay_self\" should be >= 1.")
  if (max_lag < 1) stop("\"max_lag\" should be >= 1.")
  
  # Preparation
  if (include_var == "tp0_only") {
    block_mvd <- data.frame(block[,effect_var])
    colnames(block_mvd) <- sprintf("%s_tp0", effect_var)
  } else {
    # block_mvd <- data.frame(rEDM::make_block(block[,effect_var], max_lag = E_effect_var)[,-1])
    # colnames(block_mvd) <- sprintf("%s_tp%s", effect_var, 0:(-(E_effect_var-1)))
    block_mvd <- data.frame(make_block(block[,effect_var], max_lag = max_lag)[,-1])
    colnames(block_mvd) <- sprintf("%s_tp%s", effect_var, 0:(-(max_lag-1)))
  }
  # Pre-screening (p & tp)
  if (!silent) message(sprintf("UIC results with `tp` <= 0 and `pval` <= %s are kept for further analyses.", p_threshold))
  uic_res <- uic_res[uic_res$pval <= p_threshold & uic_res$tp <= 0,]
  if (nrow(uic_res) < 1) stop("No significant causal variables were detected. Please use the univariate S-map.")
  
  # Sort block columns according to tp for each causal variable
  if (sort_tp) {
    sort_id <- 0
    for (col_i in unique(uic_res[,cause_var_colname])) {
      # sort_id_i <- order(-uic_res[uic_res[,cause_var_colname] == col_i,"tp"])
      sort_id_i <- order(-uic_res[uic_res[,cause_var_colname] == col_i,"tp", drop = TRUE]) # edited by Ohigashi
      sort_id <- c(sort_id, max(sort_id) + sort_id_i)
    }
    sort_id <- sort_id[-1]; uic_res <- uic_res[sort_id,]
  }
  
  # ---------------------------------------------------- #
  # Select variables to be included in block
  # ---------------------------------------------------- #
  if (include_var == "all_significant") {
    # ---------------------------------------------------- #
    # Select all significant variables
    # ---------------------------------------------------- #
    for (i in 1:nrow(uic_res)) {
      block_new <- dplyr::lag(block[,uic_res[i, cause_var_colname]], n = abs(uic_res[i,"tp"]))
      block_new <- data.frame(block_new)
      colnames(block_new) <- sprintf("%s_tp%s", uic_res[i, cause_var_colname], uic_res[i,"tp"])
      block_mvd <- cbind(block_mvd, block_new)
    }
  } else if (include_var == "strongest_only") {
    # ---------------------------------------------------- #
    # Select tp with the strongest influence from each causal variable
    # (te is used as a criterion)
    # ---------------------------------------------------- #
    for (cause_i in unique(uic_res[,cause_var_colname])) {
      uic_res_tmp <- uic_res[uic_res[,cause_var_colname] == cause_i,]
      if (!exists("uic_res_new")) {
        # uic_res_new <- uic_res_tmp[which.max(uic_res_tmp[,"te"]),]
        uic_res_new <- uic_res_tmp[which.max(uic_res_tmp$te),] # edited by Ohigashi
      } else {
        # uic_res_new <- rbind(uic_res_new, uic_res_tmp[which.max(uic_res_tmp[,"te"]),])
        uic_res_new <- rbind(uic_res_new, uic_res_tmp[which.max(uic_res_tmp$te),]) # edited by Ohigashi
      }
    }
    # Replace "uic_res"
    uic_res <- uic_res_new
    #uic_res <- uic_res %>% dplyr::group_by(.data$cause_var) %>% dplyr::filter(.data$te == max(.data$te))
    for (i in 1:nrow(uic_res)) {
      # block_new <- dplyr::lag(block[,uic_res[i,cause_var_colname]], n = abs(uic_res[i,"tp"]))
      block_new <- dplyr::lag(block[,unlist(uic_res[i,cause_var_colname])], n = abs(unlist(uic_res[i,"tp"]))) # edited by Ohigashi
      block_new <- data.frame(block_new)
      colnames(block_new) <- sprintf("%s_tp%s", uic_res[i,cause_var_colname], uic_res[i,"tp"])
      block_mvd <- cbind(block_mvd, block_new)
    }
  } else if (include_var == "tp0_only") {
    # ---------------------------------------------------- #
    # Select causal variables with tp = 0
    # ---------------------------------------------------- #
    block_new <- block[,unique(uic_res[,cause_var_colname])]
    block_new <- data.frame(block_new)
    colnames(block_new) <- sprintf("%s_tp0", unique(uic_res[,cause_var_colname]))
    block_mvd <- cbind(block_mvd, block_new)
  }
  # Return results
  return(block_mvd)
}



### function for uic_across but added "group" option
uic_across_TOed <- function(block,
                            effect_var,
                            cond_var = NULL,
                            E_range = 0:10,
                            tp_range = -4:0,
                            tau = 1,
                            num_surr = 1000,
                            alpha = 0.05,
                            fdr = FALSE,
                            group = NULL,
                            #random_seed = 1234,
                            silent = FALSE) {
  # Set random seed
  #set.seed(random_seed)
  
  # Retrieve colnames
  x_names <- colnames(block)
  
  # Check input arguments
  if (!all(unique(x_names) == x_names)) stop("\"block\" should have unique column names.")
  if (!is.data.frame(block)) stop("\"block\" should be data.frame.")
  if (!is.numeric(E_range) | !is.numeric(tp_range)) stop("\"E_range\" and \"tp_range\" should be numeric.")
  if (is.numeric(effect_var)) effect_var <- x_names[effect_var]
  
  # ---------------------------------------------------- #
  # Preparation for UIC
  # ---------------------------------------------------- #
  # Rearrange the column order
  block <- dplyr::select(block, tidyselect::all_of(effect_var), dplyr::everything())
  
  # ---------------------------------------------------- #
  # Identify causal relationship using rUIC
  # ---------------------------------------------------- #
  ## Perform UIC for all pairs
  # ---------------------------------------------------- #
  # uic.optimal()
  # ---------------------------------------------------- #
  # initialize uic result
  uic_res <- NULL
  
  # loop
  # y_looped <- x_names[x_names != effect_var & x_names != group & x_names != cond_var]
  y_looped <- setdiff(x_names, c(effect_var, group, cond_var))
  for (y_i in y_looped) {
    time_start <- proc.time()
    
    # Testing the effect of "y_i" on "effect_var" using uic.optimal()
    uic_xy <- rUIC::uic.optimal(block, lib_var = effect_var, tar_var = y_i, cond_var = cond_var,
                                E = E_range, tau = tau, tp = tp_range, num_surr = num_surr, alpha = alpha,
                                group = group) %>%
      dplyr::mutate(effect_var = effect_var, cause_var = y_i,
                    cond_var = if (is.null(cond_var)) NA_character_ else cond_var)
    # Combine results
    # if (y_i != x_names[x_names != effect_var][1]) { uic_res <- rbind(uic_res, uic_xy) } else { uic_res <- uic_xy }
    # if (y_i != y_looped[1]) { uic_res <- rbind(uic_res, uic_xy) } else { uic_res <- uic_xy }
    if (is.null(uic_res)) {
      uic_res <- uic_xy
    } else {
      uic_res <- rbind(uic_res, uic_xy)
    }
    # Output message
    time_used <- (proc.time() - time_start)[3]
    if (!silent) { message(sprintf("Effects from %s to %s tested by UIC: %.2f sec elapsed", y_i, effect_var, time_used)) }
  }
  
  if (fdr) {
    # Add FDR ("BH" method)
    uic_res$fdr <- stats::p.adjust(uic_res$pval, method = "BH")
  }
  
  # Message
  if (!silent) {
    message("By using `uic_across()` the effect variable is recorded in the column named `effect_var`, while potential causal variables are recorded in the column named `cause_var`.")
    message("Please note that these column names are used in the subsequent analysis if you use `make_block_mvd()`.")
  }
  
  # Return results
  return(as.data.frame(uic_res))
}

