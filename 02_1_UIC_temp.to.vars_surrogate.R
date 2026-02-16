####
#### R script for Ohigashi et al (2025)
#### Analysis of causal effect from temperature to fish using UIC
#### 2025.05.19 written by Ohigashi
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


### load data
# fish phyloseq object
ps_all <- readRDS("01_dataformatting_out/fishdata_ps.rds")

# fish count data
fishcount <- otu_table(ps_all) %>% as.data.frame()

# sample data
sampledata <- data.frame(sample_data(ps_all))
sampledata <- sampledata %>% rownames_to_column("Sample")


### format data
# set rownames as Sample ID column
fishcount <- fishcount %>% rownames_to_column("Sample")
sample_df <- sampledata %>%
  # use only bottom temperature and remove unused variables
  select(Sample, Water_temp_bottom, Total_abun, Num_species, Shannon.d, Simpson.d,
         lat.high.total, lat.low.total, lat.high.rich, lat.low.rich)

# combine fish and sample data
sample_fishcount_df <- merge(sample_df, fishcount, by = "Sample", sort = F)


## create surrogate data for bottom water temperature and their differences
# set how to generate random numbers
set.seed(123)

# generate 1000 surrogate data for bottom water temperature
surr_bottomT <- make_surrogate_seasonal(sample_df$Water_temp_bottom,
                                        num_surr = 1000, T_period = 24)
colnames(surr_bottomT) <- sprintf("surr_bottomT%04d", seq(ncol(surr_bottomT)))
rownames(surr_bottomT) <- sprintf("Sample%03d", seq(nrow(surr_bottomT)))
surr_bottomT <- as.data.frame(surr_bottomT)


### preparation for looping
# recognize fixed vars (which means these vars do not become members of looping)
fixed_vars <- c("Sample", "Water_temp_bottom")

# set looping vars
looped_varnames <- setdiff(names(sample_fishcount_df), fixed_vars)

# scale data
surr_bottomT_sc <- scale(surr_bottomT)
looped_var_df <- sample_fishcount_df |> select(all_of(looped_varnames)) # extract looped variables
looped_var_sc <- scale(looped_var_df)


###### run UIC for surrogate data of bottom temperature ######
### to fish sp count and macroecology index
# check the number of dimensions
length(looped_varnames); dim(looped_var_sc); dim(surr_bottomT_sc)
all(looped_varnames == colnames(looped_var_sc))

# parallel calculation
uic_surr_bottomT_to_vars <- pblapply(1:length(looped_varnames), function(j) {
  # initialize 
  uic_res_all <- data.frame()
  
  # loop for surrogate
  for (i in 1:ncol(surr_bottomT_sc)) {
    # set method to generate random values
    set.seed(123)
    
    # compute UIC
    uic_res_tmp <- uic.optimal(data.frame(looped_var_sc[, j], surr_bottomT_sc[, i]),
                               lib_var = 1, tar_var = 2,
                               E = 0:24, tau = 1, tp = -12:2)
    
    # make a data frame
    # check variable names !!!
    uic_res.df <- data.frame(effect_var = looped_varnames[j], cause_var = colnames(surr_bottomT_sc)[i], uic_res_tmp)
    
    # combine result for one looped var
    uic_res_all <- rbind(uic_res_all, uic_res.df)
  }
  
  return(uic_res_all)  # return result
}, cl = 32)

# transform the list to a data frame
uic_surr_bt.to.vars <- do.call(rbind, uic_surr_bottomT_to_vars)

# save the result
dir.create("02_UIC_out")
saveRDS(uic_surr_bt.to.vars, "02_UIC_out/uic_surr_bottomT_to_vars.rds")


###### run UIC for OBSERVED data of bottom temperature ######
# prepare output object
uic_bt_obs <- data.frame()
# scale bottom temperature
bottomT_sc <- scale(sampledata$Water_temp_bottom)

# Calculate UIC
for (k in 1:length(looped_varnames)){
  set.seed(123)
  uic_res_tmp <- uic.optimal(data.frame(looped_var_sc[, k], bottomT_sc),
                             lib_var = 1, tar_var = 2,
                             E=0:24, tau=1, tp = -12:2
  ) %>%
    cbind(data.frame(effect_var = looped_varnames[k],
                     cause_var = "Water_temp_bottom"), .)
  uic_bt_obs <- rbind(uic_bt_obs, uic_res_tmp)
}

# save the result
saveRDS(uic_bt_obs, "02_UIC_out/uic_obs_bottomT_to_vars.rds")


### save environment
save(list = ls(all.names = TRUE), file = "02_UIC_out/UIC_temp.to.vars_surrogate.RData")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/02_1_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))



