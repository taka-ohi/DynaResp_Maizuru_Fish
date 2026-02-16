####
#### R script for Ohigashi et al (2025)
#### Analysis of causal effect from temperature to fish using UIC
#### 2025.05.20 written by Ohigashi
#### R 4.5.0
#### summarize UIC results (calculate p-values with surrogate data)


### load packages and functions
source("Function/F1_HelperFunction_stats.R")
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(phyloseq); packageVersion("phyloseq")
library(pbapply); packageVersion("pbapply")


### load data
# uic surrogate result
uic_surr_bt_to_vars <- readRDS("02_UIC_out/uic_surr_bottomT_to_vars.rds")

# uic observed result
uic_bt_obs <- readRDS("02_UIC_out/uic_obs_bottomT_to_vars.rds")

# fish ecology information
ps_all <- readRDS("01_dataformatting_out/fishdata_ps.rds")
tax_sheet <- tax_table(ps_all) %>% as.data.frame()

### loop to calculate p-values for the observation data comparing their TE to surrogate data's E
# set variable names
varnames <- unique(uic_bt_obs$effect_var)

# get a combination of parameters (effect varnames x tp)
params_list <- expand.grid(i = varnames, j = -12:2)

# get surrogate names
surr_names_bt <- unique(uic_surr_bt_to_vars$cause_var)

# parallel calculation of p-values
results_bt <- pblapply(1:nrow(params_list), function(index) {
  params <- params_list[index, ]
  calculate_surr_p_ci(params = params,
                      uic_obs = uic_bt_obs,
                      uic_surr = uic_surr_bt_to_vars,
                      surr_names = surr_names_bt,
                      test_column = "te"
  )
}, cl = 32)

# put results in a data frame
results_bt.df <- do.call(rbind, lapply(results_bt, as.data.frame))

# prepare a data frame which includes the te and rmse to combine with the result
uic_bt_res <- uic_bt_obs |>
  select(effect_var, cause_var, obs_optE=E, tp, rmse, te) # preserve optimal E of observation data

# merge the result
uic_bt_res <- merge(uic_bt_res, results_bt.df, by = c("effect_var", "tp"), sort = F)
uic_bt_res <- uic_bt_res |>
  select(effect_var, cause_var, tp, names(uic_bt_res)[4:ncol(uic_bt_res)])

# divide the result into macrooecology and fishcount data
uic_bt_res.macro <- uic_bt_res |> filter(!grepl("Fish", effect_var)) # macroecology
uic_bt_res.fc <- uic_bt_res |> filter(grepl("Fish", effect_var)) # fish count

# put adjusted p values for fish species
# fish count
uic_bt_res.fc <- uic_bt_res.fc |>
  group_by(tp) |>
  mutate(surr_p.adj = p.adjust(surr_p, method = "BH")) |>
  ungroup()


### save data
## all result
# raw result
saveRDS(results_bt.df, "02_UIC_out/uic_bottomT_w_season.surrogate_rawres.rds")
# summarized w/ adjusted p
saveRDS(uic_bt_res, "02_UIC_out/uic_bottomT_w_season.surrogate.rds")

## save the data dividing into macroecology and fish count
# fish count
saveRDS(uic_bt_res.fc, "02_UIC_out/uic_bottomT_w_season.surrogate_fc.rds")
# macro ecology
saveRDS(uic_bt_res.macro, "02_UIC_out/uic_bottomT_w_season.surrogate_macro.rds")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/02_2_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))







