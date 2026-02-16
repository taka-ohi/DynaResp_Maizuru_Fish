####
#### R script for Ohigashi et al (2025)
#### Analysis of causal effect from temperature to fish using UIC
#### 2025.05.20 written by Ohigashi
#### R 4.5.0
#### eliminate (mitigate) seasonality by using water temperature as a condition variable in UIC


### load packages
library(rUIC); packageVersion("rUIC")
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

# add columns for differences from previous time point
# fish count data
fish_count_w_diff <- fishcount |>
  mutate(across(2:ncol(fishcount),
                .fns = list(diff = ~ . -lag(.)),
                .names = "{col}.diff")
  )

# sample data
sample_df_w_diff <- sample_df |>
  mutate(across(2:ncol(sample_df),
                .fns = list(diff = ~ . -lag(.)),
                .names = "{col}.diff")
  )

# combine fish and sample data
sample_fishcount_diff <- merge(sample_df_w_diff, fish_count_w_diff, by = "Sample", sort = F)

# remove normal count data
sample_fishcount_diff <- sample_fishcount_diff %>% 
  select(Sample, Water_temp_bottom, any_of(ends_with(".diff")))

# delete first row (because it doesn't have difference values)
sample_fishcount_diff <- sample_fishcount_diff[-1,]


### preparation for looping
# scale variables
sa_fc_diff_sc <- data.frame(sample_fishcount_diff[, 1], scale(sample_fishcount_diff[, 2:ncol(sample_fishcount_diff)]))
colnames(sa_fc_diff_sc)[1] <- "Sample"

# recognize fixed vars (which means these vars do not become members of looping)
fixed_vars <- c("Sample", "Water_temp_bottom", "Water_temp_bottom.diff")

# set looping vars
looped_vars <- setdiff(names(sa_fc_diff_sc), fixed_vars)


####### run UIC with conditioned#######
## water bottom temperature difference -> each variable's difference* (conditional data: bottom temp)
res_bt.to.vars <- pblapply(1:length(looped_vars), function(var) {
  # set different seed among variables
  set.seed(123)
  
  # UIC 
  uic_opt_yx <- uic.optimal(sa_fc_diff_sc, lib_var = looped_vars[var], 
                            tar_var = "Water_temp_bottom.diff", # check!
                            cond_var = "Water_temp_bottom", # check!
                            E = 0:24, tau = 1, tp = -12:2)
  
  # make a data frame
  # check variable names !!!
  uic_res <- data.frame(effect_var = looped_vars[var], cause_var = "Water_temp_bottom.diff", cond_var = "Water_temp_bottom",
                        uic_opt_yx)
  return(uic_res)
}, cl = 32)

# transform the list to a data frame
uic_bt.diff_to_vars.diff <- do.call(rbind, res_bt.to.vars)
# save raw UIC result
saveRDS(uic_bt.diff_to_vars.diff, "02_UIC_out/uic_bottomT.diff_to_vars.diff_cond.temp_all.rds")


#### calculate adjusted p-values for fish
# divide the result into macrooecology and fishcount data
uic_bt.diff_res.macro <- uic_bt.diff_to_vars.diff |> filter(!grepl("Fish", effect_var)) # macroecology
uic_bt.diff_res.fc <- uic_bt.diff_to_vars.diff |> filter(grepl("Fish", effect_var)) # fish count

# put adjubted p values for fish species
# fish count
uic_bt.diff_res.fc <- uic_bt.diff_res.fc |>
  group_by(tp) |>
  mutate(p.adj = p.adjust(pval, method = "BH")) |>
  ungroup()

## save the data dividing into macroecology and fish count
# fish count
saveRDS(uic_bt.diff_res.fc, "02_UIC_out/uic_bottomT.diff_to_fc.diff_cond.temp.rds")
# macro ecology
saveRDS(uic_bt.diff_res.macro, "02_UIC_out/uic_bottomT.diff_to_macro.diff_cond.temp.rds")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/02_3_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))



