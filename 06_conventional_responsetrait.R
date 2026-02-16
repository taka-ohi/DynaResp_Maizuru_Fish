####
#### R script for Ohigashi et al (2025)
#### calculate conventional environment-performance trait (linear model or GAM)
#### 2025.05.21 written by Ohigashi
#### R 4.5.0
####


### load packages
library(phyloseq); packageVersion("phyloseq")
library(dplyr); packageVersion("dplyr")
library(tibble); packageVersion("tibble")
library(mgcv); packageVersion("mgcv")
library(gratia); packageVersion("gratia")


### load data
# fish community data all time
ps_all <- readRDS("01_dataformatting_out/fishdata_ps.rds")


### calculate coefficients for lm and GAM
# extract OTU table from phyloseq
otu_mat <- as(otu_table(ps_all), "matrix")
# extract bottom temperature data
temp_vec <- sample_data(ps_all)$Water_temp_bottom

# calculate coefficients and store the result
coef_res <- lapply(colnames(otu_mat), function(fish) {
  df <- data.frame(
    Abundance = scale(otu_mat[, fish]),
    Temp = scale(temp_vec)
  )
  
  # lm and GAM
  lm_fit <- try(lm(Abundance ~ Temp, data = df), silent = TRUE)
  gam_fit <- try(gam(Abundance ~ s(Temp), data = df, family = gaussian()), silent = TRUE)
  
  # get 1st derivative of GAM model and average it across points (i.e., points with water temp)
  gam_slope <- try({
    d1 <- derivatives(gam_fit, term = "s(Temp)")
    mean(d1$.derivative, na.rm = TRUE)
  }, silent = TRUE)
  
  tibble(
    Fish_ID = fish,
    lm_coef = if (inherits(lm_fit, "lm")) coef(lm_fit)["Temp"] else NA,
    gam_slope_ave = if (!inherits(gam_slope, "try-error")) gam_slope else NA
  )
})

# convert to data frame
coef_res_df <- bind_rows(coef_res)


### save result
dir.create("06_conventional_responsetrait_out")
write.csv(coef_res_df, "06_conventional_responsetrait_out/lm_gam_coef.csv", quote = F, row.names = F)
saveRDS(coef_res_df, "06_conventional_responsetrait_out/lm_gam_coef.rds")


### save session info
writeLines(capture.output(sessionInfo()),
           # please change 0X or XX below to the script number you used.
           sprintf("00_SessionInfo/06_SessionInfo_%s.txt", substr(Sys.time(), 1, 10)))


