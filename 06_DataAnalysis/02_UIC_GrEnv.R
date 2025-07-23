####
#### Ogura Rice 2017
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2024.08.03
library(rUIC); packageVersion("rUIC") # 0.9.13, 2024.11.15

# Set random seeds (for reproduction)
set.seed(1234)

# Create output directory
(outdir <- macam::outdir_create())
dir.create("00_SessionInfo")


# ------------------------------------------------ #
# Load all data
# ------------------------------------------------ #
# All data
d_all <- readRDS("01_LoadAllDataOut/d_all.obj") %>% data.frame

# Compile sample names
d_all[,"Sample_Name2"][is.na(d_all$Sample_Name2)] <- sprintf("NoDNA%02s", 1:sum(is.na(d_all$Sample_Name2)))
rownames(d_all) <- d_all$Sample_Name2


# ------------------------------------------------ #
# Compile data for UIC between rice growth rate and environment
# ------------------------------------------------ #
uic_df <- d_all %>%
  select(date, uic_group, mean_height_diff, mean_n_stem_diff, mean_spad_diff,
         mean_height, mean_n_stem, mean_spad, starts_with("clim_")) %>% 
  filter(!is.na(uic_group)) %>% 
  select(-clim_cum_mean_temp, -clim_cum_rainfall, -clim_cum_sunlight_hours)

# Check data structure by visualizing them
uic_df %>% ggplot(aes(x = date, y = mean_height_diff, color = uic_group)) + geom_line()
uic_df %>% ggplot(aes(x = date, y = mean_n_stem_diff, color = uic_group)) + geom_line()
uic_df %>% ggplot(aes(x = date, y = mean_spad_diff, color = uic_group)) + geom_line()
uic_df %>% ggplot(aes(x = date, y = clim_mean_temp, color = uic_group)) + geom_line()
uic_df %>% ggplot(aes(x = date, y = clim_daily_rainfall, color = uic_group)) + geom_line()
uic_df %>% ggplot(aes(x = date, y = clim_sunlight_hours, color = uic_group)) + geom_line()
uic_df %>% ggplot(aes(x = date, y = clim_cum_7days_rainfall, color = uic_group)) + geom_line()

# Standardize all the variables
uic_df_std <- cbind(uic_df %>% select(date, uic_group), apply(uic_df[,3:ncol(uic_df)], 2, function(x) as.numeric(scale(x))))


# ------------------------------------------------ #
# UIC (Env => Rice)
# ------------------------------------------------ #
# Causality criteria
fdr_th <- 0.05
E_uic <- 0:5
tp_uic <- -2:0
clim_vars <- uic_df_std %>% select(starts_with("clim_")) %>% colnames

# UIC analysis
uic_height <- data.frame(NULL)
uic_spad <- data.frame(NULL)
uic_n_stem <- data.frame(NULL)
for (clim_var_i in clim_vars) {
  ## Mean height
  uic_height0 <- uic.optimal(uic_df_std, lib_var = "mean_height_diff", tar_var = clim_var_i, group = "uic_group", E = E_uic, tp = tp_uic, num_surr = 2000)
  uic_height1 <- uic_height0 %>% mutate(effect_var = "mean_height_diff", cause_var = clim_var_i, cond_var = NA)
  uic_height <- rbind(uic_height, uic_height1)
  ## SPAD
  uic_spad0 <- uic.optimal(uic_df_std, lib_var = "mean_spad_diff", tar_var = clim_var_i, group = "uic_group", E = E_uic, tp = tp_uic, num_surr = 2000)
  uic_spad1 <- uic_spad0 %>% mutate(effect_var = "mean_spad_diff", cause_var = clim_var_i, cond_var = NA)
  uic_spad <- rbind(uic_spad, uic_spad1)
  ## N. of stems
  uic_n_stem0 <- uic.optimal(uic_df_std, lib_var = "mean_n_stem_diff", tar_var = clim_var_i, group = "uic_group", E = E_uic, tp = tp_uic, num_surr = 2000)
  uic_n_stem1 <- uic_n_stem0 %>% mutate(effect_var = "mean_n_stem_diff", cause_var = clim_var_i, cond_var = NA)
  uic_n_stem <- rbind(uic_n_stem, uic_n_stem1)
}

rm(uic_height0); rm(uic_height1)
rm(uic_spad0); rm(uic_spad1)
rm(uic_n_stem0); rm(uic_n_stem1)
# Temperature, rainfall, and sunlight influence all the performance indicator of rice!


# ------------------------------------------------ #
# FDR: Total tests = 24 = 3 tp x 8 vars
# ------------------------------------------------ #
# Calculate FDR
# uic_height$fdr <- nrow(uic_height) * uic_height$pval / rank(uic_height$pval)
# uic_spad$fdr <- nrow(uic_spad) * uic_spad$pval / rank(uic_spad$pval)
# uic_n_stem$fdr <- nrow(uic_n_stem) * uic_n_stem$pval / rank(uic_n_stem$pval)
uic_height$fdr <- p.adjust(uic_height$pval, method = "BH")
uic_spad$fdr <- p.adjust(uic_spad$pval, method = "BH")
uic_n_stem$fdr <- p.adjust(uic_n_stem$pval, method = "BH")

# Visualize FDR
# plot(uic_height$pval ~ rank(uic_height$pval)); abline(h = 0.05)
# points(uic_height$fdr ~ rank(uic_height$pval), col = "red")
# qval1 <- uic_height$fdr
# qval2 <- qvalue::qvalue_truncp(uic_height$pval)$qvalue
# qval3 <- qvalue::qvalue_truncp(uic_n_stem$pval)$qvalue
# plot(qval1, qval2); abline(0,1)


# ------------------------------------------------ #
# Save results
# ------------------------------------------------ #
# R objects
saveRDS(uic_height, sprintf("%s/uic_height_res.obj", outdir))
saveRDS(uic_spad, sprintf("%s/uic_spad_res.obj", outdir))
saveRDS(uic_n_stem, sprintf("%s/uic_n_stem_res.obj", outdir))

# CSV files
write.csv(uic_height, sprintf("%s/CSV_uic_height_res.csv", outdir), row.names = F)
write.csv(uic_spad, sprintf("%s/CSV_uic_spad_res.csv", outdir), row.names = F)
write.csv(uic_n_stem, sprintf("%s/CSV_uic_n_stem_res.csv", outdir), row.names = F)

# Save workspace and session information
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))
macam::save_session_info()



# ------------------------------------------------ #
# ------------------------------------------------ #
# ------------------------------------------------ #
# Appendix (Rice => Env); All tests are insignificant
# ------------------------------------------------ #
# UIC (Rice => Env)
# ------------------------------------------------ #
# UIC analysis
uic_height_clim <- data.frame(NULL)
uic_spad_clim <- data.frame(NULL)
uic_n_stem_clim <- data.frame(NULL)
for (clim_var_i in clim_vars) {
  ## Mean height
  uic_height_clim0 <- uic.optimal(uic_df_std, tar_var = "mean_height_diff", lib_var = clim_var_i, group = "uic_group", E = E_uic, tp = tp_uic, num_surr = 2000)
  uic_height_clim1 <- uic_height_clim0 %>% mutate(cause_var = "mean_height_diff", effect_var = clim_var_i, cond_var = NA)
  uic_height_clim <- rbind(uic_height_clim, uic_height_clim1)
  ## SPAD
  uic_spad_clim0 <- uic.optimal(uic_df_std, tar_var = "mean_spad_diff", lib_var = clim_var_i, group = "uic_group", E = E_uic, tp = tp_uic, num_surr = 2000)
  uic_spad_clim1 <- uic_spad_clim0 %>% mutate(cause_var = "mean_spad_diff", effect_var = clim_var_i, cond_var = NA)
  uic_spad_clim <- rbind(uic_spad_clim, uic_spad_clim1)
  ## N. of stems
  uic_n_stem_clim0 <- uic.optimal(uic_df_std, tar_var = "mean_n_stem_diff", lib_var = clim_var_i, group = "uic_group", E = E_uic, tp = tp_uic, num_surr = 2000)
  uic_n_stem_clim1 <- uic_n_stem_clim0 %>% mutate(cause_var = "mean_n_stem_diff", effect_var = clim_var_i, cond_var = NA)
  uic_n_stem_clim <- rbind(uic_n_stem_clim, uic_n_stem_clim1)
}

rm(uic_height_clim0); rm(uic_height_clim1)
rm(uic_spad_clim0); rm(uic_spad_clim1)
rm(uic_n_stem_clim0); rm(uic_n_stem_clim1)


# ------------------------------------------------ #
# FDR: Total tests = 24 = 3 tp x 8 vars
# ------------------------------------------------ #
# Calculate FDR
uic_height_clim$fdr <- p.adjust(uic_height_clim$pval, method = "BH")
uic_spad_clim$fdr <- p.adjust(uic_spad_clim$pval, method = "BH")
uic_n_stem_clim$fdr <- p.adjust(uic_n_stem_clim$pval, method = "BH")

# CSV files
write.csv(uic_height_clim, sprintf("%s/CSV_uic_height_clim_res.csv", outdir), row.names = F)
write.csv(uic_spad_clim, sprintf("%s/CSV_uic_spad_clim_res.csv", outdir), row.names = F)
write.csv(uic_n_stem_clim, sprintf("%s/CSV_uic_n_stem_clim_res.csv", outdir), row.names = F)
