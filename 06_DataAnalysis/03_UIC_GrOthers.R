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

# Load climate data
uic_clim_hgt <- readRDS("02_UIC_GrEnvOut/uic_height_res.obj")
uic_clim_spd <- readRDS("02_UIC_GrEnvOut/uic_spad_res.obj")
uic_clim_stm <- readRDS("02_UIC_GrEnvOut/uic_n_stem_res.obj")
clim_vars <- c("clim_cum_7days_temp", # Highly causal for rice parameters
               "clim_cum_7days_sunlihg_hours") # Highly causal for rice parameters


# ------------------------------------------------ #
# Compile data for UIC between rice growth rate and environment
# ------------------------------------------------ #
uic_df <- d_all %>%
  select(date, uic_group,
         mean_height_diff, mean_n_stem_diff, mean_spad_diff,
         mean_height, mean_n_stem, mean_spad,
         starts_with("clim_"),
         starts_with("psmb_"),
         starts_with("grss_"),
         starts_with("insc_")) %>% 
  filter(!is.na(uic_group)) %>% 
  select(-clim_cum_mean_temp, -clim_cum_rainfall, -clim_cum_sunlight_hours)

# Interpolate missing insect data using linear-interpolation
uic_df <- uic_df %>%
  mutate(insc_other_planthopper_ip = approx(1:nrow(.), insc_other_planthopper, n = nrow(.))$y,
         insc_cicadellidae_spp_ip = approx(1:nrow(.), insc_cicadellidae_spp, n = nrow(.))$y,
         insc_thysanoptera_spp_ip = approx(1:nrow(.), insc_thysanoptera_spp, n = nrow(.))$y)

# Check data structure by visualizing them
uic_df %>% ggplot(aes(x = date, y = mean_height_diff, color = uic_group)) + geom_line()
uic_df %>% ggplot(aes(x = date, y = mean_n_stem_diff, color = uic_group)) + geom_line()
uic_df %>% ggplot(aes(x = date, y = mean_spad_diff, color = uic_group)) + geom_line()
uic_df %>% ggplot(aes(x = date, y = psmb_b_caryophyllene, color = uic_group)) + geom_line()
uic_df %>% ggplot(aes(x = date, y = grss_poa_1, color = uic_group)) + geom_line()
uic_df %>% ggplot(aes(x = date, y = insc_other_planthopper, color = uic_group)) + geom_line()
uic_df %>% ggplot(aes(x = date, y = insc_other_planthopper_ip, color = uic_group)) + geom_line()

# Standardize all the variables
uic_df_std <- cbind(uic_df %>% select(date, uic_group), apply(uic_df[,3:ncol(uic_df)], 2, function(x) as.numeric(scale(x))))

# Prepare time-delay for each climate variable
## Height
tp_hgt1 <- uic_clim_hgt %>% filter(cause_var == clim_vars[1]) %>% 
  filter(ete == max(ete)) %>% select(tp) %>% as.numeric
tp_hgt2 <- uic_clim_hgt %>% filter(cause_var == clim_vars[2]) %>% 
  filter(ete == max(ete)) %>% select(tp) %>% as.numeric
## SPAD
tp_spd1 <- uic_clim_spd %>% filter(cause_var == clim_vars[1]) %>% 
  filter(ete == max(ete)) %>% select(tp) %>% as.numeric
tp_spd2 <- uic_clim_spd %>% filter(cause_var == clim_vars[2]) %>% 
  filter(ete == max(ete)) %>% select(tp) %>% as.numeric
## N of stems
tp_stm1 <- uic_clim_stm %>% filter(cause_var == clim_vars[1]) %>% 
  filter(ete == max(ete)) %>% select(tp) %>% as.numeric
tp_stm2 <- uic_clim_stm %>% filter(cause_var == clim_vars[2]) %>% 
  filter(ete == max(ete)) %>% select(tp) %>% as.numeric


# ------------------------------------------------ #
# UIC
# ------------------------------------------------ #
# Causality criteria
fdr_th <- 0.05
E_uic <- 0:5
tp_uic <- -2:0
psmb_vars <- uic_df_std %>% select(starts_with("psmb_")) %>% colnames
grss_vars <- uic_df_std %>% select(starts_with("grss_")) %>% colnames
insc_vars <- uic_df_std %>% select(starts_with("insc_")) %>% select(ends_with("_ip")) %>% colnames
all_vars <- c(psmb_vars, grss_vars, insc_vars)


# ------------------------------------------------ #
# UIC: Other variables => Rice
# ------------------------------------------------ #
# Prepare df for each rice variable
uic_df_std_hgt <- uic_df_std_spd <- uic_df_std_stm <- uic_df_std
## Height
uic_df_std_hgt$clim_cum_7days_temp <- lag(uic_df_std$clim_cum_7days_temp, n = abs(tp_hgt1))
uic_df_std_hgt$clim_cum_7days_sunlihg_hours <- lag(uic_df_std$clim_cum_7days_sunlihg_hours, n = abs(tp_hgt2))
## SPAD
uic_df_std_spd$clim_cum_7days_temp <- lag(uic_df_std_spd$clim_cum_7days_temp, n = abs(tp_spd1))
uic_df_std_spd$clim_cum_7days_sunlihg_hours <- lag(uic_df_std_spd$clim_cum_7days_sunlihg_hours, n = abs(tp_spd2))
## N of stems
uic_df_std_stm$clim_cum_7days_temp <- lag(uic_df_std_stm$clim_cum_7days_temp, n = abs(tp_stm1))
uic_df_std_stm$clim_cum_7days_sunlihg_hours <- lag(uic_df_std_stm$clim_cum_7days_sunlihg_hours, n = abs(tp_stm2))

uic_others_height <- data.frame(NULL)
uic_others_spad <- data.frame(NULL)
uic_others_n_stem <- data.frame(NULL)
for (var_i in all_vars) {
  ## Mean height diff
  uic_others_height0 <- uic.optimal(uic_df_std_hgt, lib_var = "mean_height_diff", tar_var = var_i, cond_var = clim_vars, group = "uic_group", E = E_uic, tp = tp_uic, num_surr = 2000)
  uic_others_height1 <- uic_others_height0 %>% mutate(effect_var = "mean_height_diff", cause_var = var_i, cond_var = str_flatten(clim_vars, collapse = ", "))
  uic_others_height <- rbind(uic_others_height, uic_others_height1)
  ## SPAD diff
  uic_others_spad0 <- uic.optimal(uic_df_std_spd, lib_var = "mean_spad_diff", tar_var = var_i, cond_var = clim_vars, group = "uic_group", E = E_uic, tp = tp_uic, num_surr = 2000)
  uic_others_spad1 <- uic_others_spad0 %>% mutate(effect_var = "mean_spad_diff", cause_var = var_i, cond_var = str_flatten(clim_vars, collapse = ", "))
  uic_others_spad <- rbind(uic_others_spad, uic_others_spad1)
  ## N. of stems
  uic_others_n_stem0 <- uic.optimal(uic_df_std_stm, lib_var = "mean_n_stem_diff", tar_var = var_i, cond_var = clim_vars, group = "uic_group", E = E_uic, tp = tp_uic, num_surr = 2000)
  uic_others_n_stem1 <- uic_others_n_stem0 %>% mutate(effect_var = "mean_n_stem_diff", cause_var = var_i, cond_var = str_flatten(clim_vars, collapse = ", "))
  uic_others_n_stem <- rbind(uic_others_n_stem, uic_others_n_stem1)
}

## Delete temporal objects
rm(uic_others_height0); rm(uic_others_height1)
rm(uic_others_spad0); rm(uic_others_spad1)
rm(uic_others_n_stem0); rm(uic_others_n_stem1)

# Calculate FDR
uic_others_height$fdr <- p.adjust(uic_others_height$pval, method = "BH")
uic_others_spad$fdr <- p.adjust(uic_others_spad$pval, method = "BH")
uic_others_n_stem$fdr <- p.adjust(uic_others_n_stem$pval, method = "BH")


# ------------------------------------------------ #
# UIC: Rice => Other variables
# ------------------------------------------------ #
uic_others_height_rev <- data.frame(NULL)
uic_others_spad_rev <- data.frame(NULL)
uic_others_n_stem_rev <- data.frame(NULL)
for (var_i in all_vars) {
  ## Height diff
  uic_others_height_rev0 <- uic.optimal(uic_df_std_hgt, tar_var = "mean_height_diff", lib_var = var_i, cond_var = clim_vars, group = "uic_group", E = E_uic, tp = tp_uic, num_surr = 2000)
  uic_others_height_rev1 <- uic_others_height_rev0 %>% mutate(effect_var = var_i, cause_var = "mean_height_diff", cond_var = str_flatten(clim_vars, collapse = ", "))
  uic_others_height_rev <- rbind(uic_others_height_rev, uic_others_height_rev1)
  ## SPAD diff
  uic_others_spad_rev0 <- uic.optimal(uic_df_std_spd, tar_var = "mean_spad_diff", lib_var = var_i, cond_var = clim_vars, group = "uic_group",E = E_uic, tp = tp_uic, num_surr = 2000)
  uic_others_spad_rev1 <- uic_others_spad_rev0 %>% mutate(effect_var = var_i, cause_var = "mean_spad_diff", cond_var = str_flatten(clim_vars, collapse = ", "))
  uic_others_spad_rev <- rbind(uic_others_spad_rev, uic_others_spad_rev1)
  ## N. of stems diff
  uic_others_n_stem_rev0 <- uic.optimal(uic_df_std_stm, tar_var = "mean_n_stem_diff", lib_var = var_i, cond_var = clim_vars, group = "uic_group", E = E_uic, tp = tp_uic, num_surr = 2000)
  uic_others_n_stem_rev1 <- uic_others_n_stem_rev0 %>% mutate(effect_var = var_i, cause_var = "mean_n_stem_diff", cond_var = str_flatten(clim_vars, collapse = ", "))
  uic_others_n_stem_rev <- rbind(uic_others_n_stem_rev, uic_others_n_stem_rev1)
}

## Delete temporal objects
rm(uic_others_height_rev0); rm(uic_others_height_rev1)
rm(uic_others_spad_rev0); rm(uic_others_spad_rev1)
rm(uic_others_n_stem_rev0); rm(uic_others_n_stem_rev1)

# Calculate FDR
uic_others_height_rev$fdr <- p.adjust(uic_others_height_rev$pval, method = "BH")
uic_others_spad_rev$fdr <- p.adjust(uic_others_spad_rev$pval, method = "BH")
uic_others_n_stem_rev$fdr <- p.adjust(uic_others_n_stem_rev$pval, method = "BH")


# ------------------------------------------------ #
# Save results
# ------------------------------------------------ #
# R objects
saveRDS(uic_others_height, sprintf("%s/uic_others_height.obj", outdir))
saveRDS(uic_others_spad, sprintf("%s/uic_others_spad.obj", outdir))
saveRDS(uic_others_n_stem, sprintf("%s/uic_others_n_stem.obj", outdir))
saveRDS(uic_others_height_rev, sprintf("%s/uic_others_height_rev.obj", outdir))
saveRDS(uic_others_spad_rev, sprintf("%s/uic_others_spad_rev.obj", outdir))
saveRDS(uic_others_n_stem_rev, sprintf("%s/uic_others_n_stem_rev.obj", outdir))

# CSV files
write.csv(uic_others_height, sprintf("%s/CSV_uic_others_height.csv", outdir), row.names = F)
write.csv(uic_others_spad, sprintf("%s/CSV_uic_others_spad.csv", outdir), row.names = F)
write.csv(uic_others_n_stem, sprintf("%s/CSV_uic_others_n_stem.csv", outdir), row.names = F)
write.csv(uic_others_height_rev, sprintf("%s/CSV_uic_others_height_rev.csv", outdir), row.names = F)
write.csv(uic_others_spad_rev, sprintf("%s/CSV_uic_others_spad_rev.csv", outdir), row.names = F)
write.csv(uic_others_n_stem_rev, sprintf("%s/CSV_uic_others_n_stem_rev.csv", outdir), row.names = F)

# Save workspace and session information
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))
macam::save_session_info()

