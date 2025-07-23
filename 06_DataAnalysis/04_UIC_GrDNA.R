####
#### Ogura Rice 2017
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2024.08.03
library(phyloseq); packageVersion("phyloseq") # 1.46.0, 2024.11.12
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
ps_ecol <- readRDS("01_LoadAllDataOut/ps_ecol.obj")
d_ecol <- otu_table(ps_ecol) %>% data.frame

# Compile sample names
d_all[,"Sample_Name2"][is.na(d_all$Sample_Name2)] <- sprintf("NoDNA%02s", 1:sum(is.na(d_all$Sample_Name2)))
rownames(d_all) <- d_all$Sample_Name2


# ------------------------------------------------ #
# Compile rice data and eDNA data
# ------------------------------------------------ #
# Load climate data
uic_clim_hgt <- readRDS("02_UIC_GrEnvOut/uic_height_res.obj")
uic_clim_spd <- readRDS("02_UIC_GrEnvOut/uic_spad_res.obj")
uic_clim_stm <- readRDS("02_UIC_GrEnvOut/uic_n_stem_res.obj")

# Select data for UIC
uic_df <- d_all %>%
  select(date, uic_group, sample_code,
         mean_height_diff, mean_n_stem_diff, mean_spad_diff,
         mean_height, mean_n_stem, mean_spad,
         starts_with("clim_")) %>% 
  #filter(!is.na(uic_group)) %>% 
  select(-clim_cum_mean_temp, -clim_cum_rainfall, -clim_cum_sunlight_hours)
# Standardize all the variables
uic_df_std0 <- cbind(uic_df %>% select(date, uic_group, sample_code), apply(uic_df[,4:ncol(uic_df)], 2, function(x) as.numeric(scale(x))))

# Retrieve sample_code
edna_id <- ps_ecol %>% sample_data %>% pull(sample_code)


# ------------------------------------------------ #
# UIC preparation
# ------------------------------------------------ #
# Causality criteria
fdr_th <- 0.05
E_uic <- 0:5
tp_uic <- -2:0
clim_vars <- c("clim_cum_7days_temp", # Highly causal for rice parameters
               "clim_cum_7days_sunlihg_hours") # Highly causal for rice parameters
ecol_vars <- colnames(d_ecol)

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
# UIC: eDNA => Rice
# ------------------------------------------------ #
uic_edna_height <- data.frame(NULL)
uic_edna_spad <- data.frame(NULL)
uic_edna_n_stem <- data.frame(NULL)
for (var_i in ecol_vars) {
  # Collect message
  time_start <- proc.time()[3]
  i <- match(var_i, ecol_vars)
  
  # Add a temporal eDNA conc column
  uic_df_std0$edna_conc <- NA
  uic_df_std0$edna_conc[match(edna_id, uic_df_std0$sample_code)] <- d_ecol[,var_i]
  # Standardization
  uic_df_std0$edna_conc_std <- as.numeric(scale(uic_df_std0$edna_conc))
  # Remove NA samples
  uic_df_std <- uic_df_std0 %>% filter(!is.na(uic_group))
  # Prepare df for each rice variable
  uic_df_std_hgt <- uic_df_std_spd <- uic_df_std_stm <- uic_df_std
  ## Height
  uic_df_std_hgt$clim_cum_7days_temp <- lag(uic_df_std_hgt$clim_cum_7days_temp, n = abs(tp_hgt1))
  uic_df_std_hgt$clim_cum_7days_sunlihg_hours <- lag(uic_df_std_hgt$clim_cum_7days_sunlihg_hours, n = abs(tp_hgt2))
  ## SPAD
  uic_df_std_spd$clim_cum_7days_temp <- lag(uic_df_std_spd$clim_cum_7days_temp, n = abs(tp_spd1))
  uic_df_std_spd$clim_cum_7days_sunlihg_hours <- lag(uic_df_std_spd$clim_cum_7days_sunlihg_hours, n = abs(tp_spd2))
  ## N of stems
  uic_df_std_stm$clim_cum_7days_temp <- lag(uic_df_std_stm$clim_cum_7days_temp, n = abs(tp_stm1))
  uic_df_std_stm$clim_cum_7days_sunlihg_hours <- lag(uic_df_std_stm$clim_cum_7days_sunlihg_hours, n = abs(tp_stm2))
  
  ## Mean height diff
  uic_edna_height0 <- uic.optimal(uic_df_std_hgt, lib_var = "mean_height_diff", tar_var = "edna_conc_std", cond_var = clim_vars, group = "uic_group", E = E_uic, tp = tp_uic, num_surr = 2000)
  uic_edna_height1 <- uic_edna_height0 %>% mutate(effect_var = "mean_height_diff", cause_var = var_i, cond_var = str_flatten(clim_vars, collapse = ", "))
  uic_edna_height <- rbind(uic_edna_height, uic_edna_height1)
  ## SPAD diff
  uic_edna_spad0 <- uic.optimal(uic_df_std_spd, lib_var = "mean_spad_diff", tar_var = "edna_conc_std", cond_var = clim_vars, group = "uic_group", E = E_uic, tp = tp_uic, num_surr = 2000)
  uic_edna_spad1 <- uic_edna_spad0 %>% mutate(effect_var = "mean_spad_diff", cause_var = var_i, cond_var = str_flatten(clim_vars, collapse = ", "))
  uic_edna_spad <- rbind(uic_edna_spad, uic_edna_spad1)
  ## N. of stems
  uic_edna_n_stem0 <- uic.optimal(uic_df_std_stm, lib_var = "mean_n_stem_diff", tar_var = "edna_conc_std", cond_var = clim_vars, group = "uic_group", E = E_uic, tp = tp_uic, num_surr = 2000)
  uic_edna_n_stem1 <- uic_edna_n_stem0 %>% mutate(effect_var = "mean_n_stem_diff", cause_var = var_i, cond_var = str_flatten(clim_vars, collapse = ", "))
  uic_edna_n_stem <- rbind(uic_edna_n_stem, uic_edna_n_stem1)
  
  # Report progress
  message(sprintf("Process %s in %s finished: %.1f sec", i, length(ecol_vars), proc.time()[3] - time_start))
}

## Delete temporal objects
rm(uic_edna_height0); rm(uic_edna_height1)
rm(uic_edna_spad0); rm(uic_edna_spad1)
rm(uic_edna_n_stem0); rm(uic_edna_n_stem1)

# Calculate FDR
uic_edna_height$fdr <- p.adjust(uic_edna_height$pval, method = "BH")
uic_edna_spad$fdr <- p.adjust(uic_edna_spad$pval, method = "BH")
uic_edna_n_stem$fdr <- p.adjust(uic_edna_n_stem$pval, method = "BH")


# ------------------------------------------------ #
# UIC: Rice => eDNA
# ------------------------------------------------ #
uic_edna_height_rev <- data.frame(NULL)
uic_edna_spad_rev <- data.frame(NULL)
uic_edna_n_stem_rev <- data.frame(NULL)
for (var_i in ecol_vars) {
  # Collect message
  time_start <- proc.time()[3]
  i <- match(var_i, ecol_vars)
  
  # Add a temporal eDNA conc column
  uic_df_std0$edna_conc <- NA
  uic_df_std0$edna_conc[match(edna_id, uic_df_std0$sample_code)] <- d_ecol[,var_i]
  # Standardization
  uic_df_std0$edna_conc_std <- as.numeric(scale(uic_df_std0$edna_conc))
  # Remove NA samples
  uic_df_std <- uic_df_std0 %>% filter(!is.na(uic_group))
  # Prepare df for each rice variable
  uic_df_std_hgt <- uic_df_std_spd <- uic_df_std_stm <- uic_df_std
  ## Height
  uic_df_std_hgt$clim_cum_7days_temp <- lag(uic_df_std_hgt$clim_cum_7days_temp, n = abs(tp_hgt1))
  uic_df_std_hgt$clim_cum_7days_sunlihg_hours <- lag(uic_df_std_hgt$clim_cum_7days_sunlihg_hours, n = abs(tp_hgt2))
  ## SPAD
  uic_df_std_spd$clim_cum_7days_temp <- lag(uic_df_std_spd$clim_cum_7days_temp, n = abs(tp_spd1))
  uic_df_std_spd$clim_cum_7days_sunlihg_hours <- lag(uic_df_std_spd$clim_cum_7days_sunlihg_hours, n = abs(tp_spd2))
  ## N of stems
  uic_df_std_stm$clim_cum_7days_temp <- lag(uic_df_std_stm$clim_cum_7days_temp, n = abs(tp_stm1))
  uic_df_std_stm$clim_cum_7days_sunlihg_hours <- lag(uic_df_std_stm$clim_cum_7days_sunlihg_hours, n = abs(tp_stm2))
  
  ## Height diff
  uic_edna_height_rev0 <- uic.optimal(uic_df_std_hgt, tar_var = "mean_height_diff", lib_var = "edna_conc_std", cond_var = clim_vars, group = "uic_group", E = E_uic, tp = tp_uic, num_surr = 2000)
  uic_edna_height_rev1 <- uic_edna_height_rev0 %>% mutate(effect_var = var_i, cause_var = "mean_height_diff", cond_var = str_flatten(clim_vars, collapse = ", "))
  uic_edna_height_rev <- rbind(uic_edna_height_rev, uic_edna_height_rev1)
  ## SPAD diff
  uic_edna_spad_rev0 <- uic.optimal(uic_df_std_spd, tar_var = "mean_spad_diff", lib_var = "edna_conc_std", cond_var = clim_vars, group = "uic_group",E = E_uic, tp = tp_uic, num_surr = 2000)
  uic_edna_spad_rev1 <- uic_edna_spad_rev0 %>% mutate(effect_var = var_i, cause_var = "mean_spad_diff", cond_var = str_flatten(clim_vars, collapse = ", "))
  uic_edna_spad_rev <- rbind(uic_edna_spad_rev, uic_edna_spad_rev1)
  ## N. of stems diff
  uic_edna_n_stem_rev0 <- uic.optimal(uic_df_std_stm, tar_var = "mean_n_stem_diff", lib_var = "edna_conc_std", cond_var = clim_vars, group = "uic_group", E = E_uic, tp = tp_uic, num_surr = 2000)
  uic_edna_n_stem_rev1 <- uic_edna_n_stem_rev0 %>% mutate(effect_var = var_i, cause_var = "mean_n_stem_diff", cond_var = str_flatten(clim_vars, collapse = ", "))
  uic_edna_n_stem_rev <- rbind(uic_edna_n_stem_rev, uic_edna_n_stem_rev1)
  
  # Report progress
  message(sprintf("Process %s in %s finished: %.1f sec", i, length(ecol_vars), proc.time()[3] - time_start))
}


## Delete temporal objects
rm(uic_edna_height_rev0); rm(uic_edna_height_rev1)
rm(uic_edna_spad_rev0); rm(uic_edna_spad_rev1)
rm(uic_edna_n_stem_rev0); rm(uic_edna_n_stem_rev1)

# Calculate FDR
uic_edna_height_rev$fdr <- p.adjust(uic_edna_height_rev$pval, method = "BH")
uic_edna_spad_rev$fdr <- p.adjust(uic_edna_spad_rev$pval, method = "BH")
uic_edna_n_stem_rev$fdr <- p.adjust(uic_edna_n_stem_rev$pval, method = "BH")


# ------------------------------------------------ #
# Save results
# ------------------------------------------------ #
# R objects
saveRDS(uic_edna_height, sprintf("%s/uic_edna_height.obj", outdir))
saveRDS(uic_edna_spad, sprintf("%s/uic_edna_spad.obj", outdir))
saveRDS(uic_edna_n_stem, sprintf("%s/uic_edna_n_stem.obj", outdir))
saveRDS(uic_edna_height_rev, sprintf("%s/uic_edna_height_rev.obj", outdir))
saveRDS(uic_edna_spad_rev, sprintf("%s/uic_edna_spad_rev.obj", outdir))
saveRDS(uic_edna_n_stem_rev, sprintf("%s/uic_edna_n_stem_rev.obj", outdir))

# CSV files
write.csv(uic_edna_height, sprintf("%s/CSV_uic_edna_height.csv", outdir), row.names = F)
write.csv(uic_edna_spad, sprintf("%s/CSV_uic_edna_spad.csv", outdir), row.names = F)
write.csv(uic_edna_n_stem, sprintf("%s/CSV_uic_edna_n_stem.csv", outdir), row.names = F)
write.csv(uic_edna_height_rev, sprintf("%s/CSV_uic_edna_height_rev.csv", outdir), row.names = F)
write.csv(uic_edna_spad_rev, sprintf("%s/CSV_uic_edna_spad_rev.csv", outdir), row.names = F)
write.csv(uic_edna_n_stem_rev, sprintf("%s/CSV_uic_edna_n_stem_rev.csv", outdir), row.names = F)

# Save workspace and session information
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))
macam::save_session_info()

