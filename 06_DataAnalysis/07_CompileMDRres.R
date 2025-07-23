####
#### Ogura Rice 2017
#### Compile and visualize MDR S-map results
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2024.08.03
library(phyloseq); packageVersion("phyloseq") # 1.46.0, 2024.11.12
library(ggforce); packageVersion("ggforce") # 0.4.2, 2024.12.09
library(patchwork); packageVersion("patchwork") # 1.3.0, 2024.12.04

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

# Load MDR S-map results
mdr_hgt <- readRDS("06_MDRsmapOut/mdr_res_mean_height_diff.obj")
mdr_spd <- readRDS("06_MDRsmapOut/mdr_res_mean_spad_diff.obj")
mdr_stm <- readRDS("06_MDRsmapOut/mdr_res_mean_n_stem_diff.obj")


# ------------------------------------------------ #
# Compile the MDR S-map results
# ------------------------------------------------ #
## Height
mdr_hgt2 <- d_all %>%
  select(Sample_Name2, sample_code, date, n_week, treatment, sampling_point) %>% 
  bind_cols(mdr_hgt$model_output %>% select(obs, pred))
pred_tmp <- mdr_hgt$block_mvd * mdr_hgt$smap_coefficients[,2:(ncol(mdr_hgt$smap_coefficients)-1)]
pred_tmp <- apply(pred_tmp, 1, function(x) sum(x)) + mdr_hgt$smap_coefficients[,ncol(mdr_hgt$smap_coefficients)]
mdr_hgt2$man_pred <- lag(pred_tmp, 1); rm(pred_tmp)
mdr_hgt2 <- cbind(mdr_hgt2, mdr_hgt$block_mvd, mdr_hgt$smap_coefficients %>% select(starts_with("c_")))
colnames(mdr_hgt2)[(10 + ncol(mdr_hgt$block_mvd)):(9 + 2*ncol(mdr_hgt$block_mvd))] <- colnames(mdr_hgt$block_mvd) %>% paste0("smapc_from_", .)

# SPAD
mdr_spd2 <- d_all %>%
  select(Sample_Name2, sample_code, date, n_week, treatment, sampling_point) %>% 
  bind_cols(mdr_spd$model_output %>% select(obs, pred))
pred_tmp <- mdr_spd$block_mvd * mdr_spd$smap_coefficients[,2:(ncol(mdr_spd$smap_coefficients)-1)]
pred_tmp <- apply(pred_tmp, 1, function(x) sum(x)) + mdr_spd$smap_coefficients[,ncol(mdr_spd$smap_coefficients)]
mdr_spd2$man_pred <- lag(pred_tmp, 1); rm(pred_tmp)
mdr_spd2 <- cbind(mdr_spd2, mdr_spd$block_mvd, mdr_spd$smap_coefficients %>% select(starts_with("c_")))
colnames(mdr_spd2)[(10 + ncol(mdr_spd$block_mvd)):(9 + 2*ncol(mdr_spd$block_mvd))] <- colnames(mdr_spd$block_mvd) %>% paste0("smapc_from_", .)

# N of stems
mdr_stm2 <- d_all %>%
  select(Sample_Name2, sample_code, date, n_week, treatment, sampling_point) %>% 
  bind_cols(mdr_stm$model_output %>% select(obs, pred))
pred_tmp <- mdr_stm$block_mvd * mdr_stm$smap_coefficients[,2:(ncol(mdr_stm$smap_coefficients)-1)]
pred_tmp <- apply(pred_tmp, 1, function(x) sum(x)) + mdr_spd$smap_coefficients[,ncol(mdr_spd$smap_coefficients)]
mdr_stm2$man_pred <- lag(pred_tmp, 1); rm(pred_tmp)
mdr_stm2 <- cbind(mdr_stm2, mdr_stm$block_mvd, mdr_stm$smap_coefficients %>% select(starts_with("c_")))
colnames(mdr_stm2)[(10 + ncol(mdr_stm$block_mvd)):(9 + 2*ncol(mdr_stm$block_mvd))] <- colnames(mdr_stm$block_mvd) %>% paste0("smapc_from_", .)

# Convert data to a longer format
## Height
mdr_hgt_long <- mdr_hgt2 %>%
  select(Sample_Name2, sample_code, date, n_week, treatment, sampling_point, starts_with("smapc_")) %>%
  pivot_longer(cols = -c(Sample_Name2, sample_code, date, n_week, treatment, sampling_point),
               names_to = "cause_var", values_to = "smapc")
## SPAD
mdr_spd_long <- mdr_spd2 %>%
  select(Sample_Name2, sample_code, date, n_week, treatment, sampling_point, starts_with("smapc_")) %>%
  pivot_longer(cols = -c(Sample_Name2, sample_code, date, n_week, treatment, sampling_point),
               names_to = "cause_var", values_to = "smapc")
## N of stems
mdr_stm_long <- mdr_stm2 %>%
  select(Sample_Name2, sample_code, date, n_week, treatment, sampling_point, starts_with("smapc_")) %>%
  pivot_longer(cols = -c(Sample_Name2, sample_code, date, n_week, treatment, sampling_point),
               names_to = "cause_var", values_to = "smapc")

# Extract variable names
## Height
cause_var_hgt <- mdr_hgt2 %>% select(starts_with("smapc_")) %>%
  colnames %>% str_split(pattern = "smapc_from_") %>% sapply(`[`, 2)
## SPAD
cause_var_spd <- mdr_spd2 %>% select(starts_with("smapc_")) %>%
  colnames %>% str_split(pattern = "smapc_from_") %>% sapply(`[`, 2)
## N. of stems
cause_var_stm <- mdr_stm2 %>% select(starts_with("smapc_")) %>%
  colnames %>% str_split(pattern = "smapc_from_") %>% sapply(`[`, 2)

# Height
mdr_hgt_long2 <- mdr_hgt2 %>%
  select(Sample_Name2, sample_code, date, n_week, treatment, sampling_point, all_of(cause_var_hgt)) %>% 
  pivot_longer(cols = -c(Sample_Name2, sample_code, date, n_week, treatment, sampling_point),
               names_to = "cause_var", values_to = "raw_value")
id1 <- paste0(mdr_hgt_long$Sample_Name2, "_", mdr_hgt_long$sample_code, 
              mdr_hgt_long$cause_var %>% str_split(pattern = "smapc_from_") %>% sapply(`[`, 2))
id2 <- paste0(mdr_hgt_long2$Sample_Name2, "_", mdr_hgt_long2$sample_code, mdr_hgt_long2$cause_var)
all(id1 == id2); rm(id1); rm(id2)
mdr_hgt_long2$smapc <- mdr_hgt_long$smapc

# SPAD
mdr_spd_long2 <- mdr_spd2 %>%
  select(Sample_Name2, sample_code, date, n_week, treatment, sampling_point, all_of(cause_var_spd)) %>% 
  pivot_longer(cols = -c(Sample_Name2, sample_code, date, n_week, treatment, sampling_point),
               names_to = "cause_var", values_to = "raw_value")
id1 <- paste0(mdr_spd_long$Sample_Name2, "_", mdr_spd_long$sample_code, 
              mdr_spd_long$cause_var %>% str_split(pattern = "smapc_from_") %>% sapply(`[`, 2))
id2 <- paste0(mdr_spd_long2$Sample_Name2, "_", mdr_spd_long2$sample_code, mdr_spd_long2$cause_var)
all(id1 == id2); rm(id1); rm(id2)
mdr_spd_long2$smapc <- mdr_spd_long$smapc

# N of stems
mdr_stm_long2 <- mdr_stm2 %>%
  select(Sample_Name2, sample_code, date, n_week, treatment, sampling_point, all_of(cause_var_stm)) %>% 
  pivot_longer(cols = -c(Sample_Name2, sample_code, date, n_week, treatment, sampling_point),
               names_to = "cause_var", values_to = "raw_value")
id1 <- paste0(mdr_stm_long$Sample_Name2, "_", mdr_stm_long$sample_code, 
              mdr_stm_long$cause_var %>% str_split(pattern = "smapc_from_") %>% sapply(`[`, 2))
id2 <- paste0(mdr_stm_long2$Sample_Name2, "_", mdr_stm_long2$sample_code, mdr_stm_long2$cause_var)
all(id1 == id2); rm(id1); rm(id2)
mdr_stm_long2$smapc <- mdr_stm_long$smapc


# ------------------------------------------------ #
# Pred v.s. obs
# ------------------------------------------------ #
# Height diff
h1 <- mdr_hgt2 %>% ggplot(aes(x = obs, y = pred)) +
  geom_point(size = 0) +
  stat_smooth(method = "lm", se = FALSE, color = "gray") +
  geom_point(aes(x = obs, y = pred, color = treatment, shape = treatment), size = 3) +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2")) +
  labs(x = "Observed", y = "Predicted", title = "Height_diff") +
  NULL
# Height diff
h2 <- mdr_spd2 %>% ggplot(aes(x = obs, y = pred)) +
  geom_point(size = 0) +
  stat_smooth(method = "lm", se = FALSE, color = "gray") +
  geom_point(aes(x = obs, y = pred, color = treatment, shape = treatment), size = 3) +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2")) +
  labs(x = "Observed", y = "Predicted", title = "SPAD_diff") +
  NULL
# Height diff
h3 <- mdr_stm2 %>% ggplot(aes(x = obs, y = pred)) +
  geom_point(size = 0) +
  stat_smooth(method = "lm", se = FALSE, color = "gray") +
  geom_point(aes(x = obs, y = pred, color = treatment, shape = treatment), size = 3) +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2")) +
  labs(x = "Observed", y = "Predicted", title = "N_Stem_diff") +
  NULL


# ------------------------------------------------ #
# Visualization
# ------------------------------------------------ #
## Height
## Boxplot + Jitterplot
g1 <- mdr_hgt_long2 %>%
  ggplot(aes(x = treatment, y = smapc, group = interaction(cause_var, treatment), shape = treatment)) +
  facet_wrap(. ~ cause_var, scales = "free") +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2")) +
  geom_boxplot(outliers = FALSE) +
  geom_sina(alpha = 0.5, aes(color = treatment)) +
  #geom_jitter(height = 0, width = 0.2, alpha = 0.5, aes(color = treatment)) +
  ggtitle("Height") +
  NULL
g2 <- mdr_spd_long2 %>%
  ggplot(aes(x = treatment, y = smapc, group = interaction(cause_var, treatment), shape = treatment)) +
  facet_wrap(. ~ cause_var, scales = "free") +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2")) +
  geom_boxplot(outliers = FALSE) +
  geom_sina(alpha = 0.5, aes(color = treatment)) +
  #geom_jitter(height = 0, width = 0.2, alpha = 0.5, aes(color = treatment)) +
  ggtitle("SPAD") +
  NULL
g3 <- mdr_stm_long2 %>%
  ggplot(aes(x = treatment, y = smapc, group = interaction(cause_var, treatment), shape = treatment)) +
  facet_wrap(. ~ cause_var, scales = "free") +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2")) +
  geom_boxplot(outliers = FALSE) +
  geom_sina(alpha = 0.5, aes(color = treatment)) +
  #geom_jitter(height = 0, width = 0.2, alpha = 0.5, aes(color = treatment)) +
  ggtitle("N of stems") +
  NULL

## Scattered plot
## Height
g4 <- mdr_hgt_long2 %>%
  ggplot(aes(x = raw_value, y = smapc, color = treatment, shape = treatment)) +
  facet_wrap(. ~ cause_var, scales = "free") +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2")) +
  stat_smooth(method = "lm", se = FALSE) +
  geom_point() + ggtitle("Height") +
  NULL
g1_2 <- mdr_hgt_long2 %>%
  ggplot(aes(x = treatment, y = smapc*raw_value, group = interaction(cause_var, treatment), shape = treatment)) +
  facet_wrap(. ~ cause_var, scales = "free") +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2")) +
  geom_boxplot(outliers = FALSE) +
  geom_sina(alpha = 0.5, aes(color = treatment)) +
  #geom_jitter(height = 0, width = 0.2, alpha = 0.5, aes(color = treatment)) +
  ggtitle("Height") +
  NULL
## SPAD
g5 <- mdr_spd_long2 %>%
  ggplot(aes(x = raw_value, y = smapc, color = treatment, shape = treatment)) +
  facet_wrap(. ~ cause_var, scales = "free") +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2")) +
  stat_smooth(method = "lm", se = FALSE) +
  geom_point() + ggtitle("SPAD") +
  NULL
g2_2 <- mdr_spd_long2 %>%
  ggplot(aes(x = treatment, y = smapc*raw_value, group = interaction(cause_var, treatment), shape = treatment)) +
  facet_wrap(. ~ cause_var, scales = "free") +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2")) +
  geom_boxplot(outliers = FALSE) +
  geom_sina(alpha = 0.5, aes(color = treatment)) +
  #geom_jitter(height = 0, width = 0.2, alpha = 0.5, aes(color = treatment)) +
  ggtitle("SPAD") +
  NULL
## N of stems
g6 <- mdr_stm_long2 %>%
  ggplot(aes(x = raw_value, y = smapc, color = treatment, shape = treatment)) +
  facet_wrap(. ~ cause_var, scales = "free") +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2")) +
  stat_smooth(method = "lm", se = FALSE) +
  geom_point() + ggtitle("N of stems") +
  NULL
g3_2 <- mdr_stm_long2 %>%
  ggplot(aes(x = treatment, y = smapc*raw_value, group = interaction(cause_var, treatment), shape = treatment)) +
  facet_wrap(. ~ cause_var, scales = "free") +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2")) +
  geom_boxplot(outliers = FALSE) +
  geom_sina(alpha = 0.5, aes(color = treatment)) +
  #geom_jitter(height = 0, width = 0.2, alpha = 0.5, aes(color = treatment)) +
  ggtitle("N of stems") +
  NULL


# ------------------------------------------------ #
# Save figures
# ------------------------------------------------ #
# Pred v.s. Obs
ggsave("07_CompileMDRresOut/PredObs.pdf",
       (h1 + theme(legend.position = "none")) +
        (h2 + theme(legend.position = "none")) +
         (h3), width = 16, height = 6)
# Boxplot + Jitterplot
ggsave("07_CompileMDRresOut/height_boxjitter.pdf", g1, width = 14, height = 12)
ggsave("07_CompileMDRresOut/SPAD_boxjitter.pdf", g2, width = 14, height = 12)
ggsave("07_CompileMDRresOut/nstems_boxjitter.pdf", g3, width = 14, height = 12)
# Boxplot + Jitterplot (smapc * raw_value)
ggsave("07_CompileMDRresOut/height_boxjitter2.pdf", g1_2, width = 14, height = 12)
ggsave("07_CompileMDRresOut/SPAD_boxjitter2.pdf", g2_2, width = 14, height = 12)
ggsave("07_CompileMDRresOut/nstems_boxjitter2.pdf", g3_2, width = 14, height = 12)
# Scattered plot
ggsave("07_CompileMDRresOut/height_scatter.pdf", g4, width = 14, height = 12)
ggsave("07_CompileMDRresOut/SPAD_scatter.pdf", g5, width = 14, height = 12)
ggsave("07_CompileMDRresOut/nstems_scatter.pdf", g6, width = 14, height = 12)
# Save R objects
saveRDS(list(h1, h2, h3), "../07_FormatFigs/data_robj/Smap_PredObs.obj")
h_all <- (h1 + theme(legend.position = "none")) + (h2 + theme(legend.position = "none")) + h3
saveRDS(h_all, "../07_FormatFigs/data_robj/Smap_PredObs-all.obj")


# ------------------------------------------------ #
# Save compiled data
# ------------------------------------------------ #
# Wide objects
saveRDS(mdr_hgt2, "07_CompileMDRresOut/mdr_hgt.obj")
saveRDS(mdr_spd2, "07_CompileMDRresOut/mdr_spd.obj")
saveRDS(mdr_stm2, "07_CompileMDRresOut/mdr_stm.obj")
# Long objects
saveRDS(mdr_hgt_long2, "07_CompileMDRresOut/mdr_hgt_long.obj")
saveRDS(mdr_spd_long2, "07_CompileMDRresOut/mdr_spd_long.obj")
saveRDS(mdr_stm_long2, "07_CompileMDRresOut/mdr_stm_long.obj")

# Save workspace and session information
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))
macam::save_session_info()

