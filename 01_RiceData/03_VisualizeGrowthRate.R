####
#### Rice growth data
####

# Load library
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2024.10.31
library(patchwork); packageVersion("patchwork") # 1.3.0, 2024.10.31
library(ggforce); packageVersion("ggforce") # 0.4.2, 2024.11.01
theme_set(theme_gray())

# Create output folder
dir.create("00_SessionInfo")
(outdir <- macam::outdir_create())


# ------------------------------------- #
# Load data and compile data
# ------------------------------------- #
d <- read.csv("data/data_2017-ogura_rice-monitoring.csv") %>% 
  filter(treatment != "NC") %>% 
  select(date, time, n_week, treatment, sampling_point, n_stem_1, n_stem_2, n_stem_3,
         spad_1, spad_2, spad_3, height_1, height_2, height_3,
         ph, ec, tds, temperature, water_level, comment)
d$date <- ymd(d$date)
d$treatment[d$treatment == "Conv"] <- "Conventional"
d$treatment[d$treatment == "NoFert"] <- "No Fertilizer"

## Add mean SPAD, n_stem, and height, and 1st differences
d$n_stem_diff_1 <- d$n_stem_1 - lag(d$n_stem_1, n = 10) # Change rate
d$n_stem_diff_2 <- d$n_stem_2 - lag(d$n_stem_2, n = 10) # Change rate
d$n_stem_diff_3 <- d$n_stem_3 - lag(d$n_stem_3, n = 10) # Change rate
d$spad_diff_1 <- d$spad_1 - lag(d$spad_1, n = 10) # Change rate
d$spad_diff_2 <- d$spad_2 - lag(d$spad_2, n = 10) # Change rate
d$spad_diff_3 <- d$spad_3 - lag(d$spad_3, n = 10) # Change rate
d$height_diff_1 <- d$height_1 - lag(d$height_1, n = 10) # Change rate
d$height_diff_2 <- d$height_2 - lag(d$height_2, n = 10) # Change rate
d$height_diff_3 <- d$height_3 - lag(d$height_3, n = 10) # Change rate
# Mean change rates
d$mean_n_stem <- apply(d[,c("n_stem_1","n_stem_2","n_stem_3")], 1, function(x) mean(x, na.rm = TRUE))
d$mean_spad <- apply(d[,c("spad_1","spad_2","spad_3")], 1, function(x) mean(x, na.rm = TRUE))
d$mean_height <- apply(d[,c("height_1","height_2","height_3")], 1, function(x) mean(x, na.rm = TRUE))
d$mean_n_stem_diff <- d$mean_n_stem - lag(d$mean_n_stem, n = 10) # Change rate
d$mean_spad_diff <- d$mean_spad - lag(d$mean_spad, n = 10) # Change rate
d$mean_height_diff <- d$mean_height - lag(d$mean_height, n = 10) # Change rate

# Calculate mean
d_mean <- d %>% group_by(date, treatment) %>% summarize(n_stem_diff = mean(mean_n_stem_diff))
d_mean$spad_diff <- d %>% group_by(date, treatment) %>% summarize(spad = mean(mean_spad_diff)) %>% pull(spad)
d_mean$height_diff <- d %>% group_by(date, treatment) %>% summarize(height = mean(mean_height_diff)) %>% pull(height)

# Compile data
d_nst <- d %>% # The number of stems
  select(date, time, treatment, sampling_point, starts_with("n_stem_diff")) %>% 
  pivot_longer(cols = c(starts_with("n_stem_diff")),
               names_to = "ind_id", values_to = "n_stem_diff")
d_nst$ind_id[d_nst$ind_id == "n_stem_diff_1"] <- "Ind_1"
d_nst$ind_id[d_nst$ind_id == "n_stem_diff_2"] <- "Ind_2"
d_nst$ind_id[d_nst$ind_id == "n_stem_diff_3"] <- "Ind_3"

d_spd <- d %>% # SPAD
  select(date, time, treatment, sampling_point, starts_with("spad_diff")) %>% 
  pivot_longer(cols = starts_with("spad_diff"),
               names_to = "ind_id", values_to = "spad_diff")
d_spd$ind_id[d_spd$ind_id == "spad_diff_1"] <- "Ind_1"
d_spd$ind_id[d_spd$ind_id == "spad_diff_2"] <- "Ind_2"
d_spd$ind_id[d_spd$ind_id == "spad_diff_3"] <- "Ind_3"

d_hgt <- d %>% # Rice height
  select(date, time, treatment, sampling_point, starts_with("height_diff")) %>% 
  pivot_longer(cols = starts_with("height_diff"),
               names_to = "ind_id", values_to = "height_diff")
d_hgt$ind_id[d_hgt$ind_id == "height_diff_1"] <- "Ind_1"
d_hgt$ind_id[d_hgt$ind_id == "height_diff_2"] <- "Ind_2"
d_hgt$ind_id[d_hgt$ind_id == "height_diff_3"] <- "Ind_3"


# ------------------------------------- #
# Visualize data
# ------------------------------------- #
# Rice height
g1 <- d_hgt %>% 
  ggplot(aes(x = date, y = height_diff, color = treatment, shape = treatment)) +
  geom_line(data = d_mean, aes(x = date, y = height_diff, group = treatment)) +
  geom_jitter(height = 0, width = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  xlab("Date") + ylab("Changes in rice height (cm)")

# SPAD
g2 <- d_spd %>% 
  ggplot(aes(x = date, y = spad_diff, color = treatment, shape = treatment)) +
  geom_line(data = d_mean, aes(x = date, y = spad_diff, group = treatment)) +
  geom_jitter(height = 0, width = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  xlab("Date") + ylab("Changes in SPAD")

# The number of stems
g3 <- d_nst %>% 
  ggplot(aes(x = date, y = n_stem_diff, color = treatment, shape = treatment)) +
  geom_line(data = d_mean, aes(x = date, y = n_stem_diff, group = treatment)) +
  geom_jitter(height = 0, width = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  xlab("Date") + ylab("Changes in the number of stems")


# ------------------------------------- #
# Combine all the figures
# ------------------------------------- #
g_all <- (g1 / g2 / g3)
ggsave("03_VisualizeGrowthRateOut/Fig_RiceGrowthRate.pdf", g_all, height = 8, width = 10)
saveRDS(list(g1, g2, g3), "../07_FormatFigs/data_robj/Rice_Growth_diff.obj")


# ------------------------------------- #
# Save workspace
# ------------------------------------- #
macam::save_session_info("00_SessionInfo")
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))


