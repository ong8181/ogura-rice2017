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

# Calculate mean
d$mean_n_stem <- apply(d %>% select(n_stem_1, n_stem_2, n_stem_3), 1, function(x) mean(x, rm.na = T))
d$mean_spad <- apply(d %>% select(spad_1, spad_2, spad_3), 1, function(x) mean(x, rm.na = T))
d$mean_height <- apply(d %>% select(height_1, height_2, height_3), 1, function(x) mean(x, rm.na = T))
d_mean <- d %>% group_by(date, treatment) %>% summarize(n_stem = mean(mean_n_stem))
d_mean$spad <- d %>% group_by(date, treatment) %>% summarize(spad = mean(mean_spad)) %>% pull(spad)
d_mean$height <- d %>% group_by(date, treatment) %>% summarize(height = mean(mean_height)) %>% pull(height)

# Compile data
d_nst <- d %>% # The number of stems
  select(date, time, treatment, sampling_point, starts_with("n_stem")) %>% 
  pivot_longer(cols = starts_with("n_stem"),
               names_to = "ind_id", values_to = "n_stem")
d_nst$ind_id[d_nst$ind_id == "n_stem_1"] <- "Ind_1"
d_nst$ind_id[d_nst$ind_id == "n_stem_2"] <- "Ind_2"
d_nst$ind_id[d_nst$ind_id == "n_stem_3"] <- "Ind_3"

d_spd <- d %>% # SPAD
  select(date, time, treatment, sampling_point, starts_with("spad")) %>% 
  pivot_longer(cols = starts_with("spad"),
               names_to = "ind_id", values_to = "spad")
d_spd$ind_id[d_spd$ind_id == "spad_1"] <- "Ind_1"
d_spd$ind_id[d_spd$ind_id == "spad_3"] <- "Ind_3"
d_spd$ind_id[d_spd$ind_id == "spad_2"] <- "Ind_2"

d_hgt <- d %>% # Rice height
  select(date, time, treatment, sampling_point, starts_with("height")) %>% 
  pivot_longer(cols = starts_with("height"),
               names_to = "ind_id", values_to = "height")
d_hgt$ind_id[d_hgt$ind_id == "height_1"] <- "Ind_1"
d_hgt$ind_id[d_hgt$ind_id == "height_2"] <- "Ind_2"
d_hgt$ind_id[d_hgt$ind_id == "height_3"] <- "Ind_3"


# ------------------------------------- #
# Visualize data
# ------------------------------------- #
# Rice height
g1 <- d_hgt %>%
  ggplot(aes(x = date, y = height, color = treatment, shape = treatment)) +
  geom_line(data = d_mean, aes(x = date, y = height, group = treatment)) +
  #geom_sina(aes(group = date), alpha = 0.5, scale = "width", maxwidth = 1) +
  geom_jitter(height = 0, width = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  xlab("Date") + ylab("Rice height (cm)")

# SPAD
g2 <- d_spd %>% 
  ggplot(aes(x = date, y = spad, color = treatment, shape = treatment)) +
  geom_line(data = d_mean, aes(x = date, y = spad, group = treatment)) +
  #geom_sina(aes(group = date), alpha = 0.5) +
  geom_jitter(height = 0, width = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  xlab("Date") + ylab("SPAD")

# The number of stems
g3 <- d_nst %>% 
  ggplot(aes(x = date, y = n_stem, color = treatment, shape = treatment)) +
  geom_line(data = d_mean, aes(x = date, y = n_stem, group = treatment)) +
  geom_jitter(height = 0, width = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  xlab("Date") + ylab("The number of stems")


# ------------------------------------- #
# Combine all the figures
# ------------------------------------- #
g_all <- (g1 / g2 / g3)
ggsave("02_VisualizeGrowthOut/Fig_RiceGrowth.pdf", g_all, height = 8, width = 10)
saveRDS(list(g1, g2, g3), "../07_FormatFigs/data_robj/Rice_Growth_rawdata.obj")


# ------------------------------------- #
# Save workspace
# ------------------------------------- #
macam::save_session_info("00_SessionInfo")
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))


