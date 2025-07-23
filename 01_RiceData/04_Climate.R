####
#### Climate data (Kyoto)
####

# Download climate data
# https://www.data.jma.go.jp/risk/obsdl/index.php
# Kyoto City, mean/max/mix daily temperature, daily precipitation, and daily sunlight hours

# Load library
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2024.10.31
library(patchwork); packageVersion("patchwork") # 1.3.0, 2024.10.31
library(ggforce); packageVersion("ggforce") # 0.4.2, 2024.11.01
theme_set(theme_bw())

# Create output folder
dir.create("00_SessionInfo")
(outdir <- macam::outdir_create())


# ------------------------------------- #
# Load data and compile data
# ------------------------------------- #
d <- read.csv("data/data_2017-climate-kyoto.csv") 
d$date <- ymd(d$date)
d <- d %>% filter(date >= "2017-05-01" & date <= "2017-09-30")


# ------------------------------------- #
# Visualize data
# ------------------------------------- #
g1 <- d %>% 
  ggplot(aes(x = date, y = mean_temp)) +
  geom_line() +
  geom_line(aes(x = date, y = max_temp), linetype = 2, color = "gray20") +
  geom_line(aes(x = date, y = min_temp), linetype = 2, color = "gray20") +
  ylab(expression(paste("Mean daily temperature (", degree, "C)"))) +
  xlab(NULL)

g2 <- d %>% 
  ggplot(aes(x = date, y = daily_rainfall)) +
  geom_line() +
  ylab("Daily rainfall (mm)") + xlab(NULL)

g3 <- d %>% 
  ggplot(aes(x = date, y = sunlight_hours)) +
  geom_line() +
  ylab("Daily sunlight (hours)") + xlab(NULL)



# ------------------------------------- #
# Combine all the figures
# ------------------------------------- #
g_all <- (g1 / g2 / g3)
ggsave("04_ClimateOut/Fig_Climate.pdf", g_all, height = 8, width = 10)
saveRDS(list(g1, g2, g3), "../07_FormatFigs/data_robj/Climate_all.obj")


# ------------------------------------- #
# Save workspace
# ------------------------------------- #
macam::save_session_info("00_SessionInfo")
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))


