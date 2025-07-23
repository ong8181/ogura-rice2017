####
#### Rice volatile data
####

# Load library
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2024.10.31
library(patchwork); packageVersion("patchwork") # 1.3.0, 2024.10.31
library(cowplot); packageVersion("cowplot") # 1.1.3, 2024.10.31
library(ggforce); packageVersion("ggforce") # 0.4.2, 2024.11.01

# Create output folder
dir.create("00_SessionInfo")
(outdir <- macam::outdir_create())


# ------------------------------------- #
# Load data and compile data
# ------------------------------------- #
d <- read.csv("data/data_insect.csv") %>% filter(treatment != "NC")
# Edit class
d$date <- ymd(d$date)
d$treatment[d$treatment == "Conv"] <- "Conventional"
d$treatment[d$treatment == "NoFert"] <- "No Fertilizer"

# Calculate mean
d_mean <- d %>%
  group_by(date, treatment) %>%
  summarize(other_planthopper = mean(other_planthopper),
            cicadellidae_spp = mean(cicadellidae_spp),
            thysanoptera_spp = mean(thysanoptera_spp))
d_mean_long <- d_mean %>% 
  pivot_longer(cols = -c(date, treatment))

# Make data longer
d_ist <- d %>%
  select(-Insect_sheet_collected) %>% 
  select(-sample_code, -n_week) %>% 
  pivot_longer(cols = -c(date, treatment, sampling_point))


# ------------------------------------- #
# Insect community (jitter + line)
# ------------------------------------- #
# Insect community
g1 <- d_ist %>%
  ggplot(aes(x = date, y = value, color = treatment, shape = treatment)) +
  geom_line(data = d_mean_long, aes(x = date, y = value, color = treatment)) +
  geom_jitter(height = 0, width = 0.5, alpha = 0.5) +
  facet_wrap(~ name, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  xlab(NULL) + ggtitle("Insect abundance") +
  ylab("Insect abundance (Ind.)")
g1

# ------------------------------------- #
# Combine all the figures
# ------------------------------------- #
ggsave("01_VisualizeInsectOut/Fig_InsectAbn.pdf", g1, height = 4, width = 10)
saveRDS(g1, "../07_FormatFigs/data_robj/Insect_temporal-pattern.obj")


# ------------------------------------- #
# Save workspace
# ------------------------------------- #
macam::save_session_info("00_SessionInfo")
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))


