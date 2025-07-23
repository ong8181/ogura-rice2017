####
#### Rice yield data
####

# Load library
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2024.10.31
library(patchwork); packageVersion("patchwork") # 1.3.0, 2024.10.31
library(cowplot); packageVersion("cowplot") # 1.1.3, 2024.10.31
library(ggforce); packageVersion("ggforce") # 0.4.2, 2024.11.01
theme_set(theme_gray())

# Create output folder
dir.create("00_SessionInfo")
(outdir <- macam::outdir_create())


# ------------------------------------- #
# Load data and compile data
# ------------------------------------- #
d <- read.csv("data/data_2017-ogura_rice-yield.csv")

# Compile data
d_hd <- d %>% # Total N of head
  select(treatment, total_head_12ind) %>% 
  pivot_longer(cols = - treatment)
d_hd$treatment[d_hd$treatment == "Conv"] <- "Conventional"
d_hd$treatment[d_hd$treatment == "NoFert"] <- "No Fertilizer"
d_wt <- d %>% # Total weight per 12 individuals
  select(treatment, total_wet_weight_12ind, total_dry_weight_12ind) %>% 
  pivot_longer(cols = - treatment)
d_wt$treatment[d_wt$treatment == "Conv"] <- "Conventional"
d_wt$treatment[d_wt$treatment == "NoFert"] <- "No Fertilizer"
d_ct <- d %>% # Total grain count per 10 heads
  select(treatment, total_grain_10head, total_fertile_10head, total_sterile_10heads) %>% 
  pivot_longer(cols = - treatment)
d_ct$treatment[d_ct$treatment == "Conv"] <- "Conventional"
d_ct$treatment[d_ct$treatment == "NoFert"] <- "No Fertilizer"
d_gn <- d %>% # Total grain weight per 10 heads
  select(treatment, total_grain_dry_g, fertile_grain_dry_g, sterile_grain_dry_g) %>% 
  pivot_longer(cols = - treatment)
d_gn$treatment[d_gn$treatment == "Conv"] <- "Conventional"
d_gn$treatment[d_gn$treatment == "NoFert"] <- "No Fertilizer"



# ------------------------------------- #
# Visualize data
# ------------------------------------- #
# Total number of rice heads
## Change facet labels
d_hd$name[d_hd$name == "total_head_12ind"] <- "Total heads"
## ggplot2
g1 <- d_hd %>%
  ggplot(aes(x = treatment, y = value, color = treatment)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(aes(shape = treatment), height = 0, width = 0.2, size = 2) +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2"), name = NULL) + 
  scale_shape_manual(values = c(16, 17), name = NULL) + 
  facet_wrap(~name) +
  xlab("Treatment") + ylab("Total number of rice heads") +
  scale_x_discrete(breaks=c("Conv","NoFert"), labels=c("Conventional", "No Fertilizer"))

# Total weight per 12 individuals
## Change facet labels
d_wt$name[d_wt$name == "total_dry_weight_12ind"] <- "Total dry weight"
d_wt$name[d_wt$name == "total_wet_weight_12ind"] <- "Total wet weight"
d_wt$name <- factor(d_wt$name, levels = c("Total wet weight", "Total dry weight"))
## ggplot2
g2 <- d_wt %>%
  ggplot(aes(x = treatment, y = value, color = treatment)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(aes(shape = treatment), height = 0, width = 0.2, size = 2) +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2"), name = NULL) + 
  scale_shape_manual(values = c(16, 17), name = NULL) + 
  facet_wrap(~name) + xlab("Treatment") + ylab("Weight (g/12 inds.)") +
  scale_x_discrete(breaks=c("Conv","NoFert"), labels=c("Conventional", "No Fertilizer"))

# Total grain count per 10 heads
## Change facet labels
d_ct$name[d_ct$name == "total_grain_10head"] <- "Total grain"
d_ct$name[d_ct$name == "total_fertile_10head"] <- "Fertile grain"
d_ct$name[d_ct$name == "total_sterile_10heads"] <- "Sterile grain"
d_ct$name <- factor(d_ct$name, levels = c("Total grain", "Fertile grain", "Sterile grain"))
## ggplot2
g3 <- d_ct %>%
  ggplot(aes(x = treatment, y = value, color = treatment)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(aes(shape = treatment), height = 0, width = 0.2, size = 2) +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2"), name = NULL) + 
  scale_shape_manual(values = c(16, 17), name = NULL) + 
  facet_wrap(~name) + xlab("Treatment") + ylab("Grain counts (/10 heads)") +
  scale_x_discrete(breaks=c("Conv","NoFert"), labels=c("Conventional", "No Fertilizer"))

# Total grain dry weight per 10 heads
## Change facet labels
d_gn$name[d_gn$name == "total_grain_dry_g"] <- "Total grain"
d_gn$name[d_gn$name == "fertile_grain_dry_g"] <- "Fertile grain"
d_gn$name[d_gn$name == "sterile_grain_dry_g"] <- "Sterile grain"
d_gn$name <- factor(d_gn$name, levels = c("Total grain", "Fertile grain", "Sterile grain"))
## ggplot2
g4 <- d_gn %>%
  ggplot(aes(x = treatment, y = value, color = treatment)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(aes(shape = treatment), height = 0, width = 0.2, size = 2) +
  scale_color_manual(values = c("darkgreen", "darkgoldenrod2"), name = NULL) + 
  scale_shape_manual(values = c(16, 17), name = NULL) + 
  facet_wrap(~name) + xlab("Treatment") + ylab("Dry wegith (g/10 heads)") +
  scale_x_discrete(breaks=c("Conv","NoFert"), labels=c("Conventional", "No Fertilizer"))


# ------------------------------------- #
# Combine all the figures
# ------------------------------------- #
g_all <- (g1 + g2 + plot_layout(widths = c(1,2.1))) / g3 / g4
ggsave("01_VisualizeYieldOut/Fig_RiceYield.pdf", g_all, height = 12, width = 10)
saveRDS(g_all, "../07_FormatFigs/data_robj/Rice_Yield.obj")
saveRDS(list(g1, g2, g3, g4), "../07_FormatFigs/data_robj/Rice_Yield_separate.obj")


# ------------------------------------- #
# Save workspace
# ------------------------------------- #
macam::save_session_info("00_SessionInfo")
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))


