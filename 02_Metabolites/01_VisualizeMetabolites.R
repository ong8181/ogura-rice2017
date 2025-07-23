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
d1 <- read.csv("data/data_volatile.csv") %>% filter(treatment != "NC")
d2 <- read.csv("data/data_phytoalexin.csv") %>% filter(treatment != "NC")

# Edit variables
## d1
d1$date <- ymd(d1$date)
d1$treatment[d1$treatment == "Conv"] <- "Conventional"
d1$treatment[d1$treatment == "NoFert"] <- "No Fertilizer"
## d2
d2$date <- ymd(d2$date)
d2$treatment[d2$treatment == "Conv"] <- "Conventional"
d2$treatment[d2$treatment == "NoFert"] <- "No Fertilizer"
d2 <- d2 %>% mutate(diterpenoids_ug_g_FW = momilactone_A_ug_g_FW + momilactone_B_ug_g_FW +
                      phytocassane_B_ug_g_FW + phytocassane_B3_ug_g_FW)

# Calculate mean
d1_mean <- d1 %>%
  group_by(date, treatment) %>%
  summarize(b_caryophyllene = mean(b_caryophyllene, na.rm = TRUE),
            E_2_Hexenal = mean(E_2_Hexenal, na.rm = TRUE),
            linalool = mean(linalool, na.rm = TRUE),
            MeSA = mean(MeSA, na.rm = TRUE),
            X2_4_Heptadienal = mean(X2_4_Heptadienal, na.rm = TRUE),
            Z_3_Hexen_1_ol = mean(Z_3_Hexen_1_ol, na.rm = TRUE),
            Z_3_Hexenal = mean(Z_3_Hexenal, na.rm = TRUE))
d2_mean <- d2 %>%
  group_by(date, treatment) %>%
  summarize(momilactone_A_ug_g_FW = mean(momilactone_A_ug_g_FW, na.rm = TRUE),
            momilactone_B_ug_g_FW = mean(momilactone_B_ug_g_FW, na.rm = TRUE),
            phytocassane_B_ug_g_FW = mean(phytocassane_B_ug_g_FW, na.rm = TRUE),
            phytocassane_B3_ug_g_FW = mean(phytocassane_B3_ug_g_FW, na.rm = TRUE),
            sakuranetin_ug_g_FW = mean(sakuranetin_ug_g_FW, na.rm = TRUE),
            diterpenoids_ug_g_FW = mean(diterpenoids_ug_g_FW, na.rm = TRUE))
d1_mean_long <- d1_mean %>% 
  pivot_longer(cols = -c(date, treatment))
d2_mean_long <- d2_mean %>% 
  pivot_longer(cols = -c(date, treatment))

# Make data longer
d_vl <- d1 %>%
  select(-sample_code, -n_week) %>% 
  pivot_longer(cols = -c(date, treatment, sampling_point))
d_pt <- d2 %>%
  select(-sample_code, -n_week) %>% 
  pivot_longer(cols = -c(date, treatment, sampling_point))

# Combine two datasets
## Check sample_code
d3 <- d1
d3_mean <- d1_mean
all(d3$sample_code == d2$sample_code)
all(paste(d3_mean$date, d3_mean$treatment) == paste(d2_mean$date, d2_mean$treatment))
## Add diterpenoids
d3$diterpenoids_ug_g_FW <- d2$diterpenoids_ug_g_FW
d3_mean$diterpenoids_ug_g_FW <- d2_mean$diterpenoids_ug_g_FW
d_d3 <- d3 %>%
  select(-sample_code, -n_week) %>% 
  pivot_longer(cols = -c(date, treatment, sampling_point))
d3_mean_long <- d3_mean %>% 
  pivot_longer(cols = -c(date, treatment))
## Change levels and names
mtbl_levels <- d_d3$name %>% unique
d_d3$name <- factor(d_d3$name, levels = mtbl_levels)
d3_mean_long$name <- factor(d3_mean_long$name, levels = mtbl_levels)

# ------------------------------------- #
# Visualize data (jitter + line)
# ------------------------------------- #
# Volatile
g1 <- d_vl %>%
  ggplot(aes(x = date, y = value, color = treatment, shape = treatment)) +
  geom_line(data = d1_mean_long, aes(x = date, y = value, color = treatment)) +
  geom_jitter(height = 0, width = 0.5, alpha = 0.5) +
  facet_wrap(~ name, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  xlab(NULL) + ggtitle("Plant organic volatile") +
  ylab("Volatile concentration (relative peak area / leaf mass)")

# Phytoalexin
g2 <- d_pt %>%
  ggplot(aes(x = date, y = value, color = treatment, shape = treatment)) +
  geom_line(data = d2_mean_long, aes(x = date, y = value, color = treatment)) +
  geom_jitter(height = 0, width = 0.5, alpha = 0.5) +
  facet_wrap(~ name, scales = "free_y", nrow = 2) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  xlab(NULL) + ggtitle("Phytoalexin") +
  ylab(expression(paste("Phytoalexin content (", mu, "g/g)")))

# Combined
g3 <- d_d3 %>%
  ggplot(aes(x = date, y = value, color = treatment, shape = treatment)) +
  geom_line(data = d3_mean_long, aes(x = date, y = value, color = treatment)) +
  geom_jitter(height = 0, width = 0.5, alpha = 0.5) +
  facet_wrap(~ name, scales = "free_y", nrow = 2) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  xlab(NULL) + ggtitle("Plant specialized metabolites") +
  ylab("Metabolite concentration")

# ------------------------------------- #
# Combine all the figures
# ------------------------------------- #
g_all <- (g1 / g2) + plot_layout(heights = c(1.55, 1))
ggsave("01_VisualizeMetabolitesOut/Fig_VolitilePhytoalexin.pdf", g_all, height = 12, width = 10)
ggsave("01_VisualizeMetabolitesOut/Fig_SpecializedMetabolites.pdf", g3, height = 5, width = 11)
saveRDS(list(g1, g2, g3), "../07_FormatFigs/data_robj/SpecializedMetabolites.obj")


# ------------------------------------- #
# Save workspace
# ------------------------------------- #
macam::save_session_info("00_SessionInfo")
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))


