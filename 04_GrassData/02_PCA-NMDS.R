####
#### Plant-coverage data; Dimension reduction
####

# Load library
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2024.10.31
library(patchwork); packageVersion("patchwork") # 1.3.0, 2024.10.31
library(cowplot); packageVersion("cowplot") # 1.1.3, 2024.10.31
library(vegan); packageVersion("vegan") # 2.6.8, 2024.11.5, PCA and NMDS
library(phyloseq); packageVersion("phyloseq") # 1.46.0, 2024.11.5, PCA and NMDS

# Create output folder
dir.create("00_SessionInfo")
(outdir <- macam::outdir_create())


# ------------------------------------- #
# Load data and compile data
# ------------------------------------- #
d <- read.csv("data/data_grass.csv") %>% filter(treatment != "NC")
# Edit class
d$date <- ymd(d$date)
d$treatment[d$treatment == "Conv"] <- "Conventional"
d$treatment[d$treatment == "NoFert"] <- "No Fertilizer"
# Remove dates when grass pics were not taken
d <- d %>%
  filter(grass_pic_taken == "Y") %>%
  replace(is.na(.), 0)


# ------------------------------------- #
# Import data to phyloseq
# ------------------------------------- #
# Add rownames and colnames
rownames(d) <- sprintf("S%03d", 1:nrow(d))
# Omit samples with NA
d_nona <- na.omit(d)

# Import to phyloseq
ps_d <- phyloseq(sample_data(d_nona[,c("date", "n_week", "treatment", "sampling_point", "sample_code")]),
                  otu_table(d_nona[,7:21], taxa_are_rows = FALSE))
# Remove "0" samples
ps_d <- prune_samples(sample_sums(ps_d) > 0, ps_d)


# ------------------------------------- #
# Visualize data (PCA and NMDS)
# ------------------------------------- #
set.seed(1234)
# Absolute abundance
# PCA
ps_euc_d <- ordinate(ps_d, "PCoA", "euclidean")
# Visualize
plot_ordination(ps_d, ps_euc_d, color = "treatment")
# NMDS
ps_bray_d <- ordinate(ps_d, "NMDS", "bray", trymax = 100)
# Visualize
plot_ordination(ps_d, ps_bray_d, color = "treatment")

# Statistical test (adonis2)
pd_com_df <- ps_d %>% otu_table %>% data.frame
meta_df <- ps_d %>% sample_data %>% data.frame
meta_df$treatment <- factor(meta_df$treatment)
ad_res <- adonis2(pd_com_df ~ n_week * treatment, data = meta_df, method = "bray", permutations = 9999, by="terms")
ad_res

# Customize NMDS for plant volatile, absolute value
g1 <- plot_ordination(ps_d, ps_bray_d, color = "treatment")

# Customize data
g_df <- data.frame(sample_data(ps_d)) %>% 
  mutate(NMDS1 = g1$data$NMDS1,
         NMDS2 = g1$data$NMDS2)
g2 <- g_df %>% ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(alpha = n_week, fill = treatment, shape = treatment), size = 3) +
  geom_point(data = g_df %>% filter(treatment == "Conventional"), aes(x = NMDS1, y = NMDS2), shape = 21, color = "black", size = 3) +
  geom_point(data = g_df %>% filter(treatment == "No Fertilizer"), aes(x = NMDS1, y = NMDS2), shape = 24, color = "black", size = 3) +
  stat_ellipse(geom = "polygon", aes(color = treatment, fill = treatment), alpha = 0.05) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_fill_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_shape_manual(values = c(21, 24), name = NULL) +
  scale_alpha(name = "N of weeks") +
  ggtitle("Grass coverage (%), NMDS, Bray-Curtis")


# ------------------------------------- #
# Combine all the figures
# ------------------------------------- #
ggsave("02_PCA-NMDSOut/Fig_DimReduc.pdf", g2, height = 6, width = 8)
saveRDS(list(g1, g2), "../07_FormatFigs/data_robj/Grass_DimRed.obj")


# ------------------------------------- #
# Save workspace
# ------------------------------------- #
macam::save_session_info("00_SessionInfo")
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))


