####
#### Rice volatile data; multivariate analysis
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


# ------------------------------------- #
# Import data to phyloseq
# ------------------------------------- #
# Add rownames and colnames
rownames(d1) <- sprintf("S%03d", 1:nrow(d1))
rownames(d2) <- sprintf("S%03d", 1:nrow(d2))

# Omit samples with NA
d1_nona <- na.omit(d1)
d2_nona <- na.omit(d2)

# Import to phyloseq
ps_d1 <- phyloseq(sample_data(d1_nona[,c("date", "n_week", "treatment", "sampling_point", "sample_code")]),
                  otu_table(d1_nona[,6:12], taxa_are_rows = FALSE))
ps_d2 <- phyloseq(sample_data(d2_nona[,c("date", "n_week", "treatment", "sampling_point", "sample_code")]),
                  otu_table(d2_nona[,6:10], taxa_are_rows = FALSE))
# Remove "0" samples
ps_d2 <- prune_samples(sample_sums(ps_d2) > 0, ps_d2)

# Convert to the relative abundance
ps_d1_rel <- transform_sample_counts(ps_d1, function(x) x/sum(x))
ps_d2_rel <- transform_sample_counts(ps_d2, function(x) x/sum(x))


# ------------------------------------- #
# Visualize data (PCA and NMDS)
# ------------------------------------- #
set.seed(1234)
# Absolute abundance
# PCA
ps_euc_d1 <- ordinate(ps_d1, "PCoA", "euclidean")
ps_euc_d2 <- ordinate(ps_d2, "PCoA", "euclidean")
# Visualize
plot_ordination(ps_d1, ps_euc_d1, color = "treatment")
plot_ordination(ps_d2, ps_euc_d2, color = "treatment")
# NMDS
ps_bray_d1 <- ordinate(ps_d1, "NMDS", "bray", trymax = 100)
ps_bray_d2 <- ordinate(ps_d2, "NMDS", "bray", trymax = 100)
# Visualize
plot_ordination(ps_d1, ps_bray_d1, color = "treatment")
plot_ordination(ps_d2, ps_bray_d2, color = "treatment")

# Statistical test (adonis2)
pd_com_df <- ps_d1 %>% otu_table %>% data.frame
meta_df <- ps_d1 %>% sample_data %>% data.frame
meta_df$treatment <- factor(meta_df$treatment)
ad_res <- adonis2(pd_com_df ~ n_week * treatment, data = meta_df, method = "bray", permutations = 9999, by="terms")
ad_res

# Relative abundance
# PCA
ps_euc_d1r <- ordinate(ps_d1_rel, "PCoA", "euclidean")
ps_euc_d2r <- ordinate(ps_d2_rel, "PCoA", "euclidean")
# Visualize
plot_ordination(ps_d1_rel, ps_euc_d1r, color = "treatment")
plot_ordination(ps_d2_rel, ps_euc_d2r, color = "treatment")
# NMDS
ps_bray_d1r <- ordinate(ps_d1_rel, "NMDS", "bray", trymax = 100)
ps_bray_d2r <- ordinate(ps_d2_rel, "NMDS", "bray", trymax = 100)
# Visualize
plot_ordination(ps_d1_rel, ps_bray_d1r, color = "treatment")
plot_ordination(ps_d2_rel, ps_bray_d2r, color = "treatment")

# Customize NMDS for plant volatile, absolute value
g1_0 <- plot_ordination(ps_d1, ps_bray_d1, color = "treatment")
g1_df <- data.frame(sample_data(ps_d1)) %>% 
  mutate(NMDS1 = g1_0$data$NMDS1,
         NMDS2 = g1_0$data$NMDS2)

g1 <- g1_df %>% ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(alpha = n_week, fill = treatment, shape = treatment), size = 3) +
  geom_point(data = g1_df %>% filter(treatment == "Conventional"), aes(x = NMDS1, y = NMDS2), shape = 21, color = "black", size = 3) +
  geom_point(data = g1_df %>% filter(treatment == "No Fertilizer"), aes(x = NMDS1, y = NMDS2), shape = 24, color = "black", size = 3) +
  stat_ellipse(geom = "polygon", aes(color = treatment, fill = treatment), alpha = 0.05) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_fill_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_shape_manual(values = c(21, 24), name = NULL) +
  scale_alpha(name = "N of weeks") +
  ggtitle("Plant volatile (abs. conc.), NMDS, Bray-Curtis")

g2_0 <- plot_ordination(ps_d2, ps_euc_d2, color = "treatment")
g2_df <- data.frame(sample_data(ps_d2)) %>% 
  mutate(PC1 = g2_0$data$Axis.1,
         PC2 = g2_0$data$Axis.2)
g2 <- g2_df %>% ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(alpha = n_week, fill = treatment, shape = treatment), size = 3) +
  geom_point(data = g2_df %>% filter(treatment == "Conventional"), aes(x = PC1, y = PC2), shape = 21, color = "black", size = 3) +
  geom_point(data = g2_df %>% filter(treatment == "No Fertilizer"), aes(x = PC1, y = PC2), shape = 24, color = "black", size = 3) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_fill_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_shape_manual(values = c(21, 24), name = NULL) +
  scale_alpha(name = "N of weeks") +
  xlab("PC1 (90.8%)") + ylab("PC2 (8.8%)") +
  ggtitle("Phytoalexin (abs. conc.), PCA")


# ------------------------------------- #
# Combine all the figures
# ------------------------------------- #
g_all <- g1 + g2
ggsave("02_PCA-NMDSOut/Fig_DimReduc.pdf", g_all, height = 6, width = 14)
saveRDS(list(g1, g2), "../07_FormatFigs/data_robj/SpecializedMetabolites_DimRed.obj")


# ------------------------------------- #
# Save workspace
# ------------------------------------- #
macam::save_session_info("00_SessionInfo")
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))


