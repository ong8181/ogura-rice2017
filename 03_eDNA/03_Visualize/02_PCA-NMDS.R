####
#### eDNA data, dimension reduction
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
# Load compiled phyloseq data
# ------------------------------------- #
ps_all <- readRDS("../02_AllDataCompile/04_RemoveSameSpOut/ps_comb_filt.obj")
# Compile data
sample_data(ps_all)$date <- ymd(sample_data(ps_all)$date)
sample_data(ps_all)$treatment[sample_data(ps_all)$treatment == "NF"] <- "No Fertilizer"
sample_data(ps_all)$treatment[sample_data(ps_all)$treatment == "CN"] <- "Conventional"
ps_sample <- ps_all %>% subset_samples(sample_nc == "sample")


# ------------------------------------- #
# Visualize data (PCA and NMDS)
# ------------------------------------- #
set.seed(1234)
# Absolute abundance
# NMDS
ps_bray1 <- ordinate(ps_sample, "NMDS", "bray", trymax = 100)
# Customize NMDS for plant volatile, absolute value
g1 <- plot_ordination(ps_sample, ps_bray1, color = "treatment")

# Statistical test (adonis2)
pd_com_df <- ps_sample %>% otu_table %>% data.frame
meta_df <- ps_sample %>% sample_data %>% data.frame
meta_df$treatment <- factor(meta_df$treatment)
ad_res <- adonis2(pd_com_df ~ Week * treatment, data = meta_df, method = "bray", permutations = 9999, by="terms")
ad_res

# Customize data
g_df <- data.frame(sample_data(ps_sample)) %>% 
  mutate(NMDS1 = g1$data$NMDS1,
         NMDS2 = g1$data$NMDS2)
g2 <- g_df %>% ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(alpha = Week, fill = treatment, shape = treatment), size = 3) +
  geom_point(data = g_df %>% filter(treatment == "Conventional"), aes(x = NMDS1, y = NMDS2), shape = 21, color = "black", size = 3) +
  geom_point(data = g_df %>% filter(treatment == "No Fertilizer"), aes(x = NMDS1, y = NMDS2), shape = 24, color = "black", size = 3) +
  stat_ellipse(geom = "polygon", aes(color = treatment, fill = treatment), alpha = 0.05) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_fill_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_shape_manual(values = c(21, 24), name = NULL) +
  scale_alpha(name = "N of weeks") +
  ggtitle("eDNA copy numbers, NMDS, Bray-Curtis")


# Relative abundance
ps_sample_rel <- ps_sample %>% transform_sample_counts(function(x) x/sum(x))
# NMDS
ps_bray2 <- ordinate(ps_sample_rel, "NMDS", "bray", trymax = 100)
# Customize NMDS for plant volatile, absolute value
g3 <- plot_ordination(ps_sample_rel, ps_bray2, color = "treatment")
# Customize data
g_df2 <- data.frame(sample_data(ps_sample_rel)) %>% 
  mutate(NMDS1 = g3$data$NMDS1,
         NMDS2 = g3$data$NMDS2)
g4 <- g_df2 %>% ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(alpha = Week, fill = treatment, shape = treatment), size = 3) +
  geom_point(data = g_df2 %>% filter(treatment == "Conventional"), aes(x = NMDS1, y = NMDS2), shape = 21, color = "black", size = 3) +
  geom_point(data = g_df2 %>% filter(treatment == "No Fertilizer"), aes(x = NMDS1, y = NMDS2), shape = 24, color = "black", size = 3) +
  stat_ellipse(geom = "polygon", aes(color = treatment, fill = treatment), alpha = 0.05) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_fill_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_shape_manual(values = c(21, 24), name = NULL) +
  scale_alpha(name = "N of weeks") +
  ggtitle("eDNA relative abundance, NMDS, Bray-Curtis")

# Statistical test (adonis2)
pd_com_df_rel <- ps_sample_rel %>% otu_table %>% data.frame
meta_df_rel <- ps_sample_rel %>% sample_data %>% data.frame
meta_df_rel$treatment <- factor(meta_df_rel$treatment)
ad_res_rel <- adonis2(pd_com_df_rel ~ Week * treatment, data = meta_df_rel, method = "bray", permutations = 9999, by="terms")
ad_res_rel


# ------------------------------------- #
# Check influential species
# ------------------------------------- #
# Extract species score (NOT reliable)
sp_scores <- scores(ps_bray1, display="species") %>% data.frame
sp_scores$len <- sqrt(sp_scores$NMDS1^2 + sp_scores$NMDS2^2)
## Add taxa information
sp_scores <- cbind(sp_scores, ps_sample %>% tax_table %>% data.frame)
sp_scores <- sp_scores[order(sp_scores$len, decreasing = TRUE),]

# EnvFit method
# Extract NMDS scores
site_scores <- scores(ps_bray1, display = "sites")
# Evaluate species contribution
fit_species <- envfit(site_scores, pd_com_df, permutations = 999)
fit_df <- fit_species$vectors$arrows %>% data.frame
fit_df$r2 <- fit_species$vectors$r
fit_df$pval <- fit_species$vectors$pvals
fit_df <- cbind(fit_df, ps_sample %>% tax_table %>% data.frame)
fit_df <- fit_df[order(fit_df$r2, decreasing = TRUE),]

# Extract influential species based on the envfit
top10_sp <- fit_df %>% head(12) %>% rownames
ps_top10 <- prune_taxa(taxa_names(ps_sample) %in% top10_sp, ps_sample)
## Convert the data format
ps_top10_melt <- psmelt(ps_top10)
ps_top10_melt$OTU <- factor(ps_top10_melt$OTU, levels = top10_sp)
ps_mt <- ps_top10_melt %>%
  group_by(date, Week, treatment, OTU) %>% 
  summarize(mean_abundance = mean(Abundance))
t1 <- ps_top10_melt %>%
  ggplot(aes(x = date, y = Abundance + 0.5, color = treatment)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.8) +
  geom_line(data = ps_mt, aes(x = date, y = mean_abundance+0.5)) +
  facet_wrap(~ OTU, scales = "free_y") +
  scale_y_log10(label= macam::label_10_to_power) +
  labs(x = NULL, y = "Log_10(eDNA copy numbers/ml water + 0.5)")


# ------------------------------------- #
# Combine all the figures
# ------------------------------------- #
ggsave("02_PCA-NMDSOut/Fig_DimReduc.pdf", g2, height = 6, width = 8)
saveRDS(list(g2, g4), "../../07_FormatFigs/data_robj/eDNA_DimRed.obj")
ggsave("02_PCA-NMDSOut/Fig_TopTaxa.pdf", t1, height = 8, width = 12)
saveRDS(t1, "../../07_FormatFigs/data_robj/eDNA_TopTaxa.obj")
write.csv(sp_scores, sprintf("%s/species_scores.csv", outdir), row.names = TRUE)
write.csv(fit_df, sprintf("%s/envfit_res.csv", outdir), row.names = TRUE)


# ------------------------------------- #
# Save workspace
# ------------------------------------- #
macam::save_session_info("00_SessionInfo")
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))


