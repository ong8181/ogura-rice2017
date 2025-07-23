####
#### Rice ecological community data
####

# Load library
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2024.10.31
library(phyloseq); packageVersion("phyloseq") # 1.46.0, 2024.11.05
library(patchwork); packageVersion("patchwork") # 1.3.0, 2024.10.31
library(cowplot); packageVersion("cowplot") # 1.1.3, 2024.10.31
library(ggforce); packageVersion("ggforce") # 0.4.2, 2024.11.01
library(macam); packageVersion("macam") # 0.1.5
library(cols4all); packageVersion("cols4all") # 0.8

# Create output folder
dir.create("00_SessionInfo")
(outdir <- macam::outdir_create())


# ------------------------------------- #
# Load compiled phyloseq data
# ------------------------------------- #
ps_all <- readRDS("../02_AllDataCompile/04_RemoveSameSpOut/ps_comb_filt.obj")
# Compile date
sample_data(ps_all)$date <- ymd(sample_data(ps_all)$date)
sample_data(ps_all)$treatment[sample_data(ps_all)$treatment == "NF"] <- "No Fertilizer"
sample_data(ps_all)$treatment[sample_data(ps_all)$treatment == "CN"] <- "Conventional"
ps_sample <- ps_all %>% subset_samples(sample_nc == "sample")


# ------------------------------------- #
# Compile taxa information
# ------------------------------------- #
# Add representative-taxa name
new_tax_df <- tax_table(ps_sample) %>%
  data.frame %>% 
  mutate(rep_tax = tax_table(ps_sample)[,"superkingdom"])

# Add "Fungi" category
new_tax_df$rep_tax[new_tax_df$kingdom == "Fungi"] <- "Fungi"

# Add "Other Eukaryota" and "Undetermined" categories
new_tax_df[new_tax_df$rep_tax == "Eukaryota","rep_tax"] <- "Non-Fungi Eukaryota"
new_tax_df[new_tax_df$rep_tax == "","rep_tax"] <- "Undetermined"
new_tax_df$rep_tax <- factor(new_tax_df$rep_tax,
                             levels = c("Archaea","Bacteria","Fungi","Non-Fungi Eukaryota","Undetermined"))
#unique(new_tax_df$rep_tax)
tax_table(ps_sample) <- as.matrix(new_tax_df)


# ------------------------------------- #
# Visualize
# ------------------------------------- #
# psmelt
ps_m1 <- speedyseq::psmelt(ps_sample) %>% 
  select(OTU, Sample, Abundance, date, Time, Week, plot, treatment, sample_nc, rep_tax) %>% 
  filter(Abundance > 0)

# Abundance (DNA copies/ml)
g1 <- ps_m1 %>% ggplot(aes(x = date, y = Abundance, fill = rep_tax)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 0) +
  ggsci::scale_fill_startrek(name = NULL) +
  facet_wrap(~treatment) +
  scale_y_continuous(label= label_10_to_power) +
  xlab(NULL) + ylab("Abundance (eDNA copies/ml water)") +
  NULL

# Diversity (# of ASVs)
## Count the number of ASVs
ps_m2 <- ps_m1 %>%
  group_by(date, Week, plot, treatment) %>% 
  summarize(n_asv = n())
ps_m3 <- ps_m1 %>%
  group_by(date, Week, plot, treatment) %>% 
  summarize(n_asv = n()) %>% 
  group_by(date, Week, treatment) %>% 
  summarize(mean_asv = mean(n_asv))

g2 <- ps_m2 %>% 
  ggplot(aes(x = date, y = n_asv, color = treatment, shape = treatment)) +
  geom_jitter(width = 0, alpha = 0.5) +
  geom_line(data = ps_m3, aes(x = date, y = mean_asv)) +
  #scale_fill_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  xlab(NULL) + ylab("No. of ASV") +
  ggtitle("No. of ASV ") +
  ylim(0,850) +
  NULL


# ------------------------------------- #
# Check top species
# ------------------------------------- #
# Overall top
top10 <- taxa_sums(ps_sample) %>% sort(decreasing = TRUE) %>% head(12) %>% names
ps_top <- prune_taxa(taxa_names(ps_sample) %in% top10, ps_sample)
ps_top_df <- ps_top %>% tax_table %>% data.frame
ps_top_df$total_abundance <- taxa_sums(ps_top)
ps_top_df <- ps_top_df[order(ps_top_df$total_abundance, decreasing = TRUE),]
# Non-bac top
top10_nonbac <- taxa_sums(ps_sample %>% subset_taxa(superkingdom != "Bacteria")) %>% sort(decreasing = TRUE) %>% head(12) %>% names
ps_top_nobac <- prune_taxa(taxa_names(ps_sample) %in% top10_nonbac, ps_sample)
ps_top_nobac_df <- ps_top_nobac %>% tax_table %>% data.frame
ps_top_nobac_df$total_abundance <- taxa_sums(ps_top_nobac)
ps_top_nobac_df <- ps_top_nobac_df[order(ps_top_nobac_df$total_abundance, decreasing = TRUE),]

## Convert the data format
ps_top_melt <- psmelt(ps_top)
ps_top_melt$OTU <- factor(ps_top_melt$OTU, levels = top10)
ps_m4 <- ps_top_melt %>%
  group_by(date, Week, treatment, OTU) %>% 
  summarize(mean_abundance = mean(Abundance))

t1 <- ps_top_melt %>%
  ggplot(aes(x = date, y = Abundance + 0.5, color = treatment)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.8) +
  geom_line(data = ps_m4, aes(x = date, y = mean_abundance+0.5)) +
  #stat_smooth(se = FALSE) +
  facet_wrap(~ OTU, scales = "free_y") +
  scale_y_log10(label= macam::label_10_to_power) +
  labs(x = NULL, y = "Abundance (eDNA copies /ml water + 0.5)")


# ------------------------------------- #
# Check specific groups
# ------------------------------------- #
## Visualize Fungi
ps_fungi <- subset_taxa(ps_sample, kingdom == "Fungi") %>% prune_taxa(taxa_sums(.)>0, .)
ps_fungi_order <- taxa_sums(ps_fungi) %>% sort(decreasing = TRUE) %>% names
ps_fungi_melt <- psmelt(ps_fungi)
ps_fungi_melt$OTU <- factor(ps_fungi_melt$OTU, levels = ps_fungi_order)
ps_a1 <- ps_fungi_melt %>%
  group_by(date, Week, treatment, plot, phylum) %>% 
  summarize(total_abundance = sum(Abundance))
ps_a2 <- ps_a1 %>%
  group_by(date, Week, treatment, phylum) %>% 
  summarize(mean_abundance = mean(total_abundance))
ps_a1$phylum[ps_a1$phylum == ""] <- "Undetermined"
ps_a2$phylum[ps_a2$phylum == ""] <- "Undetermined"
a1 <- ps_a1 %>%
  ggplot(aes(x = date, y = total_abundance + 0.5, color = treatment)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.8) +
  geom_line(data = ps_a2, aes(x = date, y = mean_abundance+0.5, color = treatment)) +
  facet_wrap(~ phylum, scales = "free_y", ncol = 4) +
  scale_y_log10(label= macam::label_10_to_power) +
  labs(x = NULL, y = "Abundance (eDNA copies/ml water + 0.5)",
       title = "Fungi")

## Visualize Bacteria
ps_bacte <- subset_taxa(ps_sample, superkingdom == "Bacteria") %>% prune_taxa(taxa_sums(.)>0, .)
ps_bacte_order <- taxa_sums(ps_bacte) %>% sort(decreasing = TRUE) %>% names
ps_bacte_melt <- psmelt(ps_bacte)
ps_bacte_melt$OTU <- factor(ps_bacte_melt$OTU, levels = ps_bacte_order)
ps_a3 <- ps_bacte_melt %>%
  group_by(date, Week, treatment, plot, phylum) %>% 
  summarize(total_abundance = sum(Abundance))
ps_a4 <- ps_a3 %>%
  group_by(date, Week, treatment, phylum) %>% 
  summarize(mean_abundance = mean(total_abundance))
ps_a3$phylum[ps_a3$phylum == ""] <- "Undetermined"
ps_a4$phylum[ps_a4$phylum == ""] <- "Undetermined"
ps_a3$phylum <- factor(ps_a3$phylum, levels = c(unique(ps_a3$phylum)[2:11], "Undetermined"))
ps_a4$phylum <- factor(ps_a4$phylum, levels = c(unique(ps_a4$phylum)[2:11], "Undetermined"))
a2 <- ps_a3 %>%
  ggplot(aes(x = date, y = total_abundance + 0.5, color = treatment)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.8) +
  geom_line(data = ps_a4, aes(x = date, y = mean_abundance+0.5, color = treatment)) +
  facet_wrap(~ phylum, scales = "free_y", ncol = 4) +
  scale_y_log10(label= macam::label_10_to_power) +
  labs(x = NULL, y = "Abundance (eDNA copies/ml water + 0.5)",
       title = "Bacteria")

## Visualize Non-Fungi eukaryota
ps_eukar <- subset_taxa(ps_sample, superkingdom == "Eukaryota" & kingdom != "Fungi") %>% prune_taxa(taxa_sums(.)>0, .)
ps_eukar_order <- taxa_sums(ps_eukar) %>% sort(decreasing = TRUE) %>% names
ps_eukar_melt <- psmelt(ps_eukar)
ps_eukar_melt$OTU <- factor(ps_eukar_melt$OTU, levels = ps_eukar_order)
ps_a5 <- ps_eukar_melt %>%
  group_by(date, Week, treatment, plot, phylum) %>% 
  summarize(total_abundance = sum(Abundance))
ps_a6 <- ps_a5 %>%
  group_by(date, Week, treatment, phylum) %>% 
  summarize(mean_abundance = mean(total_abundance))
ps_a5$phylum[ps_a5$phylum == ""] <- "Undetermined"
ps_a6$phylum[ps_a6$phylum == ""] <- "Undetermined"
ps_a5$phylum <- factor(ps_a5$phylum, levels = c(unique(ps_a5$phylum)[2:21], "Undetermined"))
ps_a6$phylum <- factor(ps_a6$phylum, levels = c(unique(ps_a6$phylum)[2:21], "Undetermined"))
a3 <- ps_a5 %>%
  ggplot(aes(x = date, y = total_abundance + 0.5, color = treatment)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.8) +
  geom_line(data = ps_a6, aes(x = date, y = mean_abundance+0.5, color = treatment)) +
  facet_wrap(~ phylum, scales = "free_y", ncol = 4) +
  scale_y_log10(label= macam::label_10_to_power) +
  labs(x = NULL, y = "Abundance (eDNA copies/ml water + 0.5)",
       title = "Non-Fungi Eukaryota")


# ------------------------------------- #
# Check causal species
# ------------------------------------- #
# Rice height growth
cause_spp1 <- c("ITS_Taxa00138","PRO_Taxa00019","PRO_Taxa00332","ITS_Taxa00007","PRO_Taxa00168","PRO_Taxa00186",
                "PRO_Taxa00037","PRO_Taxa00024","PRO_Taxa00040","PRO_Taxa00096","COI_Taxa00221",
                "EUK_Taxa00213","PRO_Taxa00046","EUK_Taxa00443","PRO_Taxa00017")
ps_cause1 <- prune_taxa(cause_spp1, ps_sample) %>% prune_taxa(taxa_sums(.) > 0, .)
ps_cause1_melt <- psmelt(ps_cause1)
ps_c1 <- ps_cause1_melt %>%
  group_by(date, Week, treatment, plot, OTU) %>% 
  summarize(total_abundance = sum(Abundance))
ps_c2 <- ps_c1 %>%
  group_by(date, Week, treatment, OTU) %>% 
  summarize(mean_abundance = mean(total_abundance))
ps_c1$OTU <- factor(ps_c1$OTU, levels = cause_spp1)
ps_c2$OTU <- factor(ps_c2$OTU, levels = cause_spp1)
c1 <- ps_c1 %>%
  ggplot(aes(x = date, y = total_abundance + 0.5, color = treatment)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.8) +
  geom_line(data = ps_c2, aes(x = date, y = mean_abundance+0.5, color = treatment)) +
  facet_wrap(~ OTU, scales = "free_y", ncol = 4) +
  scale_y_log10(label= macam::label_10_to_power) +
  labs(x = NULL, y = "Abundance (eDNA copies/ml water + 0.5)",
       title = "Causal taxa to rice height growth")

# SPAD
cause_spp2 <- c("COI_Taxa00127","EUK_Taxa00659","PRO_Taxa00050","COI_Taxa00134","EUK_Taxa00306","PRO_Taxa00013","PRO_Taxa00011")
ps_cause2 <- prune_taxa(cause_spp2, ps_sample) %>% prune_taxa(taxa_sums(.) > 0, .)
ps_cause2_melt <- psmelt(ps_cause2)
ps_c3 <- ps_cause2_melt %>%
  group_by(date, Week, treatment, plot, OTU) %>% 
  summarize(total_abundance = sum(Abundance))
ps_c4 <- ps_c3 %>%
  group_by(date, Week, treatment, OTU) %>% 
  summarize(mean_abundance = mean(total_abundance))
ps_c3$OTU <- factor(ps_c3$OTU, levels = cause_spp2)
ps_c4$OTU <- factor(ps_c4$OTU, levels = cause_spp2)
c2 <- ps_c3 %>%
  ggplot(aes(x = date, y = total_abundance + 0.5, color = treatment)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.8) +
  geom_line(data = ps_c4, aes(x = date, y = mean_abundance+0.5, color = treatment)) +
  facet_wrap(~ OTU, scales = "free_y", ncol = 4) +
  scale_y_log10(label= macam::label_10_to_power) +
  labs(x = NULL, y = "Abundance (eDNA copies/ml water + 0.5)",
       title = "Causal taxa to SPAD")

# The number of stems
cause_spp3 <- c("PRO_Taxa00195","EUK_Taxa00423")
ps_cause3 <- prune_taxa(cause_spp3, ps_sample) %>% prune_taxa(taxa_sums(.) > 0, .)
ps_cause3_melt <- psmelt(ps_cause3)
ps_c5 <- ps_cause3_melt %>%
  group_by(date, Week, treatment, plot, OTU) %>% 
  summarize(total_abundance = sum(Abundance))
ps_c6 <- ps_c5 %>%
  group_by(date, Week, treatment, OTU) %>% 
  summarize(mean_abundance = mean(total_abundance))
ps_c5$OTU <- factor(ps_c5$OTU, levels = cause_spp3)
ps_c6$OTU <- factor(ps_c6$OTU, levels = cause_spp3)
c3 <- ps_c5 %>%
  ggplot(aes(x = date, y = total_abundance + 0.5, color = treatment)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.8) +
  geom_line(data = ps_c6, aes(x = date, y = mean_abundance+0.5, color = treatment)) +
  facet_wrap(~ OTU, scales = "free_y", ncol = 4) +
  scale_y_log10(label= macam::label_10_to_power) +
  labs(x = NULL, y = "Abundance (eDNA copies/ml water + 0.5)",
       title = "Causal taxa to the number of stems")


# ------------------------------------- #
# Combine all the figures
# ------------------------------------- #
g_all <- (g1 / g2) + plot_layout(heights = c(1.1, 1))
ggsave("01_VisualizeEcolOut/Fig_eDNAOverallPattern.pdf", g_all, height = 8, width = 10)
ggsave("01_VisualizeEcolOut/Fig_eDNAOverallPattern.jpg", g_all, height = 8, width = 10)
ggsave("01_VisualizeEcolOut/Fig_eDNATopPattern.pdf", t1, height = 8, width = 14)
saveRDS(list(g1, g2), "../../07_FormatFigs/data_robj/eDNA_temporal-pattern.obj")
saveRDS(t1, "../../07_FormatFigs/data_robj/eDNA_top_temporal-pattern.obj")
saveRDS(list(a1, a2, a3), "../../07_FormatFigs/data_robj/eDNA_major_groups.obj")
saveRDS(list(c1, c2, c3), "../../07_FormatFigs/data_robj/eDNA_caual_taxa.obj")
write.csv(ps_top_df, "01_VisualizeEcolOut/top_taxa.csv")
write.csv(ps_top_nobac_df, "01_VisualizeEcolOut/top_taxa_nobac.csv")

# ------------------------------------- #
# Save workspace
# ------------------------------------- #
macam::save_session_info("00_SessionInfo")
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))


