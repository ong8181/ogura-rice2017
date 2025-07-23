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
# psmelt_rel, Bacteria
ps_m1_rel <- ps_sample %>% subset_taxa(superkingdom == "Bacteria") %>% 
  transform_sample_counts(function(x) x/sum(x)) %>% 
  speedyseq::psmelt() %>% 
  select(OTU, Sample, Abundance, date, Time, Week, plot, treatment, sample_nc, phylum) %>% 
  filter(Abundance > 0) %>% 
  group_by(date, treatment, phylum) %>% 
  summarize(Abundance = sum(Abundance),
            n_plot = n_distinct(plot))
ps_m1_rel_plot <- ps_m1_rel %>% group_by(date, treatment) %>% summarize(max_n_plot = max(n_plot))
ps_m1_rel_plot$event_id <- paste0(ps_m1_rel_plot$date, "--", ps_m1_rel_plot$treatment)
ps_m1_rel$event_id <- paste0(ps_m1_rel$date, "--", ps_m1_rel$treatment)
ps_m1_rel$max_n_plot <- ps_m1_rel_plot$max_n_plot[match(ps_m1_rel$event_id, ps_m1_rel_plot$event_id)]
ps_m1_rel$relative_abundance <- ps_m1_rel$Abundance/ps_m1_rel$max_n_plot
ps_m1_rel$phylum[ps_m1_rel$phylum == ""] <- "Undetermined"
ps_m1_rel$phylum <- factor(ps_m1_rel$phylum, levels = c(unique(ps_m1_rel$phylum)[2:11], "Undetermined"))

# Relative abundance
g1 <- ps_m1_rel %>% ggplot(aes(x = date, y = relative_abundance, fill = phylum)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 0) +
  scale_fill_manual(values = c(c4a(palette = "cols4all.friendly13")[1:10], "gray80")) +
  facet_wrap(~treatment) +
  xlab(NULL) + ylab("Relative abundance") +
  theme_bw() + ggtitle("Relative abundance of bacteria") +
  NULL

# psmelt_rel, Fungi
ps_m2_rel <- ps_sample %>% subset_taxa(kingdom == "Fungi") %>% 
  transform_sample_counts(function(x) x/sum(x)) %>% 
  speedyseq::psmelt() %>% 
  select(OTU, Sample, Abundance, date, Time, Week, plot, treatment, sample_nc, phylum) %>% 
  filter(Abundance > 0) %>% 
  group_by(date, treatment, phylum) %>% 
  summarize(Abundance = sum(Abundance),
            n_plot = n_distinct(plot))
ps_m2_rel_plot <- ps_m2_rel %>% group_by(date, treatment) %>% summarize(max_n_plot = max(n_plot))
ps_m2_rel_plot$event_id <- paste0(ps_m2_rel_plot$date, "--", ps_m2_rel_plot$treatment)
ps_m2_rel$event_id <- paste0(ps_m2_rel$date, "--", ps_m2_rel$treatment)
ps_m2_rel$max_n_plot <- ps_m2_rel_plot$max_n_plot[match(ps_m2_rel$event_id, ps_m2_rel_plot$event_id)]
ps_m2_rel$relative_abundance <- ps_m2_rel$Abundance/ps_m2_rel$max_n_plot
ps_m2_rel$phylum[ps_m2_rel$phylum == ""] <- "Undetermined"
ps_m2_rel$phylum <- factor(ps_m2_rel$phylum, levels = c(unique(ps_m2_rel$phylum)[2:6], "Undetermined"))

# Relative abundance
g2 <- ps_m2_rel %>% ggplot(aes(x = date, y = relative_abundance, fill = phylum)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 0) +
  scale_fill_manual(values = c(c4a(palette = "cols4all.friendly13")[1:5], "gray80")) +
  facet_wrap(~treatment) +
  xlab(NULL) + ylab("Relative abundance") +
  theme_bw() + ggtitle("Relative abundance of fungi") +
  NULL

# psmelt_rel, non-fungal eukaryota
ps_m3_rel <- ps_sample %>% subset_taxa(superkingdom == "Eukaryota" & kingdom != "Fungi") %>% 
  transform_sample_counts(function(x) x/sum(x)) %>% 
  speedyseq::psmelt() %>% 
  select(OTU, Sample, Abundance, date, Time, Week, plot, treatment, sample_nc, phylum) %>% 
  filter(Abundance > 0) %>% 
  group_by(date, treatment, phylum) %>% 
  summarize(Abundance = sum(Abundance),
            n_plot = n_distinct(plot))
ps_m3_rel_plot <- ps_m3_rel %>% group_by(date, treatment) %>% summarize(max_n_plot = max(n_plot))
ps_m3_rel_plot$event_id <- paste0(ps_m3_rel_plot$date, "--", ps_m3_rel_plot$treatment)
ps_m3_rel$event_id <- paste0(ps_m3_rel$date, "--", ps_m3_rel$treatment)
ps_m3_rel$max_n_plot <- ps_m3_rel_plot$max_n_plot[match(ps_m3_rel$event_id, ps_m3_rel_plot$event_id)]
ps_m3_rel$relative_abundance <- ps_m3_rel$Abundance/ps_m3_rel$max_n_plot
ps_m3_rel$phylum[ps_m3_rel$phylum == ""] <- "Undetermined"
ps_m3_rel$phylum <- factor(ps_m3_rel$phylum, levels = c(unique(ps_m3_rel$phylum)[2:21], "Undetermined"))

# Relative abundance
c4a.friendly2 <- c("#C73F74","#D1CE63","#17AADD","#209988","#3467BD","#DF9D9A","#69DDFF","gray80")
g3 <- ps_m3_rel %>% ggplot(aes(x = date, y = relative_abundance, fill = phylum)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 0) +
  scale_fill_manual(values = c(c4a(palette = "cols4all.friendly13"), c4a.friendly2)) +
  facet_wrap(~treatment) +
  xlab(NULL) + ylab("Relative abundance") +
  theme_bw() + ggtitle("Relative abundance of non-fungi eukaryota") +
  NULL


# ------------------------------------- #
# Combine all the figures
# ------------------------------------- #
saveRDS(list(g1, g2, g3), "../../07_FormatFigs/data_robj/eDNA_relative.obj")

# ------------------------------------- #
# Save workspace
# ------------------------------------- #
macam::save_session_info("00_SessionInfo")
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))


