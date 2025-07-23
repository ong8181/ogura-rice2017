####
#### Ogura Rice 2017
#### Visualize MDR S-map results
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2024.08.03
library(phyloseq); packageVersion("phyloseq") # 1.46.0, 2024.11.12
library(ggforce); packageVersion("ggforce") # 0.4.2, 2024.12.09
library(patchwork); packageVersion("patchwork") # 1.3.0, 2024.12.04

# Set random seeds (for reproduction)
set.seed(1234)

# Create output directory
(outdir <- macam::outdir_create())
dir.create("00_SessionInfo")


# ------------------------------------------------ #
# Load all data
# ------------------------------------------------ #
# All data
d_all <- readRDS("01_LoadAllDataOut/d_all.obj") %>% data.frame
ps_ecol <- readRDS("01_LoadAllDataOut/ps_ecol.obj")
d_ecol <- otu_table(ps_ecol) %>% data.frame
d_ecol_tax <- tax_table(ps_ecol) %>% data.frame

# Compile sample names
d_all[,"Sample_Name2"][is.na(d_all$Sample_Name2)] <- sprintf("NoDNA%02s", 1:sum(is.na(d_all$Sample_Name2)))
rownames(d_all) <- d_all$Sample_Name2

# Load compiled MDR S-map results
mdr_hgt <- readRDS("07_CompileMDRresOut/mdr_hgt.obj")
mdr_hgt_long <- readRDS("07_CompileMDRresOut/mdr_hgt_long.obj")

# Extract causal variables
## Height
cause_vars1 <- mdr_hgt %>% select(starts_with("smapc_")) %>% colnames %>% str_split(pattern = "smapc_from_") %>% sapply(`[`, 2)
cause_vars1 <- cause_vars1[!str_detect(cause_vars1, "mean_") & !str_detect(cause_vars1, "clim_")]
cause_vars2 <- cause_vars1 %>% str_split(pattern = "_tp") %>% sapply(`[`, 1)
cause_vars2 <- cause_vars2[!str_detect(cause_vars2, "mean_") & !str_detect(cause_vars2, "clim_")]

# Specify data set
mdr_res_used <- mdr_hgt
mdr_res_long_used <- mdr_hgt_long
cause_var_used1 <- cause_vars1
cause_var_used2 <- cause_vars2


# ------------------------------------------------ #
# Visualize data
# ------------------------------------------------ #
# Extract taxa information from ps object
tax_id_used <- match(cause_var_used2, taxa_names(ps_ecol))
tax_name_df <- tax_table(ps_ecol)[tax_id_used,c("superkingdom","phylum","class","order","family","genus","species")]
tax_name_comb <- tax_name_df %>% apply(1, function(x) str_flatten(x, collapse = "_"))
## Detect the lowest taxa name
tax_na_none <- (is.na(tax_name_df) | (tax_name_df == "")) 
right_most_ture <- apply((!tax_na_none), 1, function(x) if (any(x)) max(which(x)) else NA)
tax_name_comb_lowest <- c(rep(NA, length(right_most_ture)))
for (i in 1:length(right_most_ture)) tax_name_comb_lowest[i] <- tax_name_df[i,right_most_ture[i]]
## Finalize taxa name data.frame
mdr_res_tax <- mdr_res_long_used %>% filter(mdr_res_long_used$cause_var %in% cause_var_used1)
mdr_res_tax$tax_info <- NA
mdr_res_tax$tax_info2 <- NA
for (i in 1:length(cause_var_used1)) {
  tax_i <- cause_var_used1[i]
  if (!is.na(tax_i)){
    mdr_res_tax[mdr_res_tax$cause_var == tax_i, "tax_info"] <- tax_name_comb[i]
    mdr_res_tax[mdr_res_tax$cause_var == tax_i, "tax_info2"] <- tax_name_comb_lowest[i]
  }
}
mdr_res_tax$tax_info2[is.na(mdr_res_tax$tax_info)] <- "Grass_Poa_1"
mdr_res_tax$tax_info[is.na(mdr_res_tax$tax_info)] <- "Grass_Poa_1"

# Summarize tax information
mdr_tax_order <- mdr_res_tax %>%
  group_by(cause_var) %>%
  summarize(median = median(smapc, na.rm = TRUE),
            tax_info = unique(tax_info),
            tax_info2 = unique(tax_info2)) %>% arrange(median) %>% 
  mutate(tax_info3 = paste0(tax_info2, " (", cause_var, ")"))
mdr_res_tax$cause_var <- factor(mdr_res_tax$cause_var, levels = mdr_tax_order$cause_var)
mdr_res_tax$treatment[mdr_res_tax$treatment == "Conv"] <- "Conventional"
mdr_res_tax$treatment[mdr_res_tax$treatment == "NoFert"] <- "No Fertilizer"
mdr_res_tax$treatment <- factor(mdr_res_tax$treatment, levels = c("No Fertilizer", "Conventional"))
# Add sequence information
mdr_tax_order$seq <- mdr_tax_order$cause_var %>% str_split(pattern = "_tp") %>%
  sapply(`[`, 1) %>% match(rownames(d_ecol_tax)) %>% d_ecol_tax[.,"seq"]


# Visualize
## S-map coefficients
g1 <- mdr_res_tax %>% filter(!is.na(smapc)) %>% 
  ggplot(aes(y = cause_var, x = smapc, color = treatment, shape = treatment)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_boxplot(outliers = FALSE) +
  geom_sina(alpha = 0.5, jitter_y = 10, seed = 1234, size = 2) +
  scale_color_manual(values = c("Conventional" = "darkgreen", "No Fertilizer" = "darkgoldenrod2"),
                     breaks = c("Conventional", "No Fertilizer")) +
  scale_shape_manual(values = c(16, 17), breaks = c("Conventional", "No Fertilizer")) +
  scale_y_discrete(labels = mdr_tax_order$tax_info3) + 
  xlab("Influence from ecological variables") + ylab(NULL) +
  ggtitle("Influence on overall rice\nheight growth per unit variable") +
  theme_bw() +
  NULL

## Raw values (Standardized abundance)
g2 <- mdr_res_tax %>% filter(!is.na(raw_value)) %>% 
  ggplot(aes(y = cause_var, x = raw_value, color = treatment, shape = treatment)) +
  geom_boxplot(outliers = FALSE) +
  geom_sina(alpha = 0.5, jitter_y = 10, seed = 1234, size = 2) +
  scale_color_manual(values = c("Conventional" = "darkgreen", "No Fertilizer" = "darkgoldenrod2"),
                     breaks = c("Conventional", "No Fertilizer")) +
  scale_shape_manual(values = c(16, 17), breaks = c("Conventional", "No Fertilizer")) +
  scale_y_discrete(labels = mdr_tax_order$tax_info3) + 
  xlab("Standardized abundance") + ylab(NULL) +
  ggtitle("Abundance of ecological variables\nthat affect rice height growth") +
  theme_bw() +
  NULL

## S-map coefficients x Raw values (Overall effects)
g3 <- mdr_res_tax %>% filter(!is.na(raw_value) & !is.na(smapc)) %>% 
  ggplot(aes(y = cause_var, x = raw_value*smapc, color = treatment, shape = treatment)) +
  geom_boxplot(outliers = FALSE) +
  geom_sina(alpha = 0.5, jitter_y = 10, seed = 1234, size = 2) +
  scale_color_manual(values = c("Conventional" = "darkgreen", "No Fertilizer" = "darkgoldenrod2"),
                     breaks = c("Conventional", "No Fertilizer")) +
  scale_shape_manual(values = c(16, 17), breaks = c("Conventional", "No Fertilizer")) +
  scale_y_discrete(labels = mdr_tax_order$tax_info3) + 
  xlab("Overall effects on rice height growth") + ylab(NULL) +
  ggtitle("Overall effect of ecological\nvariables on rice height growth") +
  theme_bw() + 
  NULL
#g3


# ------------------------------------------------ #
# Save figures
# ------------------------------------------------ #
g_all <- (g1 + theme(legend.position = "none")) +
  (g2 + theme(axis.text.y = element_blank(), legend.position = "none")) +
  (g3 + theme(axis.text.y = element_blank()))
ggsave(sprintf("%s/HeightDiff_overall.pdf", outdir), g_all, width = 14, height = 9)
saveRDS(g_all, "../07_FormatFigs/data_robj/Smap_Height-all.obj")
saveRDS(list(g1, g2, g3), "../07_FormatFigs/data_robj/Smap_Height.obj")


# ------------------------------------------------ #
# Save compiled data
# ------------------------------------------------ #
# Save workspace and session information
save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))
macam::save_session_info()

