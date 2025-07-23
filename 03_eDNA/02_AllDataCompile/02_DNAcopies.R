####
#### OGR eDNA study 2017
#### Visualize DNA concentration
####

# Load library and functions
library(tidyverse); packageVersion("tidyverse") #2.0.0, 2024.08.01
library(cowplot); packageVersion("cowplot") #1.1.3, 2024.08.01
library(reshape2); packageVersion("reshape2") #1.4.4, 2024.08.01
library(phyloseq); packageVersion("phyloseq") #1.46.0, 2024.08.02
library(ggsci); packageVersion("ggsci") #3.0.0, 2024.08.02
theme_set(theme_cowplot())

# Load workspace
load("01_PhyloseqImportOut/01_PhyloseqImportOut.RData")

# Set random seeds (for reproduction)
ran.seed <- 8181
set.seed(ran.seed)

# Create output directory
wdir <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(wdir, end = -3), "Out")); rm(wdir)
dir.create(output_folder)


# ------------------------------------------------------------ #
# Extracting data as time series
# ------------------------------------------------------------ #
ps_pro_sample2 <- prune_taxa(taxa_sums(ps_pro_sample) > 0, ps_pro_sample)
ps_fun_sample2 <- prune_taxa(taxa_sums(ps_fun_sample) > 0, ps_fun_sample)
ps_inv_sample2 <- prune_taxa(taxa_sums(ps_inv_sample) > 0, ps_inv_sample)
ps_euk_sample2 <- prune_taxa(taxa_sums(ps_euk_sample) > 0, ps_euk_sample)

# Check the number of taxa detected
(tax_n_pro <- ncol(otu_table(ps_pro_sample2))) # 7053
(tax_n_fun <- ncol(otu_table(ps_fun_sample2))) # 3379
(tax_n_inv <- ncol(otu_table(ps_inv_sample2))) # 5199
(tax_n_euk <- ncol(otu_table(ps_euk_sample2))) # 5404
sum(tax_n_pro, tax_n_fun, tax_n_inv, tax_n_euk)
# 21035 taxa in total


# ------------------------------------------------------------ #
# Visualize copy numbers
# ------------------------------------------------------------ #
# Preparations
pro_conc_df <- data.frame(conc = c(otu_table(ps_pro_sample2)))
fun_conc_df <- data.frame(conc = c(otu_table(ps_fun_sample2)))
inv_conc_df <- data.frame(conc = c(otu_table(ps_inv_sample2)))
euk_conc_df <- data.frame(conc = c(otu_table(ps_euk_sample2)))

# Histogram
g8_1 <- ggplot(pro_conc_df, aes(x = conc)) + geom_histogram() + scale_x_log10() +
  xlab("eDNA (copies/ul)") + ggtitle("16S eDNA copy numbers")
g8_2 <- ggplot(fun_conc_df, aes(x = conc)) + geom_histogram() + scale_x_log10() +
  xlab("eDNA (copies/ul)") + ggtitle("ITS eDNA copy numbers")
g8_3 <- ggplot(inv_conc_df, aes(x = conc)) + geom_histogram() + scale_x_log10() +
  xlab("eDNA (copies/ul)") + ggtitle("COI eDNA copy numbers")
g8_4 <- ggplot(euk_conc_df, aes(x = conc)) + geom_histogram() + scale_x_log10() +
  xlab("eDNA (copies/ul)") + ggtitle("18S eDNA copy numbers")
g8_all <- plot_grid(g8_1, g8_2, g8_3, g8_4, ncol = 2, labels = "auto", align = "hv")

ggsave(sprintf("%s/DNAcopy.pdf", output_folder), plot = g8_all, width = 10, height = 8)


# ------------------------------------------------------------ #
# Remove low-frequency and low-concentration DNAs
# ------------------------------------------------------------ #
# Set the lowest DNA copy numbers of standard DNAs
ll_pro <- 2500
ll_fun <- 20
ll_inv <- 3.125
ll_euk <- 25
ll_lim <- 2
freq_lim <- 0.1

# Calculate mean DNA copy numbers for each taxon (less than a half of ll)
low_tax_pro <- taxa_sums(ps_pro_sample2)/nrow(sample_data(ps_pro_sample2)) < ll_pro/ll_lim
low_tax_fun <- taxa_sums(ps_fun_sample2)/nrow(sample_data(ps_fun_sample2)) < ll_fun/ll_lim
low_tax_inv <- taxa_sums(ps_inv_sample2)/nrow(sample_data(ps_inv_sample2)) < ll_inv/ll_lim
low_tax_euk <- taxa_sums(ps_euk_sample2)/nrow(sample_data(ps_euk_sample2)) < ll_euk/ll_lim

# Calculate detection frequency for each taxon (less than 10%)
low_freq_pro <- apply(otu_table(ps_pro_sample2) > 0, 2, sum)/nrow(sample_data(ps_pro_sample2)) < freq_lim
low_freq_fun <- apply(otu_table(ps_fun_sample2) > 0, 2, sum)/nrow(sample_data(ps_fun_sample2)) < freq_lim
low_freq_inv <- apply(otu_table(ps_inv_sample2) > 0, 2, sum)/nrow(sample_data(ps_inv_sample2)) < freq_lim
low_freq_euk <- apply(otu_table(ps_euk_sample2) > 0, 2, sum)/nrow(sample_data(ps_euk_sample2)) < freq_lim

ps_pro_sample3 <- prune_taxa(!(low_tax_pro | low_freq_pro), ps_pro_sample2)
ps_fun_sample3 <- prune_taxa(!(low_tax_fun | low_freq_fun), ps_fun_sample2)
ps_inv_sample3 <- prune_taxa(!(low_tax_inv | low_freq_inv), ps_inv_sample2)
ps_euk_sample3 <- prune_taxa(!(low_tax_euk | low_freq_euk), ps_euk_sample2)

ntaxa(ps_pro_sample3) # 210
ntaxa(ps_fun_sample3) # 226
ntaxa(ps_inv_sample3) # 309
ntaxa(ps_euk_sample3) # 620
# 210 + 226 + 309 + 620 = 1365


# ------------------------------------------------------------ #
# Save results
# ------------------------------------------------------------ #
# Save phyloseq objects
saveRDS(ps_pro_sample3, sprintf("%s/ps_pro_sample.obj", output_folder))
saveRDS(ps_fun_sample3, sprintf("%s/ps_fun_sample.obj", output_folder))
saveRDS(ps_inv_sample3, sprintf("%s/ps_inv_sample.obj", output_folder))
saveRDS(ps_euk_sample3, sprintf("%s/ps_euk_sample.obj", output_folder))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/02_SessionInfo_DNACopies_%s.txt", substr(Sys.time(), 1, 10)))
