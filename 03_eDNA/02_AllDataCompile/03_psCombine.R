####
#### OGR eDNA study 2017
#### Combine all phyloseq objects
####

# Load library and functions
library(tidyverse); packageVersion("tidyverse") #2.0.0, 2024.08.01
library(cowplot); packageVersion("cowplot") #1.1.3, 2024.08.01
library(phyloseq); packageVersion("phyloseq") #1.46.0, 2024.08.02
library(ggsci); packageVersion("ggsci") #3.0.0, 2024.08.02
theme_set(theme_cowplot())

# Set random seeds (for reproduction)
ran.seed <- 8181
set.seed(ran.seed)

# Create output directory
wdir <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(wdir, end = -3), "Out")); rm(wdir)
dir.create(output_folder)


# ------------------------------------------------------------ #
# Load data
# ------------------------------------------------------------ #
ps_pro_sample <- readRDS("02_DNAcopiesOut/ps_pro_sample.obj")
ps_fun_sample <- readRDS("02_DNAcopiesOut/ps_fun_sample.obj")
ps_inv_sample <- readRDS("02_DNAcopiesOut/ps_inv_sample.obj")
ps_euk_sample <- readRDS("02_DNAcopiesOut/ps_euk_sample.obj")

# Combine all phyloseq object
sample_meta_all <- as.data.frame(sample_data(ps_pro_sample)[,c(4:5,13:26)])
sample_meta_pro <- as.data.frame(sample_data(ps_pro_sample)[,28:35])
sample_meta_fun <- as.data.frame(sample_data(ps_fun_sample)[,28:35])
sample_meta_inv <- as.data.frame(sample_data(ps_inv_sample)[,28:35])
sample_meta_euk <- as.data.frame(sample_data(ps_euk_sample)[,28:35])

meta_colnames <- colnames(as.data.frame(sample_data(ps_pro_sample)))[28:35]
colnames(sample_meta_pro) <- sprintf("PRO_%s", meta_colnames)
colnames(sample_meta_fun) <- sprintf("ITS_%s", meta_colnames)
colnames(sample_meta_inv) <- sprintf("COI_%s", meta_colnames)
colnames(sample_meta_euk) <- sprintf("EUK_%s", meta_colnames)

sample_combined <- cbind(sample_meta_all, sample_meta_pro, sample_meta_fun, sample_meta_inv, sample_meta_euk)

# Pre-combine
ps_combined0 <- merge_phyloseq(ps_pro_sample, ps_fun_sample, ps_inv_sample, ps_euk_sample)

# Re-order taxtable information
potential_taxcol_names <- c("query",
                            "superkingdom", "kingdom", "subkingdom",
                            "phylum",
                            "class", "subclass",
                            "superorder", "order", "suborder", "family",
                            "tribe",
                            "genus", "species",
                            "seq", "seqlen", "miseq_run")
tax_table(ps_combined0) <- tax_table(ps_combined0)[,potential_taxcol_names]

# Re-merge phyloseq object (this process is done to keep sample information)
ps_combined <- phyloseq(otu_table(ps_combined0),
                        sample_data(sample_combined), # Keep all sample information
                        tax_table(ps_combined0))

# Save and output results
saveRDS(ps_combined, sprintf("%s/ps_combined.obj", output_folder))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/03_SessionInfo_psCombine_%s.txt", substr(Sys.time(), 1, 10)))
