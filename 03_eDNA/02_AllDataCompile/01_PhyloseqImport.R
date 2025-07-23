####
#### OGR Rice eDNA study 2017
#### Importing all MiSeq run data to phyloseq objects
####

# Load library and functions
library(tidyverse); packageVersion("tidyverse") #2.0.0, 2024.08.01
library(cowplot); packageVersion("cowplot") #1.1.3, 2024.08.01
library(reshape2); packageVersion("reshape2") #1.4.4, 2024.08.01
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
dir.create("00_SessionInfo")

# Specify data folders
# RMR-080 Prokaryote
data_folder_16S <- "../01_SequenceProcessing/OGRrice2017_16S/08_ReadNSummaryOut/"
# RMR-080 Fungi
data_folder_ITS <- "../01_SequenceProcessing/OGRrice2017_ITS/08_ReadNSummaryOut/"
# RMR-103 Invertebrate
data_folder_COI <- "../01_SequenceProcessing/OGRrice2017_COI/08_ReadNSummaryOut/"
# CMR-007 Eukaryote
data_folder_18S <- "../01_SequenceProcessing/OGRrice2017_18S/08_ReadNSummaryOut/"


# ----------------------------------------------------- #
# RMR-080 16S
# ----------------------------------------------------- #
# Load objects
seqtab_conv_prok <- readRDS(sprintf("%s/seqtab_conv.obj", data_folder_16S))
sample_sheet_prok <- readRDS(sprintf("%s/sample_sheet.obj", data_folder_16S))
taxa_wo_std_prok <- readRDS(sprintf("%s/taxa_wo_std.obj", data_folder_16S))
# Preparetion to import to phyloseq
dim(sample_sheet_prok); dim(taxa_wo_std_prok); dim(seqtab_conv_prok)
all(rownames(sample_sheet_prok) == rownames(seqtab_conv_prok)) # sample name check
taxa_wo_std_prok$miseq_run <- "RMR-080-16S"
# Change taxa name to group-specific name
taxa_name_prok <- sprintf("PRO_%s", rownames(taxa_wo_std_prok))
colnames(seqtab_conv_prok) <- rownames(taxa_wo_std_prok) <- taxa_name_prok # change col name
all(rownames(taxa_wo_std_prok) == colnames(seqtab_conv_prok)) # taxa name check
# Import data to phyloseq
ps_pro_all <- phyloseq(otu_table(seqtab_conv_prok, taxa_are_rows=FALSE),
                       sample_data(sample_sheet_prok),
                       tax_table(as.matrix(taxa_wo_std_prok)))

# ----------------------------------------------------- #
# RMR-080 ITS
# ----------------------------------------------------- #
# Load objects
seqtab_conv_fungi <- readRDS(sprintf("%s/seqtab_conv.obj", data_folder_ITS))
sample_sheet_fungi <- readRDS(sprintf("%s/sample_sheet.obj", data_folder_ITS))
taxa_wo_std_fungi <- readRDS(sprintf("%s/taxa_wo_std.obj", data_folder_ITS))
# Preparetion to import to phyloseq
dim(sample_sheet_fungi); dim(taxa_wo_std_fungi); dim(seqtab_conv_fungi)
all(rownames(sample_sheet_fungi) == rownames(seqtab_conv_fungi)) # sample name check
taxa_wo_std_fungi$miseq_run <- "RMR-080-ITS"
# Change taxa name to group-specific name
taxa_name_fungi <- sprintf("ITS_%s", rownames(taxa_wo_std_fungi))
colnames(seqtab_conv_fungi) <- rownames(taxa_wo_std_fungi) <- taxa_name_fungi # change col name
all(rownames(taxa_wo_std_fungi) == colnames(seqtab_conv_fungi)) # taxa name check
# Import data to phyloseq
ps_fun_all <- phyloseq(otu_table(seqtab_conv_fungi, taxa_are_rows=FALSE),
                       sample_data(sample_sheet_fungi),
                       tax_table(as.matrix(taxa_wo_std_fungi)))

# ----------------------------------------------------- #
# RMR-103 COI
# ----------------------------------------------------- #
# Load objects
seqtab_conv_inv <- readRDS(sprintf("%s/seqtab_conv.obj", data_folder_COI))
sample_sheet_inv <- readRDS(sprintf("%s/sample_sheet.obj", data_folder_COI))
taxa_wo_std_inv <- readRDS(sprintf("%s/taxa_wo_std.obj", data_folder_COI))
# Preparetion to import to phyloseq
dim(sample_sheet_inv); dim(taxa_wo_std_inv); dim(seqtab_conv_inv)
all(rownames(sample_sheet_inv) == rownames(seqtab_conv_inv)) # sample name check
taxa_wo_std_inv$miseq_run <- "RMR-103-COI"
# Change taxa name to group-specific name
taxa_name_inv <- sprintf("COI_%s", rownames(taxa_wo_std_inv))
colnames(seqtab_conv_inv) <- rownames(taxa_wo_std_inv) <- taxa_name_inv # change col name
all(rownames(taxa_wo_std_inv) == colnames(seqtab_conv_inv)) # taxa name check
# Import data to phyloseq
ps_inv_all <- phyloseq(otu_table(seqtab_conv_inv, taxa_are_rows=FALSE),
                         sample_data(sample_sheet_inv),
                         tax_table(as.matrix(taxa_wo_std_inv)))

# ----------------------------------------------------- #
# CMR-007 18S
# ----------------------------------------------------- #
# Load objects
seqtab_conv_euk <- readRDS(sprintf("%s/seqtab_conv.obj", data_folder_18S))
sample_sheet_euk <- readRDS(sprintf("%s/sample_sheet.obj", data_folder_18S))
taxa_wo_std_euk <- readRDS(sprintf("%s/taxa_wo_std.obj", data_folder_18S))
# Preparetion to import to phyloseq
dim(sample_sheet_euk); dim(taxa_wo_std_euk); dim(seqtab_conv_euk)
all(rownames(sample_sheet_euk) == rownames(seqtab_conv_euk)) # sample name check
taxa_wo_std_euk$miseq_run <- "CMR-007-18S"
# Change taxa name to group-specific name
taxa_name_euk <- sprintf("EUK_%s", rownames(taxa_wo_std_euk))
colnames(seqtab_conv_euk) <- rownames(taxa_wo_std_euk) <- taxa_name_euk # change col name
all(rownames(taxa_wo_std_euk) == colnames(seqtab_conv_euk)) # taxa name check
# Import data to phyloseq
ps_euk_all <- phyloseq(otu_table(seqtab_conv_euk, taxa_are_rows=FALSE),
                       sample_data(sample_sheet_euk),
                       tax_table(as.matrix(taxa_wo_std_euk)))


# ----------------------------------------------------- #
# Extracting "sample" only
# ----------------------------------------------------- #
# Removing field negative/pcr negative/standard negative/positive control samples
ps_pro_sample <- prune_samples(sample_data(ps_pro_all)$sample_nc == "sample", ps_pro_all)
ps_fun_sample <- prune_samples(sample_data(ps_fun_all)$sample_nc == "sample", ps_fun_all)
ps_inv_sample <- prune_samples(sample_data(ps_inv_all)$sample_nc == "sample", ps_inv_all)
ps_euk_sample <- prune_samples(sample_data(ps_euk_all)$sample_nc == "sample", ps_euk_all)


# Save and output results
save(list = ls(all.names = TRUE),
     file = sprintf("%s/01_PhyloseqImportOut.RData", output_folder))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/01_SessionInfo_PhyloseqImport_%s.txt", substr(Sys.time(), 1, 10)))

