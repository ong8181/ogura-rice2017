####
#### Clustering DADA2 results
#### Illumina PE
####

# Load DADA2 results
load("02_DADA2Out/DADA2Out.RData")

# Load library and functions
library(tidyverse); packageVersion("tidyverse") #2.0.0, 2023.8.22
library(phyloseq); packageVersion("phyloseq") #1.44.0, 2024.7.26
library(ShortRead); packageVersion("ShortRead") #1.58.0, 2024.7.26
library(dada2); packageVersion("dada2") #1.28.0, 2024.7.26
library(DECIPHER); packageVersion("DECIPHER") #2.28.0, 2024.7.26
library(Biostrings); packageVersion("Biostrings") # 2.68.1, 2024.7.26
library(speedyseq); packageVersion("speedyseq") # 0.5.3.9018, 2022.7.27

# Set random seeds (for reproduction)
ran.seed <- 8181
set.seed(ran.seed)
wdir <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(wdir, end = -3), "Out")); rm(wdir)
dir.create(output_folder)


# <-----------------------------------------------------> #
#                       Check data
# <-----------------------------------------------------> #
# Check sample data
dim(seqtab_nochim); sum(seqtab_nochim)


# <-----------------------------------------------------> #
#                  Clustering by DECIPHER
# <-----------------------------------------------------> #
dna <- Biostrings::DNAStringSet(seqs) # "seqs" is an object from DADA2 output
clusters <- DECIPHER::Clusterize(dna, cutoff = 0.03, processors = NULL)
colnames(clusters) <- "cluster_id"
clusters$cluster_id <- factor(clusters$cluster_id, levels = unique(clusters$cluster_id))
length(unique(clusters$cluster_id))
## Use dplyr to merge the columns of the seqtab matrix for ASVs in the same OTU
# prep by adding sequences to the `clusters` data frame
merged_seqtab <- seqtab_nochim %>% t %>%
  rowsum(clusters$cluster_id) %>% t
dim(merged_seqtab)

# <-----------------------------------------------------> #
#             Make OTU-based phyloseq objects
# <-----------------------------------------------------> #
# Import DADA2 ASV output to phyloseq
ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = FALSE))
# Merge taxa in "ps" using cluster information
ps_otu <- speedyseq::merge_taxa_vec(ps, group = clusters$cluster_id)
otu_seqs <- colnames(otu_table(ps_otu))
otu_only <- data.frame(taxa_id = sprintf("OTU%05d", 1:length(otu_seqs)),
                       seq = otu_seqs)
write.csv(otu_only, sprintf("%s/otu_only.csv", output_folder), row.names = T)

# Check correspondence between ASV sequences and OTU representative sequences
clusters$seq_sum <- colSums(seqtab_nochim)
clusters$asv_seq <- seqs
asv_1st <- clusters %>% group_by(cluster_id) %>% summarize(otu_seqs = asv_seq[[1]])
clusters$otu_seq <- unlist(asv_1st[match(clusters$cluster_id, asv_1st$cluster_id), "otu_seqs"])
write.csv(clusters, sprintf("%s/cluster_summary.csv", output_folder), row.names = F)

# Save OTU table
otu_out <- as.matrix(c(rbind(sprintf(">OTU%05d", 1:length(otu_seqs)), otu_seqs)), ncol=1)
write.table(otu_out, sprintf("%s/OTU_seqs.fa", output_folder), col.names = FALSE, row.names = FALSE, quote = FALSE)

otu_mat <- as.matrix(otu_table(ps_otu)@.Data)
colnames(otu_mat) <- sprintf("OTU%05d", 1:length(otu_seqs))
write.csv(otu_mat, sprintf("%s/otu_table.csv", output_folder), row.names = T)


# <-----------------------------------------------------> #
#                     Save workspace
# <-----------------------------------------------------> #
# Save workspace
save(list = ls(all.names = TRUE),
     file = paste0(output_folder, "/", output_folder, ".RData"))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))


