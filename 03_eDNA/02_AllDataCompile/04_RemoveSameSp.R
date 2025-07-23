####
#### OGR eDNA study 2017
#### Remove OTUs originated from the same species
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
ps_combined <- readRDS("03_psCombineOut/ps_combined.obj")

# Calculate correlation coefficients
ps_cor <- cor(as.data.frame(otu_table(ps_combined)@.Data))
diag(ps_cor) <- 0
ps_cor[upper.tri(ps_cor)] <- 0

# Find highly correlated ASV columns
cor_set <- data.frame(NULL)
for(i in 1:nrow(ps_cor)){
  cor_set_j <- data.frame(NULL)
  j0 <- which(ps_cor[,i] > 0.98)
  if(length(j0) > 0 & length(j0) < 2){
    cor_set_j <- data.frame(taxa1 = i, taxa2 = j0, cor = ps_cor[j0,i])
  }else if(length(j0) > 1){
    for(j1 in 1:length(j0)){
      cor_set_j <- rbind(cor_set_j, data.frame(taxa1 = i, taxa2 = j0[j1], cor = ps_cor[j0[j1],i]))
    }
  }
  cor_set <- rbind(cor_set, cor_set_j)
}

# Add taxa names
cor_set_taxa1 <- tax_table(ps_combined)[cor_set$taxa1, c("superkingdom", "phylum", "family", "genus", "miseq_run")]@.Data
cor_set_taxa2 <- tax_table(ps_combined)[cor_set$taxa2, c("superkingdom", "phylum", "family", "genus", "miseq_run")]@.Data
cor_set$taxa_name1 <- sprintf("%s_%s_%s_%s_%s", cor_set_taxa1[,5], cor_set_taxa1[,1], cor_set_taxa1[,2], cor_set_taxa1[,3], cor_set_taxa1[,4])
cor_set$taxa_name2 <- sprintf("%s_%s_%s_%s_%s", cor_set_taxa2[,5], cor_set_taxa2[,1], cor_set_taxa2[,2], cor_set_taxa2[,3], cor_set_taxa2[,4])
cor_set$taxa_id1 <- taxa_names(ps_combined)[cor_set$taxa1]
cor_set$taxa_id2 <- taxa_names(ps_combined)[cor_set$taxa2]
rownames(cor_set) <- 1:nrow(cor_set)

# Compare taxa information
cor_set$assign_same_taxa <- NaN
# Check NA
for(cor_tax_i in 1:nrow(cor_set)){
  comp_id <- intersect(which(cor_set_taxa1[cor_tax_i,1:4] != ""), which(cor_set_taxa2[cor_tax_i,1:4] != ""))
  cor_set$assign_same_taxa[cor_tax_i] <- all(cor_set_taxa1[cor_tax_i,comp_id] == cor_set_taxa2[cor_tax_i,comp_id])
}

# Output CSV file
write.csv(cor_set, sprintf("%s/cor_set.csv", output_folder), row.names = F)

# Check correlation figures
pdf(sprintf("%s/cor_all.pdf", output_folder), width = 45, height = 45)
op <- par(mfrow=c(23,23))
for(i in 1:nrow(cor_set)){
    plot(otu_table(ps_combined)[,cor_set[i,1]], otu_table(ps_combined)[,cor_set[i,2]],
       xlab = cor_set[i,4], ylab = cor_set[i,5], main = cor_set[i,3])
  abline(0,1)
}
par(op)
dev.off()


# Identify taxa names with more than r = 0.98
sum(cor_set$cor > 0.98)
cor_set_filter <- cor_set[cor_set$cor > 0.98 & cor_set$assign_same_taxa == 1,]
cor_set_filter$taxa1_sum <- taxa_sums(ps_combined)[cor_set_filter$taxa1]
cor_set_filter$taxa2_sum <- taxa_sums(ps_combined)[cor_set_filter$taxa2]
cor_set_filter$keep <- NaN
cor_set_filter$remove <- NaN
keep_taxa1 <- cor_set_filter$taxa1_sum - cor_set_filter$taxa2_sum >= 0
keep_taxa2 <- cor_set_filter$taxa1_sum - cor_set_filter$taxa2_sum < 0
cor_set_filter$keep[keep_taxa1] <- taxa_names(ps_combined)[cor_set_filter$taxa1[keep_taxa1]]
cor_set_filter$keep[keep_taxa2] <- taxa_names(ps_combined)[cor_set_filter$taxa2[keep_taxa2]]
cor_set_filter$remove[!keep_taxa1] <- taxa_names(ps_combined)[cor_set_filter$taxa1[!keep_taxa1]]
cor_set_filter$remove[!keep_taxa2] <- taxa_names(ps_combined)[cor_set_filter$taxa2[!keep_taxa2]]
remove_cor_taxa <- unique(cor_set_filter$remove)
filtered_taxa <- taxa_names(ps_combined)[!(taxa_names(ps_combined) %in% remove_cor_taxa)]
ps_filt <- prune_taxa(filtered_taxa, ps_combined)
ps_filt

# Save phyloseq objects
saveRDS(ps_filt, sprintf("%s/ps_comb_filt.obj", output_folder))
write.csv(sample_data(ps_filt) %>% data.frame, sprintf("%s/sample_data.csv", output_folder))
write.csv(otu_table(ps_filt) %>% data.frame, sprintf("%s/otu_table.csv", output_folder))
write.csv(tax_table(ps_filt) %>% data.frame, sprintf("%s/tax_table.csv", output_folder))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/04_SessionInfo_RemoveSampeSp_%s.txt", substr(Sys.time(), 1, 10)))
