####
#### OGR eDNA study 2017
#### Sequence reads visualization
####

# Load library
library(tidyverse); packageVersion("tidyverse") #2.0.0, 2024.08.01
library(cowplot); packageVersion("cowplot") #1.1.3, 2024.08.01
library(reshape2); packageVersion("reshape2") #1.4.4, 2024.08.01
library(scales); packageVersion("scales") #1.3.0, 2024.08.01
theme_set(theme_cowplot())

# Set random seeds (for reproduction)
ran.seed <- 8181
set.seed(ran.seed)

# Create output directory
wdir <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(wdir, end = -3), "Out")); rm(wdir)
dir.create(output_folder)


# ------------------------------------- #
# Load data
# ------------------------------------- #
track_df <- read.csv("02_DADA2Out/track.csv", row.names = 1)
seqtab_nochim <- readRDS("07_STDCheckOut/seqtab_nochim.obj")
seqtab_conv <- read.csv("07_STDCheckOut/seqtab_conv.csv", row.names = 1)
tax_claident2 <- read.csv("07_STDCheckOut/tax_w_std.csv", row.names = 1)
taxa_wo_std <- read.csv("07_STDCheckOut/tax_wo_std.csv", row.names = 1)
sample_sheet <- read.csv("07_STDCheckOut/sample_sheet.csv", row.names = 1)
n_std_seq <- readRDS("07_STDCheckOut/n_std_seq.obj")
new_std_table <- readRDS("07_STDCheckOut/new_std_table.obj")
coef_summary <- readRDS("07_STDCheckOut/coef_summary.obj")
r2_summary <- readRDS("07_STDCheckOut/r2_summary.obj")


# ------------------------------------- #
# Visualize data
# ------------------------------------- #
# Track data
track_df$sample <- rownames(track_df)
track_df2 <- melt(track_df[,c(1:7,9)], id.vars = c("sample"))
v1 <- ggplot(track_df2, aes(x = variable, y = sample, fill = value)) +
  geom_tile(colour = "black", size = 0.001) + scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", midpoint = 10000, limits = c(0, 30000), oob = squish) +
  theme_gray(base_size = 8) + theme(axis.text.y = element_text(size = 1)) +
  geom_hline(yintercept = which(sample_sheet$sample_nc == "pcr_nc"), size = 0.1, alpha = 0.5)
ggsave(plot = v1, filename = sprintf("%s/Track.pdf", output_folder), width = 10, height = 10)

# Sequence table
n_taxa <- 1:200
seqtab01 <- data.frame(seqtab_nochim)[,n_taxa]
colnames(seqtab01) <- sprintf("%s_%s_%s", tax_claident2$query[n_taxa], tax_claident2$phylum[n_taxa], tax_claident2$species[n_taxa])
seqtab01$sample <- rownames(seqtab01)
seqtab01_df <- melt(seqtab01, id.vars = c("sample"))
v2 <- ggplot(seqtab01_df, aes(x = variable, y = sample, fill = value)) +
  geom_tile() + scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", midpoint = 2000, oob = squish) +
  theme_gray(base_size = 8) + theme(axis.text.y = element_text(size = 1),
                                    axis.text.x = element_text(angle = 90, size = 3)) +
  geom_hline(yintercept = which(sample_sheet$sample_nc == "pcr_nc"), size = 0.1, alpha = 0.5)
ggsave(plot = v2, filename = sprintf("%s/SeqtabNochim.pdf", output_folder), width = 10, height = 10)

# Sequence table without STD
seqtab02 <- data.frame(seqtab_conv)[,n_taxa]
colnames(seqtab02) <- sprintf("%s_%s_%s", tax_claident2$query[-n_std_seq][n_taxa], tax_claident2$phylum[-n_std_seq][n_taxa], tax_claident2$species[-n_std_seq][n_taxa])
seqtab02$sample <- rownames(seqtab02)
seqtab02_df <- melt(seqtab02, id.vars = c("sample"))
v3 <- ggplot(seqtab02_df, aes(x = variable, y = sample, fill = value)) +
  geom_tile() + scale_fill_gradient2(low = "royalblue", mid = "white", high = "red3", midpoint = 200000, oob = squish) +
  theme_gray(base_size = 8) + theme(axis.text.y = element_text(size = 1),
                                    axis.text.x = element_text(angle = 90, size = 3)) +
  geom_hline(yintercept = which(sample_sheet$sample_nc == "pcr_nc"), size = 0.1, alpha = 0.5)
ggsave(plot = v3, filename = sprintf("%s/SeqtabConvNoSTD.pdf", output_folder), width = 10, height = 10)


# ------------------------------------------ #
# Summary Table
# ------------------------------------------ #
sample_summary <- sample_sheet
track2 <- data.frame(matrix(NA, nrow = nrow(sample_summary), ncol = ncol(track_df)))
rownames(track2) <- sample_summary$Sample_Name2
colnames(track2) <- colnames(track_df)
track2[match(rownames(track_df), rownames(track2)),] <- track_df
sample_summary$dada2_input <- track2$input
sample_summary$dada2_nochim <- track2$nochim
sample_summary$dada2_prop <- track2$prop
sample_summary$STD_all <- rowSums(new_std_table)
sample_summary$NonSTD_all <- rowSums(seqtab_nochim) - rowSums(new_std_table)
sample_summary$STD_prop <- sample_summary$STD_all / rowSums(seqtab_nochim)
sample_summary$STD_coef <- sample_summary$STD_r2 <- NA
sample_summary$STD_coef[match(names(coef_summary), sample_summary$Sample_Name2)] <- coef_summary
sample_summary$STD_r2[match(names(r2_summary), sample_summary$Sample_Name2)] <- r2_summary
sample_summary$dna_copy_sum <- rowSums(seqtab_conv)

# Repalce sample_sheet
sample_sheet <- sample_summary

# Save important objects
write.csv(sample_summary, sprintf("%s/summary_sheet.csv", output_folder))
saveRDS(sample_sheet, sprintf("%s/sample_sheet.obj", output_folder))
saveRDS(seqtab_nochim, sprintf("%s/seqtab_nochim.obj", output_folder))
saveRDS(seqtab_conv, sprintf("%s/seqtab_conv.obj", output_folder))
saveRDS(new_std_table, sprintf("%s/new_std_table.obj", output_folder))
saveRDS(tax_claident2, sprintf("%s/taxa_w_std.obj", output_folder))
saveRDS(taxa_wo_std, sprintf("%s/taxa_wo_std.obj", output_folder))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/08_SessionInfo_ReadNSummary_%s.txt", substr(Sys.time(), 1, 10)))

