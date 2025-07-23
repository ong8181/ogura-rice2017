####
#### OGR eDNA study 2017
#### Negative control check
####

# Load library
library(tidyverse); packageVersion("tidyverse") #2.0.0, 2024.08.01
library(cowplot); packageVersion("cowplot") #1.1.3, 2024.08.01
library(reshape2); packageVersion("reshape2") #1.4.4, 2024.08.01
library(scales); packageVersion("scales") #1.3.0, 2024.08.01
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


# -------------------------------------------- #
# Load data
# -------------------------------------------- #
seqtab_nochim <- readRDS("08_ReadNSummaryOut/seqtab_nochim.obj")
seqtab_conv <- readRDS("08_ReadNSummaryOut/seqtab_conv.obj")
tax_claident2 <- readRDS("08_ReadNSummaryOut/taxa_w_std.obj")
taxa_wo_std <- readRDS("08_ReadNSummaryOut/taxa_wo_std.obj")
sample_sheet <- readRDS("08_ReadNSummaryOut/sample_sheet.obj")
sample_summary <- read.csv("08_ReadNSummaryOut/summary_sheet.csv", row.names = 1)


# -------------------------------------------- #
# Import to phyloseq
# -------------------------------------------- #
# Preparetion to import to phyloseq
dim(sample_sheet); dim(tax_claident2); dim(seqtab_nochim)
all(rownames(sample_sheet) == rownames(seqtab_nochim)) # sample name check
all(sample_sheet$Sample_Name2 == rownames(seqtab_nochim)) # sample name check
rownames(sample_sheet) <- rownames(seqtab_nochim)

# Compile taxa information
colnames(seqtab_nochim) <- rownames(tax_claident2) # change col name
all(rownames(tax_claident2) == colnames(seqtab_nochim)) # taxa name check
tax_claident2$std_or_field <- "Field DNA"
tax_claident2[substr(tax_claident2$species, 1, 4) == "STD_",]$std_or_field <- "Standard DNA"

# Generating table for manual NC check
sample_sheet_na <- matrix(rep(NaN, ncol(tax_claident2)*ncol(sample_sheet)), ncol=ncol(sample_sheet))
colnames(sample_sheet_na) <- colnames(sample_sheet)
sample_sheet_comb <- rbind(sample_sheet_na, as.matrix(sample_sheet))
tax_table_comb <- rbind(t(as.matrix(tax_claident2)), seqtab_nochim)
table_for_nc_check <- cbind(sample_sheet_comb, tax_table_comb)
table_for_nc_check[,"std_validity"][1:ncol(tax_claident2)] <- colnames(tax_claident2)
write.csv(table_for_nc_check, sprintf("%s/SummarizedTable_for_NCcheck.csv", output_folder))

# Import data to phyloseq
ps0 <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows=FALSE),
                sample_data(sample_sheet),
                tax_table(as.matrix(tax_claident2)))

# Visualize pattern
ps_m3 <- melt(sample_summary, measure.vars = c("STD_all", "NonSTD_all"),
              id.vars = c("sample_nc"))

p1 <- ggplot(ps_m3, aes(x = sample_nc, y = value, color = variable, group = interaction(sample_nc, variable))) +
  geom_boxplot(colour = "black", outlier.colour = "white", outlier.size = 0) +
  scale_color_manual(values = c("red3", "royalblue"), name = NULL, labels = c("Standard DNA", "Non-standard DNA")) +
  geom_jitter(alpha = 0.5, position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0.8)) +
  ylab("Sequence reads/sample") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "top") +
  scale_x_discrete(limits = c("sample", "std_nc", "field_nc", "pcr_nc"),
                   labels = c("field_nc" = "Field NC", "pcr_nc" = "PCR NC",
                              "sample" = "Sample", "std_nc" = "Standard NC"))
ggsave(sprintf("%s/SequenceReads_Boxplot.pdf", output_folder), plot = p1, width = 6, height = 6)
saveRDS(p1, "../../../07_FormatFigs/data_robj/SeqReadsOverview_16S.obj")

# Examination of ASVs that were detected in field NC, PCR NC and standard NC
# Extract taxa that are detected from NC samples
# Barplot of NC-detected taxa (Optional)
ps_nctax <- prune_taxa(taxa_sums(prune_samples(sample_data(ps0)$sample_nc != "sample" & sample_data(ps0)$sample_nc != "pc", ps0)) > 0, ps0)
p2 <- plot_bar(ps_nctax, x = "sample_nc", fill = "phylum") +
  geom_bar(stat = "identity", colour = NA) + scale_fill_igv()

# Boxplot for NC-detected taxa
ps_m4 <- psmelt(ps_nctax)
ps_m5 <- ps_m4[ps_m4$sample_nc != "sample",]
ps_m5 <- ps_m5[ps_m5$std_or_field == "Field DNA" & ps_m5$Abundance > 0,]
p3 <- ggplot(ps_m5, aes(x = sample_nc, y = Abundance, group = interaction(phylum, sample_nc), fill = phylum)) +
  geom_boxplot(width = 0.6, position = position_dodge(width=0.7)) + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_igv() + xlab(NULL) + ylab("Log10(Sequence reads)") +
  scale_x_discrete(limits = c("field_nc", "std_nc", "pcr_nc"),
                   labels = c("field_nc" = "Field NC",
                              "pcr_nc" = "PCR NC",
                              "std_nc" = "Standard NC"))
ggsave(sprintf("%s/SequenceReads_NCtaxaBoxplot.pdf", output_folder), plot = p3, width = 8, height = 6)


# Examination of field negative controls
# (Time series plot)
# (Individual sample plot: Field negative controls)
ps_nctax2 <- prune_taxa(taxa_sums(prune_samples(sample_data(ps0)$sample_nc == "field_nc", ps0)) > 0, ps0)
ps_m7 <- psmelt(ps_nctax2)
if(class(ps_m7$date) != "Date") ps_m7$date <- ymd(ps_m7$date)
ps_m8 <- ps_m7[ps_m7$sample_nc == "field_nc" & ps_m7$std_or_field == "Field DNA",]
ps_m9 <- ps_m7[ps_m7$sample_nc == "field_nc" & ps_m7$std_or_field == "Standard DNA",]
ps_m10 <- ps_m7[ps_m7$std_or_field == "Field DNA",]
ps_m11 <- ps_m7[ps_m7$std_or_field == "Standard DNA",]

p4 <- ggplot(ps_m8, aes(x = date, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_igv()

p5 <- ggplot(ps_m9, aes(x = date, y = Abundance, fill = species)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_igv()

field_nc_reads <- plot_grid(p4, p5, ncol = 2, labels = "auto", align = "hv")
ggsave(sprintf("%s/SequenceReads_FieldNCreads.pdf", output_folder), plot = field_nc_reads, width = 16, height = 6)


# (Individual sample plot: Positive control)
p6 <- ggplot(ps_m10, aes(x = Sample, y = Abundance, colour = phylum, fill = phylum)) +
  geom_bar(stat = "identity", colour = NA) +
  scale_fill_igv() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 1))

p7 <- ggplot(ps_m11, aes(x = Sample, y = Abundance, colour = species, fill = species)) +
  geom_bar(stat = "identity", colour = NA) +
  scale_fill_igv() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 1))

sample_reads <- plot_grid(p6, p7, ncol = 2, align = "hv")
ggsave(sprintf("%s/SequenceReads_Samplereads.pdf", output_folder), plot = sample_reads, width = 16, height = 6)


# Examination of PCR negative controls and standard negative controls
# (Individual sample plot: PCR negative controls)
ps_nctax3 <- prune_taxa(taxa_sums(prune_samples(sample_data(ps0)$sample_nc == "std_nc" | sample_data(ps0)$sample_nc == "pcr_nc", ps0)) > 0, ps0)
ps_m12 <- psmelt(ps_nctax3)
if(class(ps_m12$date) != "Date") ps_m12$date <- ymd(ps_m12$date)
ps_m13 <- ps_m12[ps_m12$sample_nc == "pcr_nc" & ps_m12$std_or_field == "Field DNA",]
ps_m14 <- ps_m12[ps_m12$sample_nc == "pcr_nc" & ps_m12$std_or_field == "Standard DNA",]
ps_m15 <- ps_m12[ps_m12$sample_nc == "std_nc" & ps_m12$std_or_field == "Field DNA",]
ps_m16 <- ps_m12[ps_m12$sample_nc == "std_nc" & ps_m12$std_or_field == "Standard DNA",]

p10 <- ggplot(ps_m13, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_igv() + ggtitle("PCR negative controls: Field DNA reads")

p11 <- ggplot(ps_m14, aes(x = Sample, y = Abundance, fill = species)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_igv() +
  ggtitle("PCR negative controls: Standard DNA reads")

pcr_nc_reads <- plot_grid(p10, p11, ncol = 2, labels = "auto", align = "hv")
ggsave(sprintf("%s/SequenceReads_PCRNCreads.pdf", output_folder), plot = pcr_nc_reads, width = 16, height = 8)


# (Individual sample plot: Standard negative controls)
p12 <- ggplot(ps_m15, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity", colour = NA) + theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_igv() +
  ggtitle("Standard negative controls: Field DNA reads")

p13 <- ggplot(ps_m16, aes(x = Sample, y = Abundance, fill = species)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = 90)) + scale_fill_igv() +
  ggtitle("Standard negative controls: Standard DNA reads")

std_nc_reads <- plot_grid(p12, p13, ncol = 2, labels = "auto", align = "hv")
ggsave(sprintf("%s/SequenceReads_STDNCreads.pdf", output_folder), plot = std_nc_reads, width = 16, height = 8)


# Calculating ratio
dna_summary <- aggregate(sample_summary$NonSTD_all, list(sample_summary$sample_nc), mean)
dna_summary$prop <- dna_summary$x/dna_summary[dna_summary$Group.1 == "sample","x"]
write.csv(dna_summary, sprintf("%s/DNAreads_ContaminationLevel.csv", output_folder))

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/09_SessionInfo_NCcheck_%s.txt", substr(Sys.time(), 1, 10)))

