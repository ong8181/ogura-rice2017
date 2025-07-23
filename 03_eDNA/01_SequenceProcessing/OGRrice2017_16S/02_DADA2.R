####
#### Rice field samples (Ogura, 2017)
#### Prokaryote 16S sequence analysis by DADA2
####

# Set random seeds (for reproduction)
ran.seed <- 8181
set.seed(ran.seed)
dir.create("00_SessionInfo")
output_folder <- "02_DADA2Out"
dir.create(output_folder)

# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2024.07.05
library(dada2); packageVersion("dada2") # 1.28.5, 2024.07.05
library(ShortRead); packageVersion("ShortRead") # 1.58.0, 2024.07.05

# Load sequence reads
# RMR-080-Pro
path <- "/XXXXX/OGRrice2017_16S/01_PrimerRemoveCheckOut/"
fnFs <- sort(list.files(path, pattern="R1.fastq.gz", full.names = T)) # Forward read files
fnRs <- sort(list.files(path, pattern="R2.fastq.gz", full.names = T)) # Reverse read files
# Get sample names, assuming files named as so: SAMPLENAME_XXX.fastq
(sample_names <- sapply(strsplit(fnFs, "_"), `[`, 8))

# Visualize quality# Visualize quality
#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])

# Perform quality filtering
filtFs <- file.path(path, "filtered", paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample_names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     # Output sequences of Claident already trimmed Ns and primers
                     truncLen = c(200, 200), # For 515F-806R
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=F,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
#plotQualityProfile(filtFs[1:2])
#plotQualityProfile(filtRs[1:2])

# Exclude 0 seq samples, rename filtFs and filtRs
if(length(sample_names[out[,2]<1 | out[,1]<1]) > 0){
  filtFs <- file.path(path, "filtered", paste0(sample_names[out[,2]>0 & out[,1]>0], "_F_filt.fastq.gz"))
  filtRs <- file.path(path, "filtered", paste0(sample_names[out[,2]>0 & out[,1]>0], "_R_filt.fastq.gz"))
}

# Learn the error rates
min_nbases <- 200 * sum(out[,2])
errF <- learnErrors(filtFs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20, nbases = min_nbases)
errR <- learnErrors(filtRs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20, nbases = min_nbases)

# Visualize errors
ggsave(sprintf("%s/errF.pdf", output_folder), plotErrors(errF, nominalQ = T), width = 10, height = 10)
ggsave(sprintf("%s/errR.pdf", output_folder), plotErrors(errR, nominalQ = T), width = 10, height = 10)

# Dereplicatin
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample_names[out[,2] > 0 & out[,1] > 0]
names(derepRs) <- sample_names[out[,2] > 0 & out[,1] > 0]

# Sample inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = TRUE)

# Merging paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE, trimOverhang = TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab); sum(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# Cutting unexpected length sequences
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(230,270)] # For 16S
#seqtab2 <- seqtab
table(nchar(getSequences(seqtab2)))

# Remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab2, method = "consensus", multithread = TRUE, verbose = TRUE)
table(nchar(getSequences(seqtab_nochim)))
dim(seqtab_nochim)
sum(seqtab_nochim)/sum(seqtab2)

# Track reads thourhg the pipeline
out2 <- out[out[,2]>0 & out[,1]>0,]
getN <- function(x) sum(getUniques(x))
track <- cbind(out2, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab2), rowSums(seqtab_nochim),  rowSums(seqtab_nochim)/out2[,1])
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "tabled2", "nonchim", "prop(last/first)")
rownames(track) <- sample_names[out[,2]>0 & out[,1]>0]
head(track); hist(track[,"prop(last/first)"])

# Taxa output for claident tax assignment
seqs <- colnames(seqtab_nochim)
seqs_out <- as.matrix(c(rbind(sprintf(">ASV%05d", 1:length(seqs)), seqs)), ncol = 1)
seqtab_nochim_for_csv <- seqtab_nochim
colnames(seqtab_nochim_for_csv) <- sprintf("ASV%05d", 1:length(seqs))

# Save outputs
write.csv(seqs, paste0(output_folder, "/seq_only.csv"), row.names = colnames(seqtab_nochim))
write.table(seqs_out, paste0(output_folder, "/ASV_seqs.fa"), col.names = FALSE, row.names = FALSE, quote = FALSE)
write.csv(seqtab_nochim_for_csv, paste0(output_folder, "/seqtab_nochim.csv"))
write.csv(track, paste0(output_folder, "/track.csv"))

# Save workspace
rm(derepFs)
rm(derepRs)
rm(dadaFs)
rm(dadaRs)
save(list = ls(all.names = TRUE),
     file = paste0(output_folder, "/DADA2Out.RData"))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("00_SessionInfo/02_DADA2Out_", substr(Sys.time(), 1, 10), ".txt"))
