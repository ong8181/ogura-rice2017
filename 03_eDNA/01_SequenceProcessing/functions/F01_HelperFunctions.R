####
#### F01. Collection of helper functions
####

# Primer orientation function
# From https://benjjneb.github.io/dada2/ITS_workflow.html
AllOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

# Count the number of primer hits
PrimerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}


# ggplot function 1
PlotStyle <-  function(ggobject){
  return(ggobject + theme_bw() + theme(axis.text.x = element_text(angle=0),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       axis.text = element_text(size=12),
                                       axis.title = element_text(size=12),
                                       panel.background=element_rect(colour="black", fill=NA, size=0.8)))
}

# ggplot function 2
PlotStyle2 <- function(ggobject){
  return(ggobject + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                          axis.title.x = element_blank(),
                          legend.position = "none") +
           geom_jitter(shape = 16, size = 2, alpha = 0.8, width = 0.1, height = 0))
}


# Merge standard DNA sequences
MergeSTD <- function(std_i, std_data = std_table){
  index_std <- which(match(colnames(std_data), std_i) == 1)
  #index_std <- which(colnames(std_data) == std_i)
  if(length(index_std) > 1){
    std_tmp <- rowSums(std_data[,index_std])
  }else{
    std_tmp <- std_data[,index_std]
  }
  return(std_tmp)
}

# Phylum sorting
PhylumSort <- function(phyloseq_object, select_n = 51){
  phylum_sort <- names(sort(taxa_sums(phyloseq_object), decreasing=TRUE))
  phylum_n <- length(phylum_sort)
  merged_phyloseq_object <- merge_taxa(phyloseq_object, phylum_sort[select_n:phylum_n])
  merged_phyloseq_object <- subset_taxa(merged_phyloseq_object, Phylum != "Not_Identified")
  return(merged_phyloseq_object)
}

