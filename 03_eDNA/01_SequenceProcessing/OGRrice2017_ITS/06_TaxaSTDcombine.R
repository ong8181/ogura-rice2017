####
#### OGR rice 2017 eDNA analysis
#### Assign STD
####

# Load library
library(tidyverse); packageVersion("tidyverse") # 2.0.0

# Set random seeds (for reproduction)
ran.seed <- 8181
set.seed(ran.seed)

# Create output directory
wdir <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(wdir, end = -3), "Out")); rm(wdir)
dir.create(output_folder)

# Load Claident taxa assignment
claident_tax <- read.delim("04_TaxaAssignmentOut/final_classigntax_otu")
head(claident_tax)

# Load Blastn taxa assignment for STD sequences
blastn_std0 <- read.table("05_ident_STD_BLASTnOut/STD_out.txt")
cond1 <- blastn_std0$V5 < 2 & blastn_std0$V4 > 220 & blastn_std0$V6 < 1 # Alignment length and No. of mismatch and gap
cond2 <- blastn_std0$V7 < 3 & blastn_std0$V8 > 220 # Start and end points of query
blastn_std <- blastn_std0[cond1 & cond2,]

# Check claident taxa assignments of the potential std sequences
potential_std_id <- match(blastn_std$V1, claident_tax$query)
claident_tax[potential_std_id, "family"] # Not close to any existing tax --> OK
claident_tax[potential_std_id, "species"] # Not close to any existing tax --> OK

# Replace STD taxa names with claident taxa assigment
claident_tax$species <- as.character(claident_tax$species)
claident_tax[potential_std_id, "species"] <- as.character(blastn_std$V2)

# Output new claident tax table
write.csv(claident_tax, sprintf("%s/claident_tax_revise.csv", output_folder), row.names = F)

#### save session info
writeLines(capture.output(sessionInfo()),
           sprintf("00_SessionInfo/06_SessionInfo_STDident_%s.txt", substr(Sys.time(), 1, 10)))

