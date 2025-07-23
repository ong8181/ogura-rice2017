####
#### Quality filtering for Illumina data
####

# -------------------------------------------- #
# Step 0. Specify parameters
# -------------------------------------------- #
# Set a working directory
WORKING_DIR=/XXXXX/OGRrice2017_16S

# Specify primer sequences (forward...reverse_revcomp, reverse...forward_revcomp)
PRIMER_FOR_READ=^GTGYCAGCMGCCGCGGTAA...ATTAGAWACCCBNGTAGTCC
PRIMER_REV_READ=^GGACTACNVGGGTWTCTAAT...TTACCGCGGCKGCTGRCAC
## ----------- 515F-Y - 806R (RevComp)
# -a ^GTGYCAGCMGCCGCGGTAA...ATTAGAWACCCBNGTAGTCC \
# -A ^GGACTACNVGGGTWTCTAAT...TTACCGCGGCKGCTGRCAC \

# Set other parameters
INPUT_DIR=${WORKING_DIR}/seqdata_demultiplexed
OUTPUT_DIR=01_PrimerRemoveCheckOut

cd ${INPUT_DIR}
mkdir ../${OUTPUT_DIR}
mkdir ../${OUTPUT_DIR}/01_temp

# Check and save MD5
md5sum *.gz > 0_md5sum.txt


# -------------------------------------------- #
# Step 1.Trim primers
# -------------------------------------------- #
# ----------------- cutadapt version
# Single Primer removal
for file in *forward.fastq.gz; do
cutadapt -j 72 \
-a $PRIMER_FOR_READ \
-A $PRIMER_REV_READ \
-n 2 \
-o ../${OUTPUT_DIR}/01_temp/${file%_515F.forward.fastq.gz}_checked_R1.fastq.gz \
-p ../${OUTPUT_DIR}/01_temp/${file%_515F.forward.fastq.gz}_checked_R2.fastq.gz \
${file} \
${file%forward.fastq.gz}reverse.fastq.gz
done


# -------------------------------------------- #
# Step 2. Get summary table 
# -------------------------------------------- #
seqkit stats -a *.gz > ${WORKING_DIR}/${OUTPUT_DIR}/01_1_demultiplexed.txt
seqkit stats -a *forward.fastq.gz > ${WORKING_DIR}/${OUTPUT_DIR}/01_1_demultiplexedR1.txt
cd ../${OUTPUT_DIR}/01_temp
seqkit stats -a *.gz > ${WORKING_DIR}/${OUTPUT_DIR}/01_2_after_cutadapt.txt
seqkit stats -a *R1.fastq.gz > ${WORKING_DIR}/${OUTPUT_DIR}/01_2_after_cutadaptR1.txt


# -------------------------------------------- #
# Step 3. Clean up 
# -------------------------------------------- #
cd ..
mv 01_temp/*.gz ./
rm -r 01_temp
