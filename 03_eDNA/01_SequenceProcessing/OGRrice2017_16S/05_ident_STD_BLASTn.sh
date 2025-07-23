##
## Standard STD identification
##

cd "/XXXXX/OGRrice2017_16S/"

DBPATH=STDseqs/ProkSTD_515F
QUERYPATH=03_OTUClusteringOut/OTU_seqs.fa
OUTPUT=05_ident_STD_BLASTnOut/STD_out.txt
EVALUE_SET=1e-100

mkdir 05_ident_ProkSTD_BLASTnOut
blastn -db ${DBPATH} -query ${QUERYPATH} -evalue ${EVALUE_SET} -outfmt 6 -out ${OUTPUT}

