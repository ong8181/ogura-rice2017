##
## Standard STD identification
##

cd "/XXXXX/OGRrice2017_COI/"

DBPATH=STDseqs/COISTD_mlCOIintF
QUERYPATH=03_OTUClusteringOut/OTU_seqs.fa
OUTPUT=05_ident_STD_BLASTnOut/STD_out.txt
EVALUE_SET=1e-167

mkdir 05_ident_STD_BLASTnOut
blastn -db ${DBPATH} -query ${QUERYPATH} -evalue ${EVALUE_SET} -outfmt 6 -out ${OUTPUT}

