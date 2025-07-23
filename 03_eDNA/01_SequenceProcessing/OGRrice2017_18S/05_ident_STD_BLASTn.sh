##
## Standard STD identification
##

cd "/XXXXX/OGRrice2017_18S/"

DBPATH=STDseqs/EukSTD_Euk1391f
QUERYPATH=03_OTUClusteringOut/OTU_seqs.fa
OUTPUT=05_ident_STD_BLASTnOut/STD_out.txt
EVALUE_SET=1e-50

mkdir 05_ident_STD_BLASTnOut
blastn -db ${DBPATH} -query ${QUERYPATH} -evalue ${EVALUE_SET} -outfmt 6 -out ${OUTPUT}

