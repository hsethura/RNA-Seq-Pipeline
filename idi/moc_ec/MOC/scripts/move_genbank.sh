#!/bin/sh

OUTPUT_PATH="/idi/moc_ec/MOC/genbank/"

mkdir -p $OUTPUT_PATH

cd $OUTPUT_PATH
#rm $OUTPUT_PATH"assembly_summary"*

ASSEMB_FILE=$OUTPUT_PATH"assembly_summary.txt"

# if [ ! -z $OUTPUT_PATH"assembly_summary.txt" ];then
# 	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
# fi


#ALL_SPEC=`cat $ASSEMB_FILE | awk -F"\t" '{print $8}' | awk '{print $1, $2}'| sed 's/ /_/g' | sort -T /broad/hptmp/ | uniq -c | sort -T /broad/hptmp/ | awk '{if($1 >10) print $2}'`

ALL_SPEC="Vibrio_parahaemolyticus Enterococcus_faecalis Mycobacterium_tuberculosis Shigella_flexneri Pseudomonas_aeruginosa Acinetobacter_baumannii Clostridioides_difficile Shigella_sonnei Enterococcus_faecium Klebsiella_pneumoniae Campylobacter_coli Streptococcus_pneumoniae Staphylococcus_aureus Listeria_monocytogenes Campylobacter_jejuni Escherichia_coli Salmonella_enterica" 

for SPEC in $ALL_SPEC
do
	
	ALL_DIRS=`cat assembly_summary.txt | sed 's/ /_/g' | grep $SPEC |  awk -F"\t" '{print $(NF-3)"/"}' | head -200`
	
	for DIR in $ALL_DIRS
	do	
		ACC=`echo $DIR | rev | cut -d"/" -f2 | rev`
		OUTPUT_DIR=$OUTPUT_PATH"/"$SPEC"/"$ACC
		echo $OUTPUT_DIR
		
		rm -r $OUTPUT_DIR
		mkdir -p $OUTPUT_DIR
		cd $OUTPUT_DIR
	
		FNA=$DIR$ACC"_genomic.fna.gz"
		GFF=$DIR$ACC"_genomic.gff.gz"
		RNA=$DIR$ACC"_rna_from_genomic.fna.gz"
		
		
# 		wget -q $FNA 2>&1 /dev/null
# 		wget -q $GFF 2>&1 /dev/null
		wget -q $RNA 2>&1 /dev/null

		echo $OUTPUT_DIR
		ls -lrt $OUTPUT_DIR
		cd
	done
done

ls -lrt $OUTPUT_PATH*