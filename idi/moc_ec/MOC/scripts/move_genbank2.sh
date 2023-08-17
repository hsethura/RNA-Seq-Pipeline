#!/bin/sh

ALL_ID_FILE=$1

source /idi/moc_ec/MOC/scripts/bash_header

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"


### run path_suff function to set RESPATH_SUFF.
##	If -moc_id N included in command line, do not moc_ID to RESPATH_SUFF.  
##	If -user_id included with Y or no string add USID to RESPATH_SUFF.  
##	If -user_id followed by userID, add userID to RESPATH_SUFF.

path_suff $@

### determining paths and headers 
### default config file is /broad/IDP-Dx_storage/MOC/config_files/PC_config.yaml
paths_and_headers $MOC_ID $@


CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/PC_config.yaml" 1 $@`
PHYLO=`extract_option -phylo S 1 $@`


wget_gb ()
{
	SPEC=$1
	shift
	OUTPUT_PATH=$1
	shift
	ALL_DIRS=$@
	
	NUM_GFF=0
	
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
	
	
		wget -q $FNA 2>&1 /dev/null
		wget -q $GFF 2>&1 /dev/null
		wget -q $RNA 2>&1 /dev/null

		GFF_FOUND=`ls -lrt $OUTPUT_DIR | grep "gff.gz" | wc -l`
		echo $GFF_FOUND
		
		if [ $GFF_FOUND -gt 0 ];then
			NUM_GFF=`expr $NUM_GFF + 1`
			ls -lrt $OUTPUT_DIR
		else 
			rm $OUTPUT_DIR
		fi
		echo $NUM_GFF
		if [ $NUM_GFF -eq 2 ];then
			
			break
		fi
	done
}


ALL_ID=`cat $ALL_ID_FILE | sed 's/ /_/g'`


echo $ALL_ID

OUTPUT_PATH="/idi/moc_ec/MOC/genbank/"

mkdir -p $OUTPUT_PATH

cd $OUTPUT_PATH
#rm $OUTPUT_PATH"assembly_summary"*

ASSEMB_FILE=$OUTPUT_PATH"assembly_summary.txt"

# if [ ! -z $OUTPUT_PATH"assembly_summary.txt" ];then
# 	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
# fi


#ALL_ID=`cat $ASSEMB_FILE | awk -F"\t" '{print $8}' | awk '{print $1, $2}'| sed 's/ /_/g' | sort -T /broad/hptmp/ | uniq -c | sort -T /broad/hptmp/ | awk '{if($1 >10) print $2}'`

SEG_FILE=$TEMP_DIR$"mvGB2_SGE_file.txt"

echo "" | sed 1d >  $SEG_FILE

for ID in $ALL_ID
do

	echo $ID
	echo "sh /idi/moc_ec/MOC/scripts/move_genbank3.sh $ID $OUTPUT_PATH -phylo G" > ~/temp.txt
	qsub ~/temp.txt >> $SEG_FILE					

done

echo "Jobs submitted - waiting for them to complete...."

echo $SEG_FILE
qstat

###testing that all jobs launched are finished
SGE_test $SEG_FILE	
echo "All Jobs completed...."


ls -lrt $OUTPUT_PATH*