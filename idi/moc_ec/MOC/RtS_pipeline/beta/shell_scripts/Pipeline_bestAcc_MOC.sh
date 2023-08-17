#!/bin/sh

MOC_ID=$1

source "/broad/software/scripts/useuse"
reuse -q UGER
reuse -q Java-1.6
reuse -q BWA
reuse -q BLAST+

### source all functions and set variables 
source "/broad/IDP-Dx_storage/MOC/scripts/MOC_functions.sh"
paths_and_headers $MOC_ID $@

 
USER=`id | awk '{print $1}' | cut -d"(" -f2 | sed 's/)//g'` 

STRAIN=`extract_option -STRAIN - 1 $@`
USERID=`extract_option -USERID $USER 1 $@`


############  Extract field IDs from key headers  #################
 

TEMP_DIR=$TEMP_PATH"/"$USERID"/best_acc/"

mkdir -p $TEMP_DIR

echo $TEMP_DIR 

ls -lrt $KEY_FILE


INDEX_F=`FIELD_CONFIG_HEADER $KEY_FILE $CONFIG_FILE P7`
GENUS_F=`FIELD_CONFIG_HEADER $KEY_FILE $CONFIG_FILE Genus`
STRAIN_F=`FIELD_CONFIG_HEADER $KEY_FILE $CONFIG_FILE Strain `
BC_F=`FIELD_CONFIG_HEADER $KEY_FILE $CONFIG_FILE bc` 
PROJID_F=`FIELD_CONFIG_HEADER $KEY_FILE $CONFIG_FILE Proj`
SAMPID_F=`FIELD_CONFIG_HEADER $KEY_FILE $CONFIG_FILE ID`
REF_F=`FIELD_CONFIG_HEADER $KEY_FILE $CONFIG_FILE Ref_accession`

ALL_STRAIN_FS=`FIND_VALUES_FOR_HEADER Strain_ID $KEY_FILE | sort | uniq`

MERGE_DIR=$TEMP_PATH"/"$MOC_ID"/"*"/mergedir"


for m in $ALL_STRAIN_FS
do
	FILES=`cat $KEY_FILE | grep -v "###"  | grep $m | head -1 | awk -F"\t" '{
	
				split($'$REF_F', ar, ";")
				REF=ar[1]
				
				FILES="'$MERGE_DIR'/"$'$SAMPID_F'"_R2.fastq"
				#FILES="'$BAM_DIR'/"$'$PROJID_F'"/"$'$SAMPID_F'"_R2.fastq"
				print REF";"$'$GENUS_F'";"$'$STRAIN_F'";"FILES				
				#print $'$BC_F'
		
			}'` 
	SPLIT_FILES=$SPLIT_FILES" "$FILES
done	 

echo $SPLIT_FILES


for m in $SPLIT_FILES
do
	echo $m

	REF_F=`echo $m | cut -d";" -f1`
	GENUS_F=`echo $m | cut -d";" -f2`
	STRAIN_F=`echo $m | cut -d";" -f3` 
	SEQ_FILE=`echo $m | cut -d";" -f4`		
	
	echo $SEQ_FILE 
		
	echo "sh $MOC_SCRIPT_PATH"Pipeline_align_test_MOC.sh" $STRAIN_F $GENUS_F 10000 $SEQ_FILE $ALN_DIR $TEMP_DIR $NCBI_PATH $MOC_SCRIPT_PATH"
	sh $MOC_SCRIPT_PATH"Pipeline_align_test_MOC.sh" $STRAIN_F $GENUS_F 10000 $SEQ_FILE $ALN_DIR $TEMP_DIR $NCBI_PATH $MOC_SCRIPT_PATH
done

 
