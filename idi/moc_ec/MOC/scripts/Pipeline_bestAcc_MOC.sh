#!/bin/sh

MOC_ID=$1

source "/broad/software/scripts/useuse"
reuse -q UGER
reuse -q Java-1.6
reuse -q BWA
reuse -q BLAST+

### source all functions and set variables 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"
paths_and_headers $MOC_ID $@

 
USER=`id | awk '{print $1}' | cut -d"(" -f2 | sed 's/)//g'` 

STRAIN=`extract_option -STRAIN - 1 $@`
USERID=`extract_option -USERID $USER 1 $@`


############  Extract field IDs from key headers  #################
 

TEMP_DIR=$TEMP_PATH"/"$USERID"/best_acc/"

mkdir -p $TEMP_DIR

echo $TEMP_DIR 

ls -lrt $KEY_FILE

exit

echo $CONFIG_FILE

INDEX_F=`FIELD_CONFIG_HEADER $KEY_FILE $CONFIG_FILE P7`
GENUS_F=`FIELD_CONFIG_HEADER $KEY_FILE $CONFIG_FILE Genus`
STRAIN_F=`FIELD_CONFIG_HEADER $KEY_FILE $CONFIG_FILE Strain `
BC_F=`FIELD_CONFIG_HEADER $KEY_FILE $CONFIG_FILE bc` 
PROJID_F=`FIELD_CONFIG_HEADER $KEY_FILE $CONFIG_FILE Proj`
SAMPID_F=`FIELD_CONFIG_HEADER $KEY_FILE $CONFIG_FILE ID`
REF_F=`FIELD_CONFIG_HEADER $KEY_FILE $CONFIG_FILE Ref_accession`

ALL_STRAIN_FS=`FIND_VALUES_FOR_HEADER Strain_ID $KEY_FILE | sort | uniq`
ALL_GENUS_FS=`FIND_VALUES_FOR_HEADER Genus $KEY_FILE | sort | uniq`

echo $ALL_STRAIN_FS

echo $REF_PATH


MERGE_DIR=$TEMP_PATH"/"$MOC_ID"/"*"/mergedir"

# check if the fastqs are gz
GZ=`ls -lrt $MERGE_DIR | awk '{print $NF}' | grep fastq | grep .gz | wc -l`

# get file paths and samp ids for all strains
for m in $ALL_STRAIN_FS
do
	FILES=`cat $KEY_FILE | grep -v "###"  | grep $m | head -1 | awk -F"\t" -v GZ=$GZ '{
	
				split($'$REF_F', ar, ";")
				REF=ar[1]
				
				if(GZ == 0)
					FILES="'$MERGE_DIR'/"$'$SAMPID_F'"_R2.fastq"
				else
					FILES="'$MERGE_DIR'/"$'$SAMPID_F'"_R2.fastq.gz"
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

	REF_V=`echo $m | cut -d";" -f1`
	GENUS_V=`echo $m | cut -d";" -f2`
	STRAIN_V=`echo $m | cut -d";" -f3` 
	SEQ_FILE=`echo $m | cut -d";" -f4`		
	
	echo $SEQ_FILE 
		
	echo "sh $MOC_SCRIPT_PATH"Pipeline_align_test_MOC.sh" $STRAIN_V $GENUS_V 10000 $SEQ_FILE $ALN_DIR $TEMP_DIR $REF_PATH $MOC_SCRIPT_PATH"
	sh $MOC_SCRIPT_PATH"Pipeline_align_test_MOC.sh" $STRAIN_V $GENUS_V 10000 $SEQ_FILE $ALN_DIR $TEMP_DIR $REF_PATH $MOC_SCRIPT_PATH
done

for GENUS in $ALL_GENUS_FS
do
	echo $GENUS
	STRAINS=`cat $KEY_FILE | grep -v "###" | awk -F"\t" -v GENUS=$GENUS '{
	
				if($'$GENUS_F'==GENUS)
					print $'$STRAIN_F'
				
		
			}' | sort | uniq` 

	sh $MOC_SCRIPT_PATH"Best_Acc_pick.sh" $GENUS $STRAINS
done



