#!/bin/sh


MOC_ID=$1

### source all functions 
source "/broad/IDP-Dx_storage/MOC/scripts/MOC_functions.sh"

### get options from command line

# -fastq: include fastq files along with bams in Aspera (default N)
# -in: internal collaborator (default N)

CONFIG_FILE=`extract_option -conf "/broad/IDP-Dx_storage/MOC/config_files/Universal_config.yaml" 1 $@`
INTERNAL=`extract_option -in N 1 $@`

### run path_suff function to set RESPATH_SUFF.
##	If -moc_id N included in command line, do not moc_ID to RESPATH_SUFF.  
##	If -user_id included with Y or no string add USID to RESPATH_SUFF.  
##	If -user_id followed by userID, add userID to RESPATH_SUFF.

path_suff $@

### get paths etc. from config file
KEY_DIR=`config_read $CONFIG_FILE Key_base`
USER_GUIDE=`config_read $CONFIG_FILE user_guide`
ASPERA_MAKER=`config_read $CONFIG_FILE apera_script`
RESULTS_PATH=`config_read $CONFIG_FILE Results_path`
TRANSFER_PATH=`config_read $CONFIG_FILE Transfer_path`
BAM_PATH=`config_read $CONFIG_FILE Bam_path`
TEMP_PATH=`config_read $CONFIG_FILE Temp_path`
DB_SCRIPT=`config_read $CONFIG_FILE DB_script`
MAN_DIR=`config_read $CONFIG_FILE manifest_dir`

### set paths to files
KEY_FILE=$KEY_DIR"/"$MOC_ID"_key.txt"
TEMP_FILE=~/"temp_"$MOC_ID"_transfer_.txt"

ALL_SMOC_INDEX1=`sh $DB_SCRIPT SMOC_ID,$SMOC_ID Index1_Seq | sed 1d | awk '{print $2}'  | sed 's/,/ /g'`

SAMP_ID_FIELD=`FIELD_HEADER $KEY_FILE Sample_ID`
PID_FIELD=`FIELD_HEADER $KEY_FILE Project_ID`
ALL_PROJ_ID=`FIND_FIELD $KEY_FILE $PID_FIELD | sort | uniq` 

for PROJ_ID in $ALL_PROJ_ID
do
	RES_DIR=$RESULTS_PATH$RESPATH_SUFF$PROJ_ID"/"
	BAM_DIR=$BAM_PATH$RESPATH_SUFF$PROJ_ID"/"
	MERGE_DIR=$TEMP_PATH$RESPATH_SUFF"/*/mergedir/"
	
	TRANS_DIR=$TRANSFER_PATH$RESPATH_SUFF$PROJ_ID"_Results/"
	TAR_DIR=$TRANSFER_PATH$MOC_ID"_"$PROJ_ID"_Results.tar.gz"

	echo $TAR_DIR
	echo $TRANS_DIR

	GV_DIR=$TRANS_DIR"/GV/"
	
	if [ ! -s $RESULTS_PATH ];then
	
		echo $RESULTS_PATH" not found!"
		exit
		
	fi
	
	echo $RESULTS_PATH
	echo $TRANSFER_PATH
	echo $BAM_DIR
	echo $GV_DIR
	
	rm -r $TRANS_DIR

	mkdir -p $TRANS_DIR
	mkdir -p $GV_DIR

	GV_file_paths=`cat $RES_DIR/*"paths_abs.txt" | grep -e tdf -e gff -e fna`


	echo "Copying results files from "$RES_DIR " to "$TRANS_DIR

	cp $RES_DIR/*".tsv" $TRANS_DIR
	cp $RES_DIR/*"AlignmentSummaryMetrics.txt" $TRANS_DIR
	cp $RES_DIR/*"etrics.txt" $TRANS_DIR
	cp $RES_DIR/*"etrics.txt" $TRANS_DIR
	cp $RES_DIR/*"_abs.txt" $TRANS_DIR
	cp $RES_DIR/*"_gv.txt" $TRANS_DIR
	cp $RES_DIR/*"session_info.txt" $TRANS_DIR
	
	if [ -s $RES_DIR/"DE" ];then
		cp -r $RES_DIR/"DE" $TRANS_DIR
	fi

	echo "Copying GV files to "$GV_DIR
	cp $GV_file_paths $GV_DIR

	echo "Copying user guide files to "$TRANS_DIR
	cp $USER_GUIDE $TRANS_DIR

	ls -lrt $TRANS_DIR
	
	echo "Compressing "$TRANS_DIR
	echo 	"-zczf $TAR_DIR $RES_DIR"
	cd $TRANSFER_PATH
	tar -zczf $TAR_DIR $RES_DIR
	
	# for external collaborators, create an aspera request email
	
	if [ $INTERNAL == "N" ];then
		echo "sh $ASPERA_MAKER $MOC_ID $BAM_DIR $TAR_DIR $MERGE_DIR bam,tdf $@"
		sh $ASPERA_MAKER $MOC_ID $BAM_DIR $TAR_DIR $MERGE_DIR bam,tdf $@
	fi 

	# for internal collaborator, generate email with paths to files 

	if [ $INTERNAL == "Y" ];then

		mkdir -p $MAN_DIR

		MESSAGE_FILE=~/"paths_file.txt"

		echo "Dear Collaborator," > $MESSAGE_FILE
		echo ""  >> $MESSAGE_FILE
		echo "Below are paths to your fastq and bam files for "$MOC_ID" "$PROJ_ID >> $MESSAGE_FILE
		echo "These files will be deleted 7 days from delivery so please move them to a permanent location before then." >> $MESSAGE_FILE
		echo ""  >> $MESSAGE_FILE
		echo "Fastq dir: "$MERGE_DIR >> $MESSAGE_FILE
		echo "Bam dir: "$BAM_DIR >> $MESSAGE_FILE
		echo ""  >> $MESSAGE_FILE
		echo "Output of the analysis pipeline can be found here: " >> $MESSAGE_FILE
		echo $RES_DIR >> $MESSAGE_FILE
		echo ""  >> $MESSAGE_FILE
		echo "A tar archive of this directory for faster downloading can be found here: " >> $MESSAGE_FILE
		echo $TAR_DIR >> $MESSAGE_FILE
		echo ""  >> $MESSAGE_FILE
		echo ""  >> $MESSAGE_FILE
		echo "Best regards," >> $MESSAGE_FILE
		echo ""  >> $MESSAGE_FILE
		echo "The MOC team" >> $MESSAGE_FILE

		cat $MESSAGE_FILE |  mailx -s "Paths to files for $MOC_ID $PROJ_ID" $USID"@broadinstitute.org"

		echo $TAR_DIR
	
		change_perms $TRANS_DIR $TAR_DIR

	fi
done