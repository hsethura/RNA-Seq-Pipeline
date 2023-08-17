#!/bin/sh

MOC_ID=$1
shift
CONFIG_FILE=$1
shift
USERID=$1
shift

ID=`id | cut -d"(" -f2 | cut -d")" -f1`

KEY_DIR="/broad/IDP-Dx_storage/MOC/Key_files/"
KEY_FILE=$KEY_DIR"/"$MOC_ID"_key.txt"
SCRIPTS_DIR="/broad/IDP-Dx_storage/MOC/scripts/"
FILE_DIR="/broad/IDP-Dx_storage/MOC/files/"
ASPERA_MAKER="/broad/IDP-Dx_storage/MOC/scripts/Aspera_manifest_maker.sh"

TEMP_FILE=~/"temp_"$MOC_ID"_transfer_.txt"


# typeset -i SAMP_ID_FIELD
# typeset -i CG_ID_FIELD


##############  FUNCTIONS  #################

# Pulling out paths, header names from config file 
	config_read ()
	{
	
		CONFIG_FILE=$1
		VAR_NAME=$2
		cat $CONFIG_FILE | grep $VAR_NAME":" | awk '{print $2}'
	
	
	
	}

FIELD_HEADER ()
{

	KEY_FILE=$1
	HEADER_NAME=$2

	cat $KEY_FILE | sed 's/ /_/g' | grep -v "#" | head -1 | awk -F"\t" -v HEADER_NAME=$HEADER_NAME '{	
																for(i=1; i < NF+1; i++)															
																	if($i==HEADER_NAME)
																		print i
															}' 
}

FIND_FIELD ()
{
	
	if [ $# == 4 ];then
		FILE=$1
		CGF=$2
		GCID=$3
		SIF=$4
				cat $FILE | sed 's/ /_/g' | grep -v "#" | sed 1d | awk -F"\t" -v GCID=$GCID -v CGF=$CGF -v SIF=$SIF '{	

														if($CGF==GCID)
															printf "%s,", $SIF
													}' | sed 's/,$//g'
	fi
	if [ $# == 2 ];then
		FILE=$1
		SIF=$2
				
				cat $FILE | sed 's/ /_/g' | grep -v "#" | sed 1d  | awk -F"\t" -v SIF=$SIF '{	

															print $SIF
													}' 
	fi 


} 

########  function for extracting options from command line
 
		extract_option ()
		{
			option_name=$1 
			shift
			default=$1
			shift
			
			echo $@ | awk -v def=$default -v name=$option_name '{
		
				for(i=1; i < NF+1; i++)
				{
					if(i==1)			# intitalize to $default in case no option found
						print def
					if($i==name)
						print $(i+1)
				}
		
			}' | tail -1 
			
			
		}

####################################################

options=$@

HOMOL_OW=`extract_option -RESULTS_PATH - $options`


SAMP_ID_FIELD=`FIELD_HEADER $KEY_FILE Sample_ID`
PID_FIELD=`FIELD_HEADER $KEY_FILE Project_ID`
ALL_PROJ_ID=`FIND_FIELD $KEY_FILE $PID_FIELD | sort | uniq` 

if [ $RESULTS_PATH == "-" ];then
	RESULTS_PATH=`config_read $CONFIG_FILE Results_path`
fi

TRANSFER_PATH=`config_read $CONFIG_FILE Transfer_path`
BAM_PATH=`config_read $CONFIG_FILE Bam_path`

for PROJ_ID in $ALL_PROJ_ID
do


	if [ $USERID == "Y" ];then
		RES_DIR=$RESULTS_PATH"/"$ID"/"$PROJ_ID"/"
		BAM_DIR=$BAM_PATH"/"$ID"/"$PROJ_ID"/"
	else
		RES_DIR=$RESULTS_PATH"/"$PROJ_ID"/"
		BAM_DIR=$BAM_PATH"/"$PROJ_ID"/"
	fi 
	TRANS_DIR=$TRANSFER_PATH"/"$PROJ_ID"_Results/"
	TAR_DIR=$TRANSFER_PATH"/"$PROJ_ID"_Results.tar.gz"

	GV_DIR=$TRANS_DIR"/GV/"
	
	echo $TRANSFER_PATH
	echo $RESULTS_PATH
	echo $BAM_DIR
	echo $GV_DIR

	rm -r $TRANS_DIR

	mkdir -p $TRANS_DIR
	mkdir -p $GV_DIR

	GV_file_paths=`cat $RES_DIR/*"paths_abs.txt" | grep -e tdf -e gff -e fna`


	echo "Copying results files to "$TRANS_DIR

	cp $RES_DIR/*".tsv" $TRANS_DIR
	cp $RES_DIR/*"AlignmentSummaryMetrics.txt" $TRANS_DIR
	cp $RES_DIR/*"etrics.txt" $TRANS_DIR
	cp -r $RES_DIR/"DE" $TRANS_DIR


	echo "Copying GV files to "$GV_DIR
	cp $GV_file_paths $GV_DIR

	echo "Copying user guide files to "$TRANS_DIR
	cp $FILE_DIR"/"*"CURRENT"* $TRANS_DIR

	ls -lrt $TRANS_DIR

	tar czf $TAR_DIR $TRANS_DIR

	sh $ASPERA_MAKER $BAM_DIR bam $TAR_DIR

	echo $TAR_DIR
done