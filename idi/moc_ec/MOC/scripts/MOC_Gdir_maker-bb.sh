#!/bin/sh


TYPE=$1
shift
ALL_IDS=$@ 


### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

### determining paths and headers 
paths_and_headers $ID $@

ALL_REFS_OPT=`extract_option -all_refs - 1 $@`
GREF=`extract_option -gref Y 1 $@`

### set paths to directories, files, and scripts from config file
### default config file is /idi/moc_ec/MOC/config_files/PC_config.yaml

CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/PC_config.yaml" 1 $@`

GPC_PATH=`config_read $CONFIG_FILE gdrivePC_path`
GPD_PATH=`config_read $CONFIG_FILE gdriveDEV_path`
GCID_PATH=`config_read $CONFIG_FILE gdriveGCIDTC_path`

GTRACK_PATH=`config_read $CONFIG_FILE gdrivetrack_path`
GPUNCH_PATH=`config_read $CONFIG_FILE gdrivepunch_path`
GLOCAL_PATH=`config_read $CONFIG_FILE gdrivelocal_path`
FILE_PATH=`config_read $CONFIG_FILE file_path`

GDRIVE_SCRIPT=`config_read $CONFIG_FILE gdrive_script`
TEMP_DIR=`config_read $CONFIG_FILE Temp_path`

check_and_copy ()
{
	TRAGET_NAME=$1
	GID=$2
	
	FLAG=`$GDRIVE_SCRIPT ls | grep $TRAGET_NAME | wc -l`
	if [ $FLAG == "0" ];then
		$GDRIVE_SCRIPT copy -id $GID $TRAGET_NAME
	fi

}

IMPORT_GS="Y"

########### Import google sheet databases to server ###############=
if [ $IMPORT_GS == "Y" ];then
	echo "Running $GSIMPORT_SCRIPT to import google sheet databases to server"
	echo "sh $GSIMPORT_SCRIPT"
	sh $GSIMPORT_SCRIPT
fi
#####################################################


for ID in $ALL_IDS
do
### Make local directory

	if [ $TYPE == "DP" ] || [ $TYPE == "TCP" ];then
		PROC=`echo $ID | cut -d"-" -f1`
	fi
	if [ $TYPE == "DE" ] || [ $TYPE == "TCE" ];then
		PROC=`echo $ID | cut -d"-" -f1`
		DPID=`echo $ID | cut -d"." -f1`
	fi
	
	
	if [ $TYPE == "P" ];then
		GLOCAL_DIR=$GLOCAL_PATH"/"$GPC_PATH"/"$ID
		GD_DIR=$GPC_PATH"/"$ID
		mkdir -p $GLOCAL_DIR
		cd $GLOCAL_DIR
		pwd
	fi
	
	if [ $TYPE == "TCP" ];then
		GLOCAL_DIR=$GLOCAL_PATH"/"$GCID_PATH"/"$PROC"/"
		mkdir -p $GLOCAL_DIR
		cd $GLOCAL_DIR
		DE_DB_NAME=$PROC"_Project_DB"
		GID="1k3FGDSBikNnMud0AMo2jqDA-2ZvuVBfUYFSaPqu8UfU"
		check_and_copy $DE_DB_NAME $GID
	
		GLOCAL_DIR=$GLOCAL_PATH"/"$GCID_PATH"/"$PROC"/"$ID
		mkdir -p $GLOCAL_DIR
		cd $GLOCAL_DIR

		mkdir -p $GLOCAL_DIR
		cd $GLOCAL_DIR
		pwd
	
		mkdir "Experiments"

		### Copy and rename Project Plan template to Gdrive dir (if not already there)
		echo "Copying and renaming Project Plan file template"
	
		PP_NAME=$ID"_Project_Plan"
		GID="1owfxuYJ0ev_WEOy496MNHn0J9XhoeBYwPwMsqTemlzI"
		check_and_copy $PP_NAME $GID

		DE_DB_NAME=$ID"_Expt_DB"
		GID="1ovPfuals4qsPfQgRvDXyT1uNCX81ohuTPFVSdAeAB1M"
		check_and_copy $DE_DB_NAME $GID


	fi

	if [ $TYPE == "DE" ] || [ $TYPE == "TCE" ];then
		if [ $TYPE == "DE" ];then
			GLOCAL_DIR=$GLOCAL_PATH"/"$GPD_PATH"/"$PROC"/Experiments/"$ID
		fi
		if [ $TYPE == "TCE" ];then
			GLOCAL_DIR=$GLOCAL_PATH"/"$GCID_PATH"/"$PROC"/Experiments/"$ID
		fi

	
		mkdir -p $GLOCAL_DIR
		cd $GLOCAL_DIR
	
		DE_FORM_NAME=$ID"_ExpDescription"
		GID="1Gjt_6fcK7XMqOyam_LSoJF1B1-nuhMkst0TE1eMl52g"
		check_and_copy $DE_FORM_NAME $GID
		mkdir "Documents"

	fi

	if [ $TYPE != "DP" ] && [ $TYPE != "TCP" ];then

		mkdir "Reference_files"
		mkdir "Documents"

		### Copy and rename Key_file template to Gdrive dir (if not already there)
		echo "Copying and renaming Key file template"
		KEY_NAME=$ID"_Key"
		KEY_GID="1clmGHbFzKm-fw2EPtLbX52v4MTT8gul_lJ6ehtEIbm0"
		KEY_FLAG=`$GDRIVE_SCRIPT ls | grep $KEY_NAME | wc -l`
		if [ $KEY_FLAG == "0" ];then
			$GDRIVE_SCRIPT copy -id $KEY_GID $KEY_NAME
		fi
		
		### Copy and rename GFF_annotation_key template to Gdrive dir (if not already there)
		
		cd "Reference_files"
		echo "Copying and renaming GFF_key file template"
		GFF_KEY_NAME="gff_ann_key"
		GFF_KEY_GID="1dYkANuEocgXVhMxXjicQzxz0t7ieMwiggEh4d0hcH70"
		KEY_FLAG=`$GDRIVE_SCRIPT ls | grep $GFF_KEY_NAME | wc -l`
		if [ $KEY_FLAG == "0" ];then
			$GDRIVE_SCRIPT copy -id $GFF_KEY_GID $GFF_KEY_NAME
		fi
		cd -
	fi 

	if [ $TYPE == "P" ];then

		### Make summary file
		echo "Making summary file"
		PC_DB=`ls -lrt $PCDB_DIR* | tail -1 | awk '{print $9}'`
		
		ls -lrt $PC_DB
		MOC_FILE=$ID"_summary.tsv"
		
		
		NUM_ID=`cat $PC_DB | grep $ID | wc -l`
		
		if [ $NUM_ID == 1 ];then
			cat $PC_DB | grep -e "Project_Status" -e $ID > $MOC_FILE
		else
			echo $NUM_ID "entries found for "$ID
			exit
		fi
		
		### Make instruction files
		echo "Making instructions file"
		PROJ_INST=$TEMP_DIR"/"$ID"_instructions.txt"
		COL_EMAIL=`FIND_HEADER_AND_FIELD $ID "MOCP_ID" "Collaborator_email" $PC_DB | awk '{print $1}'`	
		FIND_HEADER_AND_FIELD $ID "MOCP_ID" "Collaborator_email" $PC_DB
		echo "Collaborator email:	"$COL_EMAIL > $PROJ_INST
		echo "" >> $PROJ_INST
		cat $FILE_PATH"/collab_instructions.txt" | sed 's/UNIQID/'$ID'/g' >> $PROJ_INST
		$GDRIVE_SCRIPT pub | awk '{print $NF}' >> $PROJ_INST
	
		### Send email with instructions and summary file
		export TMPDIR=$TEMP_PATH
		echo " $PROJ_INST |  mailx -a $MOC_FILE -s "Instructions for MOC project_"$ID $USID"@broadinstitute.org""

		cat $PROJ_INST |  mailx -a $MOC_FILE -s "Instructions for MOC project_"$ID $USID"@broadinstitute.org"

	fi

	
	### Push local dir and contents to Gdrive
	echo "Pushing local files to Gdrive"
	$GDRIVE_SCRIPT push -ignore-name-clashes -quiet
	$GDRIVE_SCRIPT file-id


	### Make files accessible to everyone with link
	ALL_FILES=`ls`

# 	for FILE in $ALL_FILES
# 	do
# 		echo "yes Y | $GDRIVE_SCRIPT share -quiet -with-link -role writer $FILE"
# 		yes Y | $GDRIVE_SCRIPT share -quiet -with-link -role writer $FILE 
# 
# 
# 	done
# # 	echo "yes Y | $GDRIVE_SCRIPT share -quiet -with-link -role writer $KEY_NAME"
# # 	yes Y | $GDRIVE_SCRIPT share -quiet -with-link -role writer $KEY_NAME 


	### change permissions on local files
	change_perms $GLOCAL_DIR $TEMP_DIR 
	ls -lrt $TEMP_DIR

done 

