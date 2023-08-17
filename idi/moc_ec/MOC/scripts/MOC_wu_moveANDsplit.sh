#!/bin/sh

SMOCID=$1 

source /idi/moc_ec/MOC/scripts/bash_header

### source all functions 
source "/broad/IDP-Dx_storage/MOC/scripts/MOC_functions.sh"

### determining paths and headers 
### default config file is /broad/IDP-Dx_storage/MOC/config_files/Universal_config.yaml
paths_and_headers $SMOCID $@

echo $USID

MOVE_KEY=`extract_option -move_key Y 1 $SCRIPT_OPTIONS`
IMPORT_GS=`extract_option -import_gs N 1 $SCRIPT_OPTIONS`
RAW_SEQ_PATH=`extract_option -raw_seq_path $SEQ_PATH 1 $@`
CHECKFILES=`extract_option -checkfiles N 1 $@`
MOC_ID=`extract_option -moc_id - 1 $@`

SEQ_PATH=$RAW_SEQ_PATH
echo $SEQ_PATH

SMOCDIR_ID=`echo $SMOCID | sed 's/-//g'`

SEQ_DIR=$SEQ_PATH"/"$SMOCDIR_ID"/"
echo "Making "$SEQ_DIR
mkdir -p $SEQ_DIR
chmod 777 $SEQ_DIR
cd $SEQ_DIR

WGET_DIR=$SEQ_DIR"/get.broadinstitute.org/pkgs/"$USER"/"
VAL_FILE=$SEQ_DIR"/Download_validation.txt"

if [ $USER != "N" ] && [ $PSSWD != "N" ];then
	MOVE_DATA="Y" 
fi	

echo $MOVE_DATA, $USER, $PSSWD

if [ $MOVE_DATA == "Y" ];then
	if [ $USER == "N" ] || [ $PSSWD == "N" ];then
		echo "Missing userID and/or psswd"
		exit 
	fi
	
	echo "$WGET_DIR"
	mkdir -p $WGET_DIR
	
	echo "wget --tries=10 --continue --mirror --user $USER --password $PSSWD --no-check-certificate https://get.broadinstitute.org/pkgs/$USER/"

	wget --tries=10 --continue --mirror --user $USER --password $PSSWD --no-check-certificate https://get.broadinstitute.org/pkgs/$USER/
	
	MANIFEST_FILE=$WGET_DIR"MANIFEST"
	NEW_MANIFEST_NAME=$SEQ_DIR$USER"_MANIFEST"
	
	cp $MANIFEST_FILE $NEW_MANIFEST_NAME

	ALL_FASTQS=`ls -lrt $WGET_DIR*fastq* | awk '{print $9}'`
	
	for FASTQ in $ALL_FASTQS
	do
			FILE=`echo $FASTQ | rev | cut -d"/" -f1 | rev`
			ln -s $FASTQ $SEQ_DIR$FILE 2>/dev/null
	done
 
fi

if [ $CHECKFILES == "Y" ] || [ $MOVE_DATA == "Y" ];then
	######### compare file list to manifest ####################

	NOTFOUND="N"
	echo "" | sed 1d > $VAL_FILE
	
	echo "Checking all files were downloaded"
	NEW_MANIFEST_NAME=$SEQ_DIR$USER"_MANIFEST"
	ALL_FASTQS=`cat $NEW_MANIFEST_NAME | grep fastq | awk '{print $1}'`
	
	for FASTQ in $ALL_FASTQS
	do
		    DL_FASTQ=`ls -lrt $SEQ_DIR/""$FASTQ | awk '{print $NF}'`
		    DL_size=`ls -lrt $DL_FASTQ | awk '{printf "%0.2fG",  $5/1e9}' `
		    if [ ! -s $DL_FASTQ ];then
		    	NOTFOUND="Y"
		    	echo $FASTQ" not found in $SEQ_DIR!" >> $VAL_FILE
		    else
		    	echo $FASTQ" "$DL_size >> $VAL_FILE
		    fi
	done
	cat $VAL_FILE | sort -k2nr

	export TMPDIR=$TEMP_PATH

	if [ $NOTFOUND == "Y" ];then
		cat $VAL_FILE | sort -k2nr |  mailx -s "Download NOT complete for $SMOCID" $USID"@broadinstitute.org"
	else
		cat $VAL_FILE | sort -k2nr |  mailx -s "Download COMPLETE for $SMOCID" $USID"@broadinstitute.org"
	fi	
fi

############## Run metrics script ###############=

if [ $METRICS == "Y" ];then
	echo "Running metrics script"
	echo "sh $WU_SCRIPT $SEQ_DIR $SMOCID"
	sh $WU_SCRIPT $SEQ_DIR $SMOCID
fi
########################################################

############## Import google sheet databases to server ###############=
if [ $IMPORT_GS == "Y" ];then
	echo "Running $GSIMPORT_SCRIPT to import google sheet databases to server"
	sh $GSIMPORT_SCRIPT
fi
########################################################

################ split fastq file in SEQ_DIR by all indexes and write the split files into SEQ_DIR ################################

if [ $INDEX1_SPLIT == "Y" ];then

	INDEX_SPLIT_INDIR=$TEMP_PATH"/"$SMOCID"/symlinks/"
	INDEX_SPLIT_OUTDIR=$INDEX_SPLIT_PATH"/"$SMOCID"/"
	mkdir -p $INDEX_SPLIT_INDIR
	mkdir -p $INDEX_SPLIT_OUTDIR

	############## determining ALL_SMOCINDEX_SETS for each MOC_ID ###############
	if [ $MOC_ID == "-" ];then
		sh $DB_SCRIPT SMOC_ID,$SMOCID MOC_ID
		ALL_MOCIDS=`sh $DB_SCRIPT SMOC_ID,$SMOCID MOC_ID | sed 1d | grep -v SMOC | sed 's/,/ /g' | awk -F":" '{for(i=2; i < NF+1; i++)print $i}'`
	else 
		ALL_MOCIDS=$MOC_ID
	fi
	
	if [ -z $ALL_MOCIDS ];then
		echo "ALL_MOCIDS empty!:" 
		exit
	fi
		
	for MOC_ID in $ALL_MOCIDS
	do
		### get GID from MOC_DB or use GID included as command line option
		# if -GID_OPT not included in command line, get it from DB

		if [ $GID_OPT == "0" ];then
		G_ID=`sh $DB_SCRIPT $Q_HEAD,$MOC_ID Google_ID | grep "Google ID" | sed 's/ /_/g' | awk '{print $2}'`
		else
			G_ID=$GID_OPT
		fi
		if [ -z $G_ID ];then

			echo "GID not found in DB or entered in command line with -gid option"
			exit
		fi

		############## moving key file to server ###############
		if [ $MOVE_KEY == "Y" ];then 
			echo "Moving key file to server..."
			echo "$KEY_SCRIPT -s $G_ID -t \"Sample Information\" -p $MOC_ID --Key_dir $KEY_DIR"
			$KEY_SCRIPT -s $G_ID -t "Sample Information" -p $MOC_ID --Key_dir $KEY_DIR
			KEY_FILE=$KEY_DIR$MOC_ID"_key.txt"
			### if key file not found or empty, stop pipeline
			if [ ! -s $KEY_FILE ];then
				ls -lrt $KEY_FILE
				exit
			fi
		fi
		########################################################
		
		
		echo "Determining ALL_SMOCINDEX_SETS for $MOC_ID..."
		all_smocindex_set $MOC_ID $@
		if [ -z $ALL_SMOCINDEX_SETS ];then
			echo "ALL_SMOCINDEX_SETS empty!:" 
			exit
		else
			echo "ALL_SMOCINDEX_SETS:" $ALL_SMOCINDEX_SETS
		fi
	done
	########################################################
	
 	SMOC_P7_BARCODES=$INDEX_SPLIT_INDIR"/"$SMOCID"_P7.txt"
	SMOC_P5_BARCODES=$INDEX_SPLIT_INDIR"/"$SMOCID"_P5.txt"
	SMOC_P7_LIB=$INDEX_SPLIT_INDIR"/"$SMOCID"_P7_lib.dat"


	echo "" | sed 1d > $SMOC_P7_BARCODES
	echo "" | sed 1d > $SMOC_P5_BARCODES

	
	for SMOCINDEX_SET in $ALL_SMOCINDEX_SETS
	do
		SMOCID_SET=`echo $SMOCINDEX_SET | cut -d";" -f1`

		if [ $SMOCID_SET == $SMOCID ];then
			INDEX1=`echo $SMOCINDEX_SET | cut -d";" -f2`
			INDEX2=`echo $SMOCINDEX_SET | cut -d";" -f3`

			cat $P7_BARCODES | grep $INDEX1 >> $SMOC_P7_BARCODES
			cat $P5_BARCODES | grep $INDEX2 >> $SMOC_P5_BARCODES
		fi
	done
	
	cat $SMOC_P7_BARCODES | sort -r | uniq > $INDEX_SPLIT_INDIR"/temp.txt" && echo "i7_name	i7_seq" > $SMOC_P7_BARCODES && cat $INDEX_SPLIT_INDIR"/temp.txt" >> $SMOC_P7_BARCODES
	cat $SMOC_P5_BARCODES | sort -r | uniq > $INDEX_SPLIT_INDIR"/temp.txt" && echo "i5_name	i5_seq" > $SMOC_P5_BARCODES && cat $INDEX_SPLIT_INDIR"/temp.txt" >> $SMOC_P5_BARCODES
	
	cat $SMOC_P7_BARCODES
	cat $SMOC_P5_BARCODES


	ALL_FASTQS=`ls -lrt $SEQ_DIR*.fastq* | grep unmatched | awk '{print $9}'`
	
	echo $ALL_FASTQS
	
	for FASTQ in $ALL_FASTQS
	do
			FILE=`echo $FASTQ | rev | cut -d"/" -f1 | rev | sed 's/unmatched.//g'`
			ln -s $FASTQ $INDEX_SPLIT_INDIR$FILE 2>/dev/null
	done

	### single split
	if [ $INDEX2_SPLIT == "N" ];then
		echo "Splitting un-demultiplexed files in $INDEX_SPLIT_INDIR by P7 indexes and outputting split files to $INDEX_SPLIT_OUTDIR..."
		echo "Building dictionary"
		$DICT_BUILD -i $P7_BARCODES -o $SMOC_P7_LIB
				
		if [ $QSUB == "Y" ];then
			echo "$QSUB_SCRIPT $SINGLE_INDEX_SCRIPT -i $INDEX_SPLIT_INDIR -o $INDEX_SPLIT_OUTDIR -d $SMOC_P7_LIB"
		 	sh $QSUB_SCRIPT $SINGLE_INDEX_SCRIPT -i $INDEX_SPLIT_INDIR -o $INDEX_SPLIT_OUTDIR -d $SMOC_P7_LIB
		else
			echo "$SINGLE_INDEX_SCRIPT -i $INDEX_SPLIT_INDIR -o $INDEX_SPLIT_OUTDIR -d $SMOC_P7_LIB"
			$SINGLE_INDEX_SCRIPT -i $INDEX_SPLIT_INDIR -o $INDEX_SPLIT_OUTDIR -d $SMOC_P7_LIB

		fi
	fi
	
	### dual split
	if [ $INDEX2_SPLIT == "Y" ];then
		echo "Splitting un-demultiplexed files in $INDEX_SPLIT_INDIR by P7 and P5 indexes and outputting split files to $INDEX_SPLIT_OUTDIR..."
		if [ $QSUB == "Y" ];then
		 	sh $QSUB_SCRIPT $DUAL_INDEX_SCRIPT -i $INDEX_SPLIT_INDIR -o $INDEX_SPLIT_OUTDIR --dict1_file $SMOC_P7_BARCODES --dict2_file $SMOC_P5_BARCODES
		else
			echo "$DUAL_INDEX_SCRIPT -i $INDEX_SPLIT_INDIR -o $INDEX_SPLIT_OUTDIR --dict1_file $SMOC_P7_BARCODES --dict2_file $SMOC_P5_BARCODES"
			$DUAL_INDEX_SCRIPT -i $INDEX_SPLIT_INDIR -o $INDEX_SPLIT_OUTDIR --dict1_file $SMOC_P7_BARCODES --dict2_file $SMOC_P5_BARCODES
		fi
	fi

	ALL_SPLIT_FASTQS=`ls -lrt $INDEX_SPLIT_OUTDIR | grep fastq | awk '{print $9}'`
	
	# make symlinks in SEQDIR to index split fastqs in $INDEX_SPLIT_OUTDIR
	for FASTQ in $ALL_FASTQS
	do
			FILE=`echo $FASTQ | rev | cut -d"/" -f1 | rev | sed 's/unmatched.//g'`
			ln -s $FASTQ $SEQ_DIR | grep fastq.gz 2>/dev/null
	done


fi

### change permissions for Results and temp dirs
change_perms $INDEX_SPLIT_OUTDIR $SEQ_DIR $WGET_DIR
