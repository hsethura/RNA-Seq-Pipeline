#!/bin/sh

MOC_ID="none"

source /broad/software/scripts/useuse

use Python-2.7

source /idi/moc_ec/MOC/scripts/bash_header

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

fastq_link ()
{
	SMOC_IND=$1
	RAWDATA_PATH=$2
	SYM_DIR=$3
	
	
	echo "working on $SMOC_IND"

	SMOC=`echo $SMOC_IND | cut -d";" -f1`
	I1=`echo $SMOC_IND | cut -d";" -f2`
	I2=`echo $SMOC_IND | cut -d";" -f3`
	ALL_FASTQS=`ls -lrt $RAWDATA_PATH"/"$SMOC"/"* | grep $I1"_"$I2 | awk '{print $NF}'`
	#echo $ALL_FASTQS
	for FASTQ in $ALL_FASTQS
	do
		
		FC=`echo $FASTQ | rev | cut -d"/" -f1 | rev | cut -d"." -f1`
		SYM_FASTQ=`echo $FASTQ | rev | cut -d"/" -f1 | rev | sed 's/'$I1'/'$POOL'/g' | sed 's/'$I2'/'$POOL',/g' | sed 's/'$FC'/'$FC'-'$SMOC'-'$I1'-'$I2'/g' `			
		echo "ln -s $FASTQ $SYM_DIR"/"$SYM_FASTQ"
		ln -s $FASTQ $SYM_DIR"/"$SYM_FASTQ
	done 
}
get_bcs_seq ()
{
	FILE=$1
	ID=$2
	COL=$3
		cat $FILE | awk -v ID=$ID -v COL=$COL '{
										for(i=1;i<NF+1;i++)
											if($i == ID)
												print $(i+COL)
									}' 
}

combine_and_link_fastq ()
{
		RAWDATA_PATH=$1
		shift
		TEMP_DIR=$1
		shift
		SYM_DIR=$1
		shift
		SEG_FILE=$1
		shift
		COMB_FASTQ_FILE=$1
		shift
		COMBINE=$1
		shift
		ALL_REP_SMOC_INDEXES=$@
		
		echo $SYM_DIR
		

		echo "" | sed 1d > $SEG_FILE
		echo "" | sed 1d > $COMB_FASTQ_FILE
		echo "" | sed 1d > $TEMP_DIR/temp3.txt
		
		for SMOC_IND in $ALL_REP_SMOC_INDEXES
		do
			echo $SMOC_IND
			SMOC=`echo $SMOC_IND | cut -d";" -f1`
			I1=`echo $SMOC_IND | cut -d";" -f2`
			I2=`echo $SMOC_IND | cut -d";" -f3`
			ls -l $RAWDATA_PATH"/"$SMOC"/"* | grep $I1"_"$I2 | awk '{print $NF}' > $TEMP_DIR/temp.txt
		
			cat $TEMP_DIR/temp3.txt > $TEMP_DIR/temp2.txt
			paste $TEMP_DIR/temp.txt $TEMP_DIR/temp2.txt > $TEMP_DIR/temp3.txt
		done		
		
		while read line 
		do
			FASTQ=`echo $line | awk '{print $1}'` 
			FC=`echo $FASTQ | rev | cut -d"/" -f1 | rev | cut -d"." -f1`
			DIR_FASTQ=`dirname $FASTQ`
			FASTQ_FILE=`basename $FASTQ`
			COMB_FASTQ_NAME=`echo $FASTQ_FILE | sed 's/'$I1'/'$POOL'/g' | sed 's/'$I2'/'$POOL',/g' | sed 's/'$FC'/'$FC'_cmb/g' `			
			DIR_COMB=$DIR_FASTQ"/combined_fastqs/"
			if [ ! -d $DIR_COMB ];then
				mkdir -p $DIR_COMB 
			fi
			COMB_FASTQ=$DIR_COMB$COMB_FASTQ_NAME
			echo $COMB_FASTQ
			echo line
			echo "Launching UGER job for "$COMB_FASTQ
			echo "cat $line > $COMB_FASTQ" > $TEMP_DIR/temp.txt
			echo $COMB_FASTQ >> $COMB_FASTQ_FILE
			cat $TEMP_DIR/temp.txt
			if [ $COMBINE == "Y" ];then
				qsub $TEMP_DIR/temp.txt >> $SEG_FILE					
			fi
		done < $TEMP_DIR/temp3.txt

		echo "Jobs submitted - waiting for them to complete...."

		echo $SEG_FILE
		qstat

		###testing that all jobs launched are finished
		SGE_test $SEG_FILE	
		echo "All Jobs completed...."
	 
		ls -lrt $COMB_FASTQ_FILE
		while read line
		do	
			FASTQ_NAME=`echo $line | rev | cut -d"/" -f1 | rev`
			SYM_FASTQ=$SYM_DIR"/"$FASTQ_NAME
			echo $SYM_DIR
			echo "ln -s $line $SYM_FASTQ"
			ln -s $line $SYM_FASTQ
		done < $COMB_FASTQ_FILE
}

find_SMOC_IND ()
{
	
	SMOC=$1
	shift
	ALL_SMOCS_INDEXES=$@
	
	echo $ALL_SMOCS_INDEXES | awk -v SMOC=$SMOC '{
															for(i=1; i < NF+1; i++)
															{	
																split($i, ar, ";")
																if(ar[1]==SMOC)
																	print $i
															}
														}'

}


CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/PC_config.yaml" 1 $@`
COMBINE=`extract_option -combine Y 1 $@`
LIMIT_MOCS=`extract_option -limit MOCS 1 $@`
LIMIT=`echo $LIMIT_MOCS | sed 's/\-//g'`
OW_SYM=`extract_option -ow_sym N 1 $@`


### determining paths and headers 
### default config file is /idi/moc_ec/MOC/config_files/Universal_config.yaml
paths_and_headers $MOC_ID $@

mkdir -p MOCSDB_DIR

### set paths to directories, files, and scripts from config file

TEMP_PATH=`config_read $CONFIG_FILE Temp_path`
GLOCAL_PATH=`config_read $CONFIG_FILE gdrivelocal_path`
GDRIVE_SCRIPT=`config_read $CONFIG_FILE gdrive_script`
RAWDATA_PATH=`config_read $CONFIG_FILE Seq_base`
TEMP_DIR=$TEMP_PATH"/MOCS_move/"
mkdir -p $TEMP_DIR


DATA_DIR=`extract_option -data_dir $PARSE_DIR 1 $@`
PARSE=`extract_option -parse Y 1 $@`
IMPORT_GS=`extract_option -import_gs Y 1 $@`

GMOCPM_PATH=`config_read $CONFIG_FILE gdrivepm_path`
GWB_DIR=$GMOCPM_PATH"/Templates_and_workbooks/RtS/Pooling_and_submission_WBs/"
LOCAL_DIR=$GLOCAL_PATH"/"$GWB_DIR

RAW_SYM_PATH=`extract_option -symlink $RAWSYM_PATH 1 $@`
RAW_SEQ_PATH=`extract_option -raw_seq_path $RAWSYM_PATH 1 $@`

echo "RAW_SEQ_PATH: "$RAW_SEQ_PATH
echo "$RAW_SYM_PATH: "$RAW_SYM_PATH

mkdir -p $LOCAL_DIR

echo "cd $GLOCAL_PATH"
cd $GLOCAL_PATH
pwd 


### pull desktop files for all MOCS submission wbs from $GWB_DIR
echo "Pulling all files from "$GWB_DIR
echo "$GDRIVE_SCRIPT pull -no-prompt -files $GWB_DIR"
$GDRIVE_SCRIPT pull -no-prompt -files $GWB_DIR


### pull MOCS submission wbs from $GWB_DIR
cd $LOCAL_DIR
pwd 

#ALL_FILES=`ls -lrt $LOCAL_DIR* | grep $LIMIT_MOCS | grep desk | awk '{print $NF}'`
ALL_FILES=`ls -lrt $LOCAL_DIR* | grep desk | awk '{print $NF}'`


for FILE in $ALL_FILES
do

	GID=`cat $FILE | grep "URL" | cut -d"/" -f6`
	NAME=`cat $FILE | grep "Name=" | cut -d"=" -f2 | sed 's/_Pool_Sub_WB//g'`

	echo "$SCRIPTS_DIR"GS_import.py" -s $GID -t "Pooling" -p $NAME --Key_dir $MOCSDB_DIR -S "_Pool_Sub_WB.txt""
	$SCRIPTS_DIR"GS_import.py" -s $GID -t "Pooling" -p $NAME --Key_dir $MOCSDB_DIR -S "_Pool_Sub_WB.txt" 

done

mkdir -p $BCSDB_DIR
echo "$SCRIPTS_DIR"GS_import.py" -s $BCS_GID -t "NEW barcodes and indexes" -p $BCS_FILE_NAME --Key_dir $BCSDB_DIR"
$SCRIPTS_DIR"GS_import.py" -s $BCS_GID -t "NEW barcodes and indexes" -p $BCS_FILE_NAME --Key_dir $BCSDB_DIR -S $BCS_FILE_SUFF

MOVE_DB=$TEMP_DIR"/move_DB.txt"

## pull all pool IDs from MOCS submission wbs and make move DB
echo mocs_sub_to_db $MOCSDB_DIR "-" $MOCS_DB
mocs_sub_to_db $MOCSDB_DIR "-" $MOCS_DB

ls -lrt $MOCS_DB
cat $MOCS_DB | grep $LIMIT > $MOVE_DB
cat $MOCS_DB | grep $LIMIT

ls -lrt $MOVE_DB
ls -lrt $MOCS_DB

ALL_SYM_DIRS=""

ALL_PROJ=`cat $MOVE_DB | awk '{
															y=split($1, ar, "p")
															for(i=1;i<y-1;i++)
															{
																printf "%sp", ar[i]
															}
															print ar[y-1]

													}' | sort -T /broad/hptmp/ | uniq`


echo "Making symlinks for "$ALL_PROJ 


for PROJ in $ALL_PROJ
do
	SYM_DIR=$RAW_SYM_PATH"/"$PROJ
	if [ -d $SYM_DIR ];then
		if [ $LIMIT == "-" ] || [ $OW_SYM == "Y" ];then
			echo "removing $SYM_DIR"
			rm -r $SYM_DIR
			echo "making $SYM_DIR"
		fi
	fi
	mkdir -p $SYM_DIR
done


ALL_COMBOS=`cat $MOVE_DB`
while read line
do
	POOL=`echo $line | awk '{print $1}'`
	PROJ=`echo $POOL | awk '{
									y=split($1, ar, "p")
									for(i=1;i<y-1;i++)
									{
										printf "%sp", ar[i]
									}
									print ar[y-1]

							}'`
	echo $POOL
	echo $PROJ
	SYM_DIR=$RAW_SYM_PATH"/"$PROJ
	ALL_SYM_DIRS=$ALL_SYM_DIRS" "$SYM_DIR
	
	ALL_SMOCS_INDEXES=`echo $line | awk '{for(i=2;i<NF+1;i++) print $i}'`
	ALL_SMOCS=`echo $line | awk '{
										for(i=2;i<NF+1;i++) 
										{	
											split($i, ar, ";")
											print ar[1]
										}
									}'`


	echo $ALL_SMOCS_INDEXES
	
	for SMOC_IND in $ALL_SMOCS_INDEXES
	do
		echo "	fastq_link $SMOC_IND $RAWDATA_PATH $SYM_DIR"
		fastq_link $SMOC_IND $RAWDATA_PATH $SYM_DIR
	done
# 	NUM_REP_SMOC=`echo $ALL_SMOCS | tr ' ' '\n' | sort | uniq -c | awk '{if($1 > 1) print $2}' | wc -l`
# 	ALL_REP_SMOC=`echo $ALL_SMOCS | tr ' ' '\n' | sort | uniq -c | awk '{if($1 > 1) print $2}'`
# 	ALL_NOREP_SMOC=`echo $ALL_SMOCS | tr ' ' '\n' | sort | uniq -c | awk '{if($1 == 1) print $2}'`
# 	
# 	
# 	echo "ALL_REP_SMOC:" $ALL_REP_SMOC
# 	echo "ALL_NOREP_SMOC:" $ALL_NOREP_SMOC
# 
# 	if [ $NUM_REP_SMOC -gt 0 ];then  
# 		SEG_FILE=$TEMP_DIR$POOL"_SGE_file.txt"
# 		COMB_FASTQ_FILE=$TEMP_DIR$POOL"_combFASTQ_file.txt"
# 
# 		for REP_SMOC in $ALL_REP_SMOC
# 		do
# 			ALL_SMOC_INDEXES=`find_SMOC_IND $REP_SMOC $ALL_SMOCS_INDEXES`
# 			echo $ALL_SMOC_INDEXES
# 			
# 			
# 			### combine fastqs for same pool, same MOCS, diff indexes
# 		
# 			combine_and_link_fastq $RAWDATA_PATH $TEMP_DIR $SYM_DIR $SEG_FILE $COMB_FASTQ_FILE $COMBINE $ALL_SMOC_INDEXES 
# 		done
# 	
# 	fi	
# 			
# 	for NOREP_SMOC in $ALL_NOREP_SMOC
# 	do
# 			
# 		echo $NOREP_SMOC
# 		echo $ALL_SMOCS_INDEXES
# 		ALL_SMOC_INDEXES=`find_SMOC_IND $NOREP_SMOC $ALL_SMOCS_INDEXES`
# 		echo $ALL_SMOC_INDEXES
# 
# 		### make sym link for pools w/o rep indexes in same MOCS
# 		for SMOC_IND in $ALL_SMOC_INDEXES
# 		do
# 			echo "fastq_link $SMOC_IND $RAWDATA_PATH $SYM_DIR"
# 			fastq_link $SMOC_IND $RAWDATA_PATH $SYM_DIR
# 		done
# 	done
	
	
done < $MOVE_DB


echo $MOVE_DB
ls -lrt $MOCS_DB

### change permissions for Results and temp dirs
change_perms $LOCAL_DIR $MOVE_DB $BCSDB_DIR $SYM_DIR $MOCSDB_DIR
