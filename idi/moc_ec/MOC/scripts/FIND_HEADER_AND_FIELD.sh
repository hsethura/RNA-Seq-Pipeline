#!/bin/sh

Q=$1

###########################    Find values in column with R_HEAD for entries whose value is Q_VAL in column with header Q_HEAD ##################################
FIND_HEADER_AND_FIELD_TEST ()
{
	###########################   FIELD_HEADER  ##################################

	FIELD_HEADER ()
	{

		FILE=$1
		Q_HEAD=$2

		cat $FILE | awk '{if($1 !~ /#/) print $0}' | head -1 | sed 's/ /_/g' | awk -F"\t" -v Q_HEAD=$Q_HEAD '{	
																	for(i=1; i < NF+1; i++)															
																		if($i==Q_HEAD)
																			print i
																}' | tail -1
	}

	
	# Look in FILE, find line in which the value in Q_COL equals Q_VAL and return the values in R_COL 
	
	Q_VAL=$1
	shift
	Q_HEAD=$1
	shift
	R_HEAD=$1
	shift
		
	for FILE in $@
	do

		Q_COL=`FIELD_HEADER $FILE $Q_HEAD` 
		R_COL=`FIELD_HEADER $FILE $R_HEAD` 
#		echo $Q_COL,$R_COL,$FILE
#		
		cat $FILE | awk -F"\t" '{if($1 !~ /#/) print $0}' | sed 's/ /_/g' | awk -F"\t" -v Q_COL=$Q_COL -v R_COL=$R_COL -v Q_VAL=$Q_VAL -v RETURN_COL="-" '{	
																			
																			for(i=1; i < NF+1; i++)
																			{
																									#print $i
																									gsub(/^[_\t]+/, "", $i)
																									
																								
																				#print Q_COL_TRIMMED, Q_VAL, $Q_COL,$i
																				if($i==Q_VAL && R_COL != "")
																					RETURN_COL=$R_COL
																			}
																				if(RETURN_COL != "-")
																					print RETURN_COL

	}'	
	done		
														
} 

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

### determining paths and headers 
### default config file is /idi/moc_ec/MOC/config_files/Universal_config.yaml
paths_and_headers $MOC_ID $@

CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/Universal_config.yaml" 1 $@`
KEY_DIR=`config_read $CONFIG_FILE Key_base`
KEY_ONLY=`extract_option -key_only N 1 $@`
IMPORT_GS=`extract_option -import_gs N 1 $@`

NUM_Q=`echo $Q | sed 's/,/ /g' | wc -w`

Q_VAL=`echo $Q | cut -d"," -f2`
Q_HEAD=`echo $Q | cut -d"," -f1`
 
if [ $NUM_Q == 1 ];then
	Q_HEAD="MOC_ID"
	Q_VAL=$Q
fi

SUB_FILE=`ls -lrt $SUB_DIR/*$SUB_SUFF | tail -1 | awk '{print $9}'`
POOL_FILE=`ls -lrt $POOL_DIR/*$POOL_SUFF | tail -1 | awk '{print $9}'`
DB_FILE=`ls -lrt $MOCDB_DIR/*$MOCDB_SUFF | tail -1 | awk '{print $9}'`
RTS_FILE=`ls -lrt $RTSDB_DIR/*$RTSDB_SUFF | tail -1 | awk '{print $9}'`

ALL_DBS=$SUB_FILE" "$POOL_FILE" "$DB_FILE" "$RTS_FILE

ALL_DB_HEADERS=`cat $RTS_FILE | grep -v "###" | head -1 | sed 's/ /_/g'`
ALL_RETURN_HEADERS="$ALL_DB_HEADERS Index1_Seq Index2_Seq SMOC_ID"


Q_VAL=`echo $Q | cut -d"," -f2`
Q_HEAD=`echo $Q | cut -d"," -f1`

#echo $ALL_DBS, 

# echo $Q_VAL
# echo $Q_HEAD
for R_HEAD in $ALL_RETURN_HEADERS
do
	printf "%s:\t" $R_HEAD | sed 's/_/ /g'		
	#echo "FIND_HEADER_AND_FIELD_TEST $Q_VAL $Q_HEAD $R_HEAD $RTS_FILE"
 	FIND_HEADER_AND_FIELD_TEST $Q_VAL $Q_HEAD $R_HEAD $RTS_FILE | tr ';' ',' | tr ',' '\n' | sort | uniq | tr '\n' ',' | sed 's/^,//g' | sed 's/,$//g' | sed 's/_/ /g'
	echo ""
done

