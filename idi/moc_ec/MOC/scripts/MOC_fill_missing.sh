#!/bin/sh

MOC_ID=$1
shift

### source all functions 
source "/broad/IDP-Dx_storage/MOC/scripts/MOC_functions.sh"

### determining paths and headers 
### default config file is /broad/IDP-Dx_storage/MOC/config_files/Universal_config.yaml
paths_and_headers $MOC_ID $@

### run path_suff function to set RESPATH_SUFF.
##	If -moc_id N included in command line, do not moc_ID to RESPATH_SUFF.  
##	If -user_id included with Y or no string add USID to RESPATH_SUFF.  
##	If -user_id followed by userID, add userID to RESPATH_SUFF.

path_suff $MOC_ID $@
echo $RESPATH_SUFF

USID=`USID`

QUERY_HEADER="Sample"
KEY_HEADER="Sample_ID"

CONFIG_FILE=`extract_option -conf "/broad/IDP-Dx_storage/MOC/config_files/Universal_config.yaml" 1 $@`
USERID_OPT=`extract_option -USER_ID N 1 $@`
MOCID_OPT=`extract_option -MOC_ID Y 1 $@`
KEY_DIR=`config_read $CONFIG_FILE Key_base`
RESULTS_PATH=`config_read $CONFIG_FILE Results_path`
TEMP_PATH=`config_read $CONFIG_FILE Temp_path`
JOIN_PATH=`config_read $CONFIG_FILE join_path`

TEMP_DIR=$TEMP_PATH"/"$RESPATH_SUFF"/"
JOIN_DIR=$JOIN_PATH$MOC_ID"/"

mkdir -p $JOIN_DIR
mkdir -p $TEMP_DIR

TEMP_KEY_FILE=$TEMP_DIR"/"$MOC_ID"_key_temp.txt"

cat $KEY_FILE | sed 's/ /_/g' | grep -v "###" | sed 1d > $TEMP_KEY_FILE

SID_HNAME=`config_read $CONFIG_FILE ID`
BREF_HNAME=`config_read $CONFIG_FILE Ref_accession`
PID_HNAME=`config_read $CONFIG_FILE Proj`
SEQP_HNAME=`config_read $CONFIG_FILE Seq_file`
P7_HNAME=`config_read $CONFIG_FILE P7`
P5_HNAME=`config_read $CONFIG_FILE P5`
BC_HNAME=`config_read $CONFIG_FILE bc`

KEY_F=`FIELD_HEADER $KEY_FILE $SID_HNAME`
KEY_REF_F=`FIELD_HEADER $KEY_FILE $BREF_HNAME`
KEY_PID_F=`FIELD_HEADER $KEY_FILE $PID_HNAME`
PATH_F=`FIELD_HEADER $KEY_FILE $SEQP_HNAME`
INDEX1_F=`FIELD_HEADER $KEY_FILE $P7_HNAME`
INDEX2_F=`FIELD_HEADER $KEY_FILE $P5_HNAME`
BC_F=`FIELD_HEADER $KEY_FILE $BC_HNAME`

typeset -i NL
typeset -i i

NL=`cat $TEMP_KEY_FILE | wc -l`
p=1

PROJID_HEAD=`config_read $CONFIG_FILE Proj`
ALL_PROJIDS=`FIND_VALUES_FOR_HEADER $PROJID_HEAD $KEY_FILE | sort | uniq`
for PROJ_ID in $ALL_PROJIDS
do
	RESULTS_DIR=$RESULTS_PATH"/"$RESPATH_SUFF"/"$KEY_PROJID"/"
	METRICS_FILE=$RESULTS_DIR"/"$KEY_PROJID"_metrics.txt"
	JOIN_FILE=$JOIN_PATH"/"$MOC_ID"/"$PROJ_ID"_KeyMetrics.txt"
	ALL_FIELDS=`FIELD_HEADER $JOIN_FILE Sample_ID Pcnt_bc_in_pool Total_reads pcnt_aligned CDS_total_counts_for_replicon rRNA_pcnt_of_counted | tr '\n' ','`
	TEMP_FILE=$TEMP_DIR$MOC_ID"_"$PROJ_ID"_temp.txt"
	
	echo "Results directory: $RESULTS_DIR"/"$PROJ_ID"/"" > $TEMP_FILE
	echo "" >> $TEMP_FILE
	KEY_SAMPID=`cat $JOIN_FILE | awk -F"\t" -v ALL_FIELDS=$ALL_FIELDS '{
														y=split(ALL_FIELDS,ar,",")
														for(i=1; i < y; i++)
														{
															if(i == 1)
																printf"%-25s\t", $(ar[i])
															else
																printf"%-10s\t", $(ar[i])
														}
														print ""
											
												}' | awk '{
																if(NF != 6) 
																	print $1
																}'`
											
				echo $KEY_SAMPID
# 		cat $KEY_FILE | grep "###"						 			
# 		cat $TEMP_KEY_FILE | awk -F"\t" -v KEY_SAMPID=$KEY_SAMPID '{
# 														
# 															if($'$KEY_F' == KEY_SAMPID)
# 																print $0
# 													
# 													}'
# 
# 		echo $KEY_SAMPID
done



