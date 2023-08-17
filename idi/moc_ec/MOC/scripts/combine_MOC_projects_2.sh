#!/bin/sh


MOC_ID1=$1
MOC_ID2=$2
TAG=$3

##############  FUNCTIONS  #################

	# Pulling out paths, header names from config file 
	config_read ()
	{

		CONFIG_FILE=$1
		VAR_NAME=$2
		cat $CONFIG_FILE | grep $VAR_NAME":" | awk '{print $2}'

	}

	###############  function for extracting options ###############

	# Extracting options from the command line
		extract_option () 
		{
			option_name=$1
			shift
			default=$1
			shift 
			max=$1
			shift 
	
			typeset -i num
			num=`echo $@ | wc -w`
			if [ $num -eq 0 ];then
		
				echo $default
				exit
			
			fi
		
			echo $@ | awk -v def=$default -v name=$option_name -v flag="N" '{
	
				for(i=1; i < NF+1; i++)
				{
					if($i==name)
					{
						flag="Y"
						if(i==NF)
						{
							print "Y"
							exit
						}
						for(n=1; n < (NF+1)-i; n++)
						{
							y=substr($(n+i),1,1)
							##// if next option encountered, print Y and exit
							if(y=="-")
							{
								print "Y"
								exit
							}
							print $(i+n)
						}
					}
			
					if(flag=="N" && i==NF)
						print def
			
				}
	
			}' | tail -$max
		}

###########################################################

### set options
# -conf: sets path to config file (default /broad/IDP-Dx_storage/MOC/config_files/Universal_config.yaml)
# -q: sets name of header for Q_VAL (default MOC_ID)

CONFIG_FILE=`extract_option -conf "/broad/IDP-Dx_storage/MOC/config_files/Universal_config.yaml" 1 $@`
Q_HEAD=`extract_option -q MOC_ID 1 $@`

### set paths to directories, files, and scripts from config file

KEY_DIR=`config_read $CONFIG_FILE Key_base`
RESULTS_PATH=`config_read $CONFIG_FILE Results_path`
TEMP_PATH=`config_read $CONFIG_FILE Temp_path`
SEQ_PATH=`config_read $CONFIG_FILE Seq_base`
INDEX1_BARCODES=`config_read $CONFIG_FILE P7_barcodes`
DB_SCRIPT=`config_read $CONFIG_FILE DB_script`
INDEX1_SPLIT_SCRIPT=`config_read $CONFIG_FILE Index1_split_script`
PIPE_SCRIPT=`config_read $CONFIG_FILE pipeline`
QSUB_PIPE_SCRIPT=`config_read $CONFIG_FILE qsub_pipe`
KEY_SCRIPT=`config_read $CONFIG_FILE keyfile_importer`

### get SMOC and Google IDs associated with $MOC_ID 

ALL_SMOC_ID=`sh $DB_SCRIPT $Q_HEAD,$MOC_ID SMOC_ID | sed 1d | awk '{print $2}' | sed 's/,/ /g' | sed 's/-//g'`
ALL_INDEXES=`sh $DB_SCRIPT $Q_HEAD,$MOC_ID Index1_Seq| sed 1d | awk '{print $2}' | sed 's/,/ /g'`
G_ID=`sh $DB_SCRIPT $Q_HEAD,$MOC_ID Google_ID | sed 1d | awk '{print $2}'`

### set path to key file 

KEY_FILE1=$KEY_DIR$MOC_ID1"_key.txt"
KEY_FILE2=$KEY_DIR$MOC_ID2"_key.txt"
COMB_KEY_FILE=$KEY_DIR$TAG"_key.txt"

ls -lrt $KEY_FILE1
ls -lrt $KEY_FILE2

cat $KEY_FILE1 COMB_KEY_FILE | grep -v "###" > $COMB_KEY_FILE
cat $KEY_FILE2 COMB_KEY_FILE | grep -v "###" | sed 1d >> $COMB_KEY_FILE

ls -lrt $COMB_KEY_FILE


mkdir -p $OUT_DIR

echo $PROJ1_DIR

exit

ALL_KEY="bcLogfile.txt:c+1 counts.tsv:p-1 fpkm.tsv:p-1 corr.txt:c+1 AlignmentSummaryMetrics.txt:c-h metrics.txt:c-h gv.txt:c abs.txt:c info.txt:c+1"

for KEY in $ALL_KEY
do

	SUFF=`echo $KEY | cut -d":" -f1`
	ACT=`echo $KEY | cut -d":" -f2`
	
	FILE1=`ls -lrt $PROJ1_DIR/*$SUFF | awk '{print $9}'`
	FILE2=`ls -lrt $PROJ2_DIR/*$SUFF | awk '{print $9}'`
	OUTFILE=$OUT_DIR"/"$TAG"_"$SUFF


	if [ $ACT == "c" ];then
		cat $FILE1 $FILE2 > $OUTFILE
	fi
	if [ $ACT == "c+1" ];then
		cat $FILE1 > $OUTFILE && echo "" >> $OUTFILE && cat $FILE2 >> $OUTFILE
	fi
	if [ $ACT == "c-h" ];then
		cat $FILE1 > $OUTFILE && cat $FILE2 | sed 1d >> $OUTFILE
	fi
	if [ $ACT == "p-1" ];then
		cat $FILE2 | awk '{for(i=2; i < NF; i++)printf ("%s\t"), $i; print $NF}'> temp.txt
		paste $FILE1 temp.txt > $OUTFILE
	fi
	if [ $ACT == "p" ];then
		paste $FILE1 $FILE2 > $OUTFILE
	fi

done 