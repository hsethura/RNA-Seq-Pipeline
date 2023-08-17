#!/bin/sh

source /idi/moc_ec/MOC/scripts/bash_header

MOC_ID=$1

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

CONFIG_FILE=`extract_option -conf "/idi/moc_ec/MOC/config_files/Universal_config.yaml" 1 $@`
ALL_PAIRS_OPT=`extract_option -all_pairs "-" 1 $@`
RUN_DESEQ=`extract_option -run_deseq Y 1 $@`
RUN_EDGER=`extract_option -run_edger Y 1 $@`

# Get paths to dirs scripts from config file
read_config $CONFIG_FILE 

### determining paths and headers 
### default config file is /idi/moc_ec/MOC/config_files/Universal_config.yaml
paths_and_headers $MOC_ID $@

### run path_suff function to set RESPATH_SUFF.
##	If -moc_id N included in command line, do not moc_ID to RESPATH_SUFF.  
##	If -user_id included with Y or no string add USID to RESPATH_SUFF.  
##	If -user_id followed by userID, add userID to RESPATH_SUFF.
path_suff $@

KEY_FILE=$KEY_DIR$MOC_ID"_key.txt"
RESULTS_DIR=$RESULTS_PATH"/"$RESPATH_SUFF"/"
TEMP_DIR=$TEMP_PATH"/"$RESPATH_SUFF"/"

echo $RESPATH_SUFF

mkdir -p $TEMP_DIR
TEMP_FILE=$TEMP_DIR"/temp_"$MOC_ID"_DSlaunch.txt"

# typeset -i SAMP_ID_FIELD
# typeset -i CG_ID_FIELD

EXCLUDE=`extract_option -EXC - 100 $@`
INCLUDE=`extract_option -INC - 100 $@`

EDGR_PATH=`config_read $CONFIG_FILE DE_path`
RESULTS_PATH=`config_read $CONFIG_FILE Results_path`

SAMP_ID_FIELD=`FIELD_HEADER $KEY_FILE Sample_ID`
CG_ID_FIELD=`FIELD_HEADER $KEY_FILE CG_ID`
PID_FIELD=`FIELD_HEADER $KEY_FILE Project_ID`
ACC_FIELD=`FIELD_HEADER $KEY_FILE Bacterial_reference`
PAIRS_FIELD=`FIELD_HEADER $KEY_FILE CG_Pairs`

echo $CG_ID_FIELD

### get GID from MOC_DB or use GID included as command line option
# if -GID_OPT not included in command line, get it from DB

if [ $GID_OPT == "0" ];then
	G_ID=`sh $DB_SCRIPT $Q_HEAD,$MOC_ID Google_ID | grep "Google ID" | awk -F":" '{print $2}' | sed 's/ /_/g'`
else
	G_ID=$GID_OPT
fi

echo $G_ID

if [ -z $G_ID ];then

	echo "GID not found in DB or entered in command line with -gid option"
	exit
fi

echo $MOVE_KEY

############## moving key file to server ###############
if [ $MOVE_KEY == "Y" ];then 
	echo "Moving key file to server..."
	echo "$KEY_SCRIPT -s $G_ID -t \"$KEY_SHEET\" -p $MOC_ID --Key_dir $KEY_DIR"
	
	
	$KEY_SCRIPT -s $G_ID -t "$KEY_SHEET" -p $MOC_ID --Key_dir $KEY_DIR
	### if key file not found or empty, stop pipeline
	if [ ! -s $KEY_FILE ];then
		ls -lrt $KEY_FILE
		exit
	fi
fi
######################################################## 
 

cat $KEY_FILE | sed 's/, /,/g' | sed 's/ /_/g' |  sed 's*/*_*g' | sed 's/	_/	/g' > $TEMP_FILE

ALL_PROJID=`cat $TEMP_FILE | awk -F"\t" '{if(NR>5) print $0}' | awk -F"\t" -v PID_FIELD=$PID_FIELD '{print $PID_FIELD}' | sort | uniq`

echo $KEY_FILE
echo $ALL_PROJID

SEG_FILE=$TEMP_DIR$MOC_ID"_DE_SGE_file.txt"

echo "" | sed 1d > $SEG_FILE


for PROJ_ID in $ALL_PROJID
do

	if [ $ALL_PAIRS_OPT == "-" ];then
	 	ALL_PAIRS=`cat $TEMP_FILE | awk -F"\t" '{if(NR>5) print $0}' | awk -F"\t" -v PAIRS_FIELD=$PAIRS_FIELD -v PID_FIELD=$PID_FIELD -v PROJ_ID=$PROJ_ID '{if($PID_FIELD==PROJ_ID && $PAIRS_FIELD !~ /#/ ) print $PAIRS_FIELD}' | sed '/^$/d' | sed '/^$/d' | sort | uniq`
	else
		ALL_PAIRS=`echo $ALL_PAIRS_OPT | sed 's/:/ /g'`
	fi
	
	NUM_PAIRS=`echo $ALL_PAIRS | wc -w`
	echo "ALL_PAIRS:" $ALL_PAIRS
	echo "NUM_PAIRS:" $NUM_PAIRS

	for PAIR in $ALL_PAIRS
	do
		NUM_PAIR=`echo $PAIR | awk '{y=split($1, ar, ","); print y}'`
		if [ $NUM_PAIR -ne 2 ];then
			echo $PAIR" does not have the proper number of CG_IDs!"
			exit
		fi
	
	done
	
	if [ $PAIRS_FIELD == "0" ];then

		echo "No CG_Pairs field found..."
		exit
	fi


	MISSING_EDGER=""		
	MISSING_DESEQ=""	
	MISSING_EDGER_NUM=0	
	MISSING_DESEQ_NUM=0	

	for PAIR in $ALL_PAIRS
	do
		
		echo $PAIR

		GCIDA=`echo $PAIR | cut -d"," -f1 | sed 's/:/,/g'`
		GCIDB=`echo $PAIR | cut -d"," -f2 | sed 's/:/,/g'`

		A=`DE_FIND_FIELD $TEMP_FILE $PID_FIELD $PROJ_ID $GCIDA $SAMP_ID_FIELD $CG_ID_FIELD $SAMP_ID_FIELD | sed 's/,$//g' | tr ',' '\n' | sed '/^$/d' | sort -u | tr '\n' ','` 
		B=`DE_FIND_FIELD $TEMP_FILE $PID_FIELD $PROJ_ID $GCIDB $SAMP_ID_FIELD $CG_ID_FIELD $SAMP_ID_FIELD | sed 's/,$//g' | tr ',' '\n' | sed '/^$/d' | sort -u | tr '\n' ','` 
		ACC=`DE_FIND_FIELD $TEMP_FILE $PID_FIELD $PROJ_ID $GCIDA $ACC_FIELD $CG_ID_FIELD $CG_ID_FIELD $SAMP_ID_FIELD | sed 's/,,/,/g' | tr ',' '\n' | sed '/^$/d'  | sort -u` 
					
		echo $TEMP_FILE
				
		PROJ_DIR=$RESULTS_DIR"/"$PROJ_ID"/"
		DE_DIR=$PROJ_DIR"/DE/"
		EDGR_DIR=$DE_DIR"/EdgeR/" 
		DESEQ_DIR=$DE_DIR"/DESeq2/"
		READ_FILE=$PROJ_DIR"/"$PROJ_ID"_"$ACC"_counts.tsv"
		ALL_READ_FILE=$PROJ_DIR"/"$PROJ_ID"_"$ACC"_all_counts.tsv"
		GENE_READ_FILE=$PROJ_DIR"/"$PROJ_ID"_"$ACC"_genes_counts.tsv"
		EDGR_QSUB_FILE=$TEMP_DIR"/EdgeR_qsub.txt"
		DESEQ_QSUB_FILE=$TEMP_DIR"/DESeq_qsub.txt"

		
		cat $READ_FILE | tr -d \'\"  > $ALL_READ_FILE
		cat $READ_FILE | tr -d \'\" | awk '{if($1 ~ /^CDS:/ || $1 ~ /^ncRNA:/ || NR==1) print $0}' > $GENE_READ_FILE

		if [ $RUN_DESEQ == "Y" ] || [ $RUN_EDGER == "Y" ];then
					
			#echo "cat $READ_FILE | tr -d \'\" > $EDGER_READ_FILE"
			
			ls -lrt $GENE_READ_FILE
			ls -lrt  $READ_FILE			

			echo $TEMP_FILE
	
			mkdir -p $EDGR_DIR
			mkdir -p $DESEQ_DIR
	
			echo "#########################"
			echo $PID_FIELD
			echo $PID_FIELD
			echo "A" $A 
			echo "B" $B 
			echo "GCIDA  "$GCIDA
			echo "GCIDB  "$GCIDB
			echo "ACC  "$ACC
			echo "READ_FILE	"$READ_FILE	
			echo "ALL_READ_FILE "$ALL_READ_FILE
			echo "GENE_READ_FILE "$GENE_READ_FILE
			echo "EDGR_DIR "$EDGR_DIR
			echo "#########################"
		fi

		
		if [ $RUN_DESEQ == "Y" ];then
			if [ ! -s $READ_FILE ];then
				echo $READ_FILE "missing!"
			fi
			
			echo "source /idi/moc_ec/MOC/scripts/bash_header" > $DESEQ_QSUB_FILE
			echo "sh $DESEQ_SCRIPT $MOC_ID $PROJ_ID $ALL_READ_FILE $DESEQ_DIR $ACC "All" $A $B $GCIDA $GCIDB -test LRT $@" >> $DESEQ_QSUB_FILE
			echo "sh $DESEQ_SCRIPT $MOC_ID $PROJ_ID $GENE_READ_FILE $DESEQ_DIR $ACC "Genes" $A $B $GCIDA $GCIDB -test LRT $@" >> $DESEQ_QSUB_FILE

			echo $DESEQ_QSUB_FILE
		
			qsub -l h_rt=24:00:00 -l h_vmem=8g -l os=RedHat7 $DESEQ_QSUB_FILE >> $SEG_FILE	
			echo  "qsub -l h_rt=24:00:00 -l h_vmem=8g -l os=RedHat7 $DESEQ_QSUB_FILE" 
			change_perms $DESEQ_DIR 	
		fi
	
		
		if [ $RUN_EDGER == "Y" ];then
			if [ ! -s $EDGER_READ_FILE ];then
				echo $EDGER_READ_FILE "missing!"
			fi
			ER_A=`echo $A | sed 's/,/;/g'`
			ER_B=`echo $B | sed 's/,/;/g'`
			
			echo "source /idi/moc_ec/MOC/scripts/bash_header" > $EDGR_QSUB_FILE
			echo "$EDGER_SCRIPT -i $ALL_READ_FILE -c \"$ER_A\" -n \"$ER_B\" -t lr -o $EDGR_DIR --tag $GCIDA"_"$GCIDB"_All"" >> $EDGR_QSUB_FILE
			echo "$EDGER_SCRIPT -i $EDGER_READ_FILE -c \"$ER_A\" -n \"$ER_B\" -t lr -o $EDGR_DIR --tag $GCIDA"_"$GCIDB"_Genes"" >> $EDGR_QSUB_FILE
			
			qsub -l h_rt=24:00:00 -l h_vmem=8g -l os=RedHat7 $EDGR_QSUB_FILE >> $SEG_FILE	
			echo  "qsub -l h_rt=24:00:00 -l h_vmem=8g -l os=RedHat7 $EDGR_QSUB_FILE"

			change_perms $EDGR_DIR  	
		fi

	done
	
	echo "Jobs submitted - waiting for them to complete...."
	echo $SEG_FILE

	# testing that all jobs launched are finished
	SGE_test $SEG_FILE	
	qstat 

	for PAIR in $ALL_PAIRS
	do
		
		echo $PAIR

		GCIDA=`echo $PAIR | cut -d"," -f1 | sed 's/:/,/g'`
		GCIDB=`echo $PAIR | cut -d"," -f2 | sed 's/:/,/g'`
		echo "Calculating number of files for "$GCIDA","$GCIDB
		NUM_FILES=`ls -lrt $DESEQ_DIR | grep $GCIDA | grep $GCIDB | grep "All.txt" | wc -l`		
		if [ $NUM_FILES == 0 ];then
			echo $GCIDA","$GCIDB "pair is missing from DESeq2 comparisons"
			MISSING_DESEQ=$GCIDA","$GCIDB":"$MISSING_DESEQ
			MISSING_DESEQ_NUM=`expr $MISSING_DESEQ_NUM + 1`

		fi				
		
		NUM_FILES=`ls -lrt $EDGR_DIR | grep $GCIDA | grep $GCIDB | grep "_lr_comp.tsv" | wc -l`
		if [ $NUM_FILES == 0 ];then
			echo $GCIDA","$GCIDB "pair is missing from EdgeR comparisons"
			MISSING_EDGER=$GCIDA","$GCIDB":"$MISSING_EDGER
			MISSING_EDGER_NUM=`expr $MISSING_EDGER_NUM + 1`
		fi				
	done
	
	echo $MISSING_DESEQ_NUM " DESeq2 pairs are missing!"
	echo $MISSING_EDGER_NUM " edgeR pairs are missing!"
	USER_ID=`echo $@ | awk '{
								for(i=1; i < NF+1; i++)
								{
									if($i == "-user_id")
										print $i
								}
							}'`
								
	echo sh $DE_PIPE $MOC_ID -run_deseq N -all_pairs $MISSING_EDGER $@
	echo sh $DE_PIPE $MOC_ID -run_edger N -all_pairs $MISSING_DESEQ $@
done

