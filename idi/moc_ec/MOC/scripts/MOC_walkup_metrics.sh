#!/bin/sh

SEQ_DIR=$1
shift
SMOC_ID=$1
shift


source /idi/moc_ec/MOC/scripts/bash_header

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

### determining paths and headers 
### default config file is /idi/moc_ec/MOC/config_files/Universal_config.yaml
paths_and_headers $SMOCID $@

### set path to RtS SMOC DB

RTSDB_FILE=`ls -lrt $RTSDB_DIR | tail -1 | awk '{print $9}'`

RTSDB_FILE=$RTSDB_DIR"/"$RTSDB_FILE

echo $RTSDB_FILE

 
############################  stats_parse ###############################

stats_parse ()
{	
	ALL_STATS=$1
	shift 
	TYPE=$1
	shift
	ALL=$@


	all_total=`cat $ALL_STATS |  awk -v total=0 '{total=total+$(NF-1); print total}' | tail -1`

	echo "" | sed 1d > temp_stats_parse.txt
	echo $TYPE"	TOTAL_M_READS	PCNT_OF_TOTAL"
	for m in $ALL
	do	
		m_total=`cat $ALL_STATS | grep $m | awk -v total=0 '{total=total+$NF; print total}' | tail -1`
		cat $ALL_STATS | grep $m | awk -v all_total=$all_total -v total=0 -v m=$m '{total=total+$(NF-1); print $1"\t"total"\t"total/all_total*100}' | tail -1 | sed 's/,//g' >> temp_stats_parse.txt
	done
	cat temp_stats_parse.txt | sort -k3nr

}

########################################################################

CONFIG_FILE=`extract_option -conf "/idi/moc_ec//MOC/config_files/Universal_config.yaml" 1 $SCRIPT_OPTIONS`
TEMP_PATH=`config_read $CONFIG_FILE Temp_path`
DB_SCRIPT=`config_read $CONFIG_FILE DB_script`
echo "Getting all pool ids for $SMOC_ID from $RTSDB_FILE"
ALL_POOL_ID=`sh $DB_SCRIPT SMOC_ID,$SMOC_ID Pool_ID | grep "Pool ID" | awk -F":" '{print $2}'  | sed 's/,/ /g'`

########################################################################

for POOL in $ALL_POOL_ID
do
	echo "sh $DB_SCRIPT Pool_ID,$POOL Pool_Index1_seq -db $RTSDB_FILE "
	I1=`sh $DB_SCRIPT Pool_ID,$POOL Pool_Index1_seq -db $RTSDB_FILE | head -2 | tail -1 | awk -F"\t" '{print $2}' `
	I2=`sh $DB_SCRIPT Pool_ID,$POOL Pool_Index2_seq -db $RTSDB_FILE | head -2 | tail -1 | awk -F"\t" '{print $2}' `
	echo $POOL";"$I1";"$I2
done


USID=`USID`

OUT_DIR=$TEMP_PATH$SMOC_ID"/"
TEMP_DIR=$TEMP_PATH$SMOC_ID"/"$USID"/"

mkdir -p $OUT_DIR
mkdir -p $TEMP_DIR

ALL_STATS=$OUT_DIR"/"$SMOC_ID"_allstats.txt"
METRICS_FILE=$OUT_DIR"/"$SMOC_ID"_walkup_metrics.txt"

SEG_FILE=$TEMP_DIR$SMOC_ID"SGE_file.txt"

echo "" | sed 1d > $SEG_FILE

ALL_FC=`ls -lrt $SEQ_DIR  | grep fastq | grep -e "unmapped" | awk '{print $9}' | cut -d"_" -f1  | cut -d"." -f1 | sort | uniq`
ALL_WU_INDEXES=`ls -lrt $SEQ_DIR  | grep -e "unmapped" -e "unmatched" | awk '{print $9}' | cut -d"." -f3 | sort | uniq`
ALL_LANES=`ls -lrt $SEQ_DIR  | grep -e "unmapped" | awk '{print $9}' | cut -d"." -f2 | sort | uniq`
#ALL_SUFFIX=".1.fastq .2.fastq"
ALL_SUFFIX=".1.fastq"

NUM_INDEX=`echo $ALL_WU_INDEXES | wc -w`

echo "SEQ_DIR:" $SEQ_DIR
echo "ALL_WU_INDEXES:" $ALL_WU_INDEXES 
echo "NUM_INDEX:"	$NUM_INDEX 
echo "ALL_FC :"	$ALL_FC 
echo "ALL_LANES :"	$ALL_LANES 



for FC in $ALL_FC
do
	echo "Total files for FC "$FC"..."
	for LANE in $ALL_LANES
	do 
		printf "LANE %s:\t" $LANE
		ls -lrt $SEQ_DIR"/"$FC*$LANE* | wc -l
	done
	
done > $METRICS_FILE	

 
for FC in $ALL_FC
do	
	echo "Total files for FC "$FC"..."
	for INDEX in $ALL_WU_INDEXES
	do
		printf "INDEX %s:\t" $INDEX
		ls -lrt $SEQ_DIR"/"*$FC*$INDEX* | wc -l
	done
	
done >> $METRICS_FILE	

typeset -i NUM_FILES

echo "" | sed 1d > ~/temp.txt

# Removing old counts file for each fastq

rm -f $OUT_DIR"/"*"_counts.txt"

# Getting counts for each fastq

for FC in $ALL_FC
do
	for INDEX in $ALL_WU_INDEXES
	do
		FLAG=""
		for SMOC_INDEX in $ALL_WU_INDEXES
		do
			echo $SMOC_INDEX
			
			if [ $INDEX == $SMOC_INDEX ];then
				FLAG="*"
			fi
		done
		
		for LANE in $ALL_LANES
		do
			for SUFF in $ALL_SUFFIX
			do
				NUM_FILES=`ls -lrt $SEQ_DIR | grep -v barcode | grep -v md5 | grep fastq | grep $SUFF | awk -v INDEX=$INDEX -v FC=$FC -v SUFF=$SUFF -v LANE=$LANE '{
														
														split($9, ar, ".")
														if(ar['$FC_F']==FC && ar['$INDEX_F']==INDEX && ar['$LANE_F']==LANE)
															print $NF, FC, INDEX, LANE

														}' | wc -l`
										
				if [ $NUM_FILES == 1 ];then
				
					FASTQ=`ls -lrt $SEQ_DIR | grep -v barcode | grep -v md5 | grep fastq | grep $SUFF | awk -v INDEX=$INDEX -v FC=$FC -v SUFF=$SUFF -v LANE=$LANE '{
														
														split($9, ar, ".")
														if(ar['$FC_F']==FC && ar['$INDEX_F']==INDEX && ar['$LANE_F']==LANE)
															print $NF

														}'`
					ROOT=`basename $FASTQ .fastq.gz`
					COUNT_FILE=$OUT_DIR"/"$ROOT"_counts.txt"
					echo "Launching UGER job for "$FASTQ
					
					echo "echo $FC, $INDEX, $LANE, $SUFF, $FASTQ | tr '\n' ' ' > $COUNT_FILE" >~/temp.txt
					echo "zcat $FASTQ | grep @ | wc -l >> $COUNT_FILE" >> ~/temp.txt
					
					qsub ~/temp.txt >> $SEG_FILE					
				
				fi
			done
		done
	done
done

echo "Jobs submitted - waiting for them to complete...."

echo $SEG_FILE
qstat

# testing that all jobs launched are finished
SGE_test $SEG_FILE	
echo "All Jobs completed...."

total=`cat $OUT_DIR"/"*"_counts.txt" | awk -v total=0 '{total=total+$NF; print total}' | tail -1`
echo $total 
echo "INDEX	LANE	FILE	TOTAL_M_READS (read 1 only) PCNT_OF_TOTAL_READS" > $ALL_STATS
cat $OUT_DIR"/"*"_counts.txt" | awk -v total=$total '{print $2, $3, $5, $NF/1000000, $NF/total*100}' >> $ALL_STATS

printf "%s	%s	%s	%s	" "Flowcell" "TotalMReads"	"MReadsMatchedToIndex:(pcnt)" "MReadsNotMatchedToIndex:(pcnt)" >> $METRICS_FILE	
for LANE in $ALL_LANES
do
	printf "	Lane_%s" $LANE >> $METRICS_FILE	

done 

echo ""  >> $METRICS_FILE	


for FC in $ALL_FC
do
	
	TOTAL=`cat $ALL_STATS | grep $FC | awk -v total_reads=0 -v total_pcnt=0 -v i=0 '{total_reads=total_reads+$4;total_pcnt=total_pcnt+$5;i++; print total_reads}' | tail -1`
	MATCHED=`cat $ALL_STATS | grep $FC | grep -v "unmatched" | awk -v total_reads=0 -v TOTAL=$TOTAL -v i=0 '{total_reads=total_reads+$4;total_pcnt=total_pcnt+$5;i++; print total_reads":("total_reads/TOTAL*100")"}' | tail -1`
	UNMATCHED=`cat $ALL_STATS | grep $FC | grep "unmatched" | awk -v total_reads=0 -v TOTAL=$TOTAL -v i=0 '{total_reads=total_reads+$4;total_pcnt=total_pcnt+$5;i++; print total_reads":("total_reads/TOTAL*100")"}' | tail -1`
	
	printf "%s	%s	%s	%s	" $FC $TOTAL $MATCHED $UNMATCHED
	for LANE in $ALL_LANES
	do
		TOTAL_LANE=`cat $ALL_STATS | grep $FC"."$LANE | awk -v TOTAL=$TOTAL -v total_pcnt=0 -v i=0 '{total_reads=total_reads+$4; print total_reads/TOTAL*100}' | tail -1`
		printf "%s\t" $TOTAL_LANE
	done
	echo ""


done >> $METRICS_FILE

echo "" >> $METRICS_FILE

echo $ALL_FC
echo $METRICS_FILE

stats_parse $ALL_STATS "Index" $ALL_WU_INDEXES >> $METRICS_FILE

# echo "* indexes associated with SMOC in MOC DBs" >> $METRICS_FILE
# echo "" >> $METRICS_FILE


echo "Top_10_nmatched_Index	%_of_total" >> $METRICS_FILE
zcat $SEQ_DIR/*unmatched*barcode*fastq.gz | head -100000 | awk '{if(NR % 4 == 2)print $1}' | sort | uniq -c | sort -k1nr | head -10| awk '{print $2, $1/100000*100}' >> $METRICS_FILE

cp $METRICS_FILE $SEQ_DIR

chmod 777 $SEQ_DIR/*_walkup_metrics*

export TMPDIR=$TEMP_PATH

cat $METRICS_FILE |  mailx -a $METRICS_FILE -s "Metrics for $SMOC_ID" $USID"@broadinstitute.org"

change_perms $METRICS_FILE 