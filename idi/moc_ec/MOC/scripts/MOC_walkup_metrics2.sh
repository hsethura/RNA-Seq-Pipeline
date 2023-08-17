#!/bin/sh

SEQ_DIR=$1
MOCS_ID=$2

source /idi/moc_ec/MOC/scripts/bash_header

### source all functions 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"

### determining paths and headers 
### default config file is /idi/moc_ec/MOC/config_files/PC_config.yaml
paths_and_headers $SMOCID $@


echo $IMPORT_GS

############################  stats_parse ###############################

stats_parse ()
{	
	
	ALL_STATS=$1
	shift 
	TYPE=$1
	shift
	ALL_POOLS=$1
	shift
	ALL=$@

	all_total=`cat $ALL_STATS |  awk -v total=0 '{total=total+$(NF-1); print total}' | tail -1`

	echo "" | sed 1d > temp_stats_parse.txt
	echo $TYPE"	TOTAL_M_READS	PCNT_OF_TOTAL	POOL"
	for m in $ALL
	do	
		m_total=`cat $ALL_STATS | grep $m | awk -v total=0 '{total=total+$NF; print total}' | tail -1`
		pool=`cat $ALL_POOLS | grep $m  | cut -d";" -f1`
		cat $ALL_STATS | grep $m | awk -v all_total=$all_total -v total=0 -v m=$m -v pool=$pool '{total=total+$(NF-1); printf "%s\t%s\t%s\t%s\t%s\n", $1, $2, total, total/all_total*100, pool}' | tail -1 | sed 's/,//g' >> temp_stats_parse.txt
	done
	cat temp_stats_parse.txt | sort -k3nr

}

rev_comp ()
{

	I2=$1
	
	echo $I2 | awk '{
						y=length($1)
						for(i=y+1; i>0; i--)
						{
							b=substr($1, i, 1)
							if(b=="A")
								p="T"							
							if(b=="T")
								p="A"							
							if(b=="C")
								p="G"							
							if(b=="G")
								p="C"	
							printf "%s", p		
						}			
					}' 

}
########################################################################

CONFIG_FILE=`extract_option -conf "/idi/moc_ec//MOC/config_files/PC_config.yaml" 1 $SCRIPT_OPTIONS`
TEMP_PATH=`config_read $CONFIG_FILE Temp_path`
DB_SCRIPT=`config_read $CONFIG_FILE DB_script`
ORPHAN_SCRIPT=`config_read $CONFIG_FILE orphan_script`

########################################################################


DB_MOCS_ID=`echo $MOCS_ID  | sed 's/MOCS\-/MOCS/g'`

echo $DB_MOCS_ID

USID=`USID`

OUT_DIR=$TEMP_PATH$MOCS_ID"/"
TEMP_DIR=$TEMP_PATH$MOCS_ID"/"$USID"/"
TEMP_RTSDB_FILE=$TEMP_DIR"/temp_RtSDB.txt"
ALLPOOL_FILE=$TEMP_DIR"/allpool.txt"

echo "" | sed 1d > $ALLPOOL_FILE

mkdir -p $OUT_DIR
mkdir -p $TEMP_DIR


ALL_STATS=$OUT_DIR"/"$MOCS_ID"_allstats.txt"
METRICS_FILE=$OUT_DIR"/"$MOCS_ID"_walkup_metrics.txt"

SEG_FILE=$TEMP_DIR$MOCS_ID"SGE_file.txt"

echo "" | sed 1d > $SEG_FILE

ALL_FC=`ls -lrt $SEQ_DIR  | grep fastq | grep -e "unmapped" | awk '{print $9}' | cut -d"_" -f1  | cut -d"." -f1 | sort | uniq`
ALL_WU_INDEXES=`ls -lrt $SEQ_DIR  | grep -e "unmapped" -e "unmatched" | awk '{print $9}' | cut -d"." -f3 | sort | uniq`
ALL_LANES=`ls -lrt $SEQ_DIR  | grep -e "unmapped" | awk '{print $9}' | cut -d"." -f2 | sort | uniq`
ALL_SUFFIX=".1.fastq .2.fastq"

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
		ls -lrt $SEQ_DIR"/"$FC*$LANE*fastq* | wc -l
	done
	
done > $METRICS_FILE	


 
for FC in $ALL_FC
do	
	echo "Total files for FC "$FC"..."
	for INDEX in $ALL_WU_INDEXES
	do
		printf "INDEX %s:\t" $INDEX
		ls -lrt $SEQ_DIR"/"*$FC*$INDEX*fastq* | wc -l
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
														if(ar[1]==FC && ar[3]==INDEX && ar[2]==LANE)
															print $NF, FC, INDEX, LANE

														}' | wc -l`
										
								
				
				if [ $NUM_FILES == 1 ];then
				
					FASTQ=`ls -lrt $SEQ_DIR* | grep -v barcode | grep -v md5 | grep fastq | grep $SUFF | awk -v INDEX=$INDEX -v FC=$FC -v SUFF=$SUFF -v LANE=$LANE '{
														
														y=split($9, pr, "/")
														split(pr[y], ar, ".")
														if(ar[1]==FC && ar[3]==INDEX && ar[2]==LANE)
															print $NF

														}'`
								
					ROOT=`basename $FASTQ .fastq.gz`
					COUNT_FILE=$OUT_DIR"/"$ROOT"_counts.txt"
					echo "Launching UGER job for "$FASTQ
					
					echo "echo $FC, $INDEX, $LANE, $SUFF, $FASTQ | tr '\n' ' ' > $COUNT_FILE" > ~/temp.txt
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

#testing that all jobs launched are finished
SGE_test $SEG_FILE	
echo "All Jobs completed...."



total=`cat $OUT_DIR"/"*"_counts.txt" | awk -v total=0 '{total=total+$NF; print total}' | tail -1`
echo $total 
echo "POOL	INDEX	LANE	FILE	TOTAL_M_READS (read 1 only) PCNT_OF_TOTAL_READS" > $ALL_STATS
ALL_COUNTS_FILES=`ls -lrt $OUT_DIR"/"*"_counts.txt" | awk '{print $NF}'`

# 
# ############## Import google sheet databases to server ###############=
# echo "Running $GSIMPORT_SCRIPT to import google sheet databases to server"
# echo "sh $GSIMPORT_SCRIPT"
# sh $GSIMPORT_SCRIPT
# ########################################################################
# 

ls -lrt $MOCS_DB 


for COUNT_FILE in $ALL_COUNTS_FILES
do
	I1=`cat $COUNT_FILE | awk '{print $2}' | sed 's/,//g' | cut -d"_" -f1`
	I2=`cat $COUNT_FILE | awk '{print $2}' | sed 's/,//g' | cut -d"_" -f2`
	
	POOL_NUM=`cat $MOCS_DB | grep $DB_MOCS_ID | grep $I1 | grep $I2 | awk '{print $1}' | sort | uniq | wc -l`
	
	if [ $POOL_NUM == 0 ];then
		POOL_NAME="-"
	else
		POOL_NAME=`cat $MOCS_DB | grep $DB_MOCS_ID | awk -v DB_MOCS_ID=$DB_MOCS_ID -v  I1=$I1 -v I2=$I2  '{
										  	
										  	for(i=2; i < NF+1; i++)
										  	{
										  		y=split($i, ar, ";")
										  		for(t=1; t<y+1;t++)
										  		{
										  			if(ar[1]==DB_MOCS_ID && ar[2]==I1 && ar[3]==I2)
										  				print $1
										  		}
										  
										  
										  	}
										}' | sort | uniq`
	fi
		
	cat $COUNT_FILE | awk -v total=$total -v POOL_NAME=$POOL_NAME '{print POOL_NAME, $2, $3, $5, $NF/1000000, $NF/total*100}' >> $ALL_STATS
done
echo $ALL_COUNTS_FILES
echo $ALL_STATS


ls -lrt $MOCS_DB
ls -lrt $ALL_STATS

printf "%s	%s	%s	%s	" "Flowcell" "TotalMReads"	"MReadsMatchedToIndex:(pcnt)" "MReadsNotMatchedToIndex:(pcnt)" >> $METRICS_FILE	
for LANE in $ALL_LANES
do
	printf "	Lane_%s" $LANE >> $METRICS_FILE	

done 

echo ""  >> $METRICS_FILE	 

READS_F=5
PCNT_F=6

for FC in $ALL_FC
do
	TOTAL=`cat $ALL_STATS | grep $FC | awk -v total_reads=0 -v total_pcnt=0 -v i=0 '{total_reads=total_reads+$'$READS_F';total_pcnt=total_pcnt+$'$PCNT_F';i++; print total_reads}' | tail -1`
	MATCHED=`cat $ALL_STATS | grep $FC | grep -v "unmatched" | awk -v total_reads=0 -v TOTAL=$TOTAL -v i=0 '{total_reads=total_reads+$'$READS_F';total_pcnt=total_pcnt+$'$PCNT_F';i++; print total_reads":("total_reads/TOTAL*100")"}' | tail -1`
	UNMATCHED=`cat $ALL_STATS | grep $FC | grep "unmatched" | awk -v total_reads=0 -v TOTAL=$TOTAL -v i=0 '{total_reads=total_reads+$'$READS_F';total_pcnt=total_pcnt+$'$PCNT_F';i++; print total_reads":("total_reads/TOTAL*100")"}' | tail -1`
	
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


stats_parse $ALL_STATS "Index" $ALLPOOL_FILE $ALL_WU_INDEXES >> $METRICS_FILE

echo "" >> $METRICS_FILE

UNMATCHED_FASTQ=`ls -lrt $SEQ_DIR*.1.unmatched.barcode_1.fastq.gz | awk '{print $9}'`


echo "Top_20_nmatched_Index	%_of_total" >> $METRICS_FILE
echo "sh $ORPHAN_SCRIPT $SEQ_DIR | tail -20 | sort -k1nr >> $METRICS_FILE"
sh $ORPHAN_SCRIPT $SEQ_DIR | tail -20 | sort -k1nr >> $METRICS_FILE


cp $METRICS_FILE $SEQ_DIR

chmod 777 $SEQ_DIR/*_walkup_metrics*

export TMPDIR=$TEMP_PATH

cat $METRICS_FILE |  mailx -a $METRICS_FILE -s "Metrics for $MOCS_ID" $USID"@broadinstitute.org"

change_perms $METRICS_FILE $TEMP_DIR $OUT_DIR