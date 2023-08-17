#!/bin/sh



for SEQ_DIR in $@
do
	OUT_DIR=$SEQ_DIR

	SEQ_ID=`basename $SEQ_DIR`

	ALL_STATS=$OUT_DIR"/"$SEQ_ID"all_stats.txt"
	METRICS_FILE=$OUT_DIR"/"$SEQ_ID"_walkup_metrics.txt"

	echo $SEQ_ID

	mkdir -p $SEQ_DIR
	SEG_FILE="SGE_file.txt"

	echo "" | sed 1d > $SEG_FILE

	############################  SGE_test ###############################

	SGE_test ()
	{	
		SGE_OUT=$1
	
		jobs_submit=`cat $SGE_OUT`
	
		typeset -i NUM
	
		while :
		do
			
			job_running=`qstat | awk '{if(NR > 2 ) print $1}'`
			NUM=`echo $job_running $jobs_submit | tr ' ' '\n' | sort | uniq -c | awk '{if($1==2)print $0}' | wc -l`
			if [ $NUM == 0 ];then
		
				echo "qstat jobs are complete"
				break
			fi
		
		done
	}

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
			cat $ALL_STATS | grep $m | awk -v all_total=$all_total -v total=0 -v m=$m '{total=total+$(NF-1); print m"\t"total"\t"total/all_total*100}' | tail -1 >> temp_stats_parse.txt
		done
		cat temp_stats_parse.txt | sort -k2nr

	}

	########################################################################



	ALL_FC=`ls -lrt $SEQ_DIR  | grep -e "unmapped" -e "unmatched" | awk '{print $9}' | cut -d"_" -f2  | cut -d"." -f1 | sort | uniq`
	ALL_INDEX=`ls -lrt $SEQ_DIR  | grep -e "unmapped" -e "unmatched" | awk '{print $9}' | cut -d"." -f3 | sort | uniq`
	ALL_LANES=`ls -lrt $SEQ_DIR  | grep -e "unmapped" -e "unmatched" | awk '{print $9}' | cut -d"_" -f1 | sort | uniq`
	#ALL_SUFFIX=".1.fastq .2.fastq"
	ALL_SUFFIX=".1.fastq"

	NUM_INDEX=`echo $ALL_INDEX | wc -w`

	echo $ALL_INDEX 
	echo $NUM_INDEX 
	echo $ALL_FC 


	typeset -i NUM_FILES

	echo "" | sed 1d > temp.txt

	for FC in $ALL_FC
	do
		for INDEX in $ALL_INDEX
		do
			for LANE in $ALL_LANES
			do
				for SUFF in $ALL_SUFFIX
				do
					NUM_FILES=`ls -lrt $SEQ_DIR"/"$LANE* | grep $INDEX | grep $FC | grep $SUFF | grep -v barcode |  awk '{print $5}' | wc -l`
					if [ $NUM_FILES == 1 ];then
					
						FASTQ=`ls -lrt $SEQ_DIR"/"$LANE* | grep $INDEX | grep $FC | grep $SUFF | grep -v barcode |  awk '{print $9}' `
						ROOT=`basename $FASTQ .fastq.gz`
						COUNT_FILE=$SEQ_DIR"/"$ROOT"_counts.txt"
						echo "Launching UGER job for "$FASTQ
					
						echo "echo $FC, $INDEX, $LANE, $SUFF, $FASTQ | tr '\n' ' ' > $COUNT_FILE" >temp.txt
						echo "zcat $FASTQ | grep @ | wc -l >> $COUNT_FILE" >> temp.txt

						qsub temp.txt >> $SEG_FILE					
					
					fi
				done
			done
		done
	done

	echo $SEG_FILE
	qstat
	echo "Jobs submitted - waiting for them to complete...."
	SGE_test $SEG_FILE	


	total=`cat $SEQ_DIR"/"*"_counts.txt" | awk -v total=0 '{total=total+$NF; print total}' | tail -1`
	echo $total 
	echo "INDEX	LANE	FILE	TOTAL_M_READS PCNT_OF_TOTAL_READS" > $ALL_STATS
	cat $SEQ_DIR"/"*"_counts.txt" | awk -v total=$total '{print $2, $3, $5, $NF/1000000, $NF/total*100}' >> $ALL_STATS

	echo $ALL_STATS


	printf "%s	%s	%s	%s	" "Flowcell" "TotalMReads"	"MReadsMatchedToIndex:(pcnt)" "MReadsNotMatchedToIndex:(pcnt)" > $METRICS_FILE	
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


	stats_parse $ALL_STATS "Index" $ALL_INDEX

	echo "" 



	echo $ALL_STATS 
	echo $METRICS_FILE
done