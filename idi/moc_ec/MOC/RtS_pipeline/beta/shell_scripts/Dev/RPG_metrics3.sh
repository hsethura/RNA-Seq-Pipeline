#!/bin/sh

PROJ_ID=$1
TEMP_PATH=$2
RESULT_PATH=$3

COUNT_DIR=$TEMP_PATH"/"$PROJ_ID"/patho_result/"$PROJ_ID"/"
RESULT_DIR=$RESULT_PATH"/"$PROJ_ID"/"
TEMP_DIR=$TEMP_PATH"/"$PROJ_ID"/temp_files/"

mkdir $TEMP_DIR

echo $COUNT_DIR
echo $RESULT_DIR
echo $TEMP_DIR

mkdir $TEMP_DIR

READ_COUNT()
{

	READ_FILE=$1
	shift
	OUTPUT_FILE=$1
	shift
	TOTAL_COUNTS=$1
	shift
	ACC=$1
	shift
	
	for m in $@
	do
	
		KEY=$m
		AS_KEY="AS_"$m
		
		typeset -i NL
		
		NL=`cat $READ_FILE | wc -l`
				
		if [ $NL == 0 ];then
			echo "0" >> $READ_FILE
		fi
		
		cat $READ_FILE | awk -v total=0 -v AS_total=0 -v KEY=$KEY -v T=$TOTAL_COUNTS -v NL=$NL -v ACC=$ACC '{
											split($1, ar, ":")
											
											AS_KEY="AS_"KEY
											if($2==ACC || ACC=="ALL")
											{
												if(ar[1]==KEY)	
													total=total+$NF
												if(ar[1]==AS_KEY)	
													AS_total=AS_total+$NF
											}
											if(NR==NL)
											{
												if(KEY!="CDS" && KEY!="IGR")
												{
													printf "%s_pcnt_of_counted:\t", KEY
													if(T > 0)
														printf "%s\n", total/T*100
													else
														printf "0\n"
													printf "%s_pcnt_sense:\t", KEY
													if(total+AS_total > 0)
														printf "%s\n", total/(total+AS_total)*100
													else 
														printf "0\n"
												
												}
												if(KEY=="CDS")
												{
													printf "%s_total_counts:\t%s\n", KEY, total
													
													printf "%s_pcnt_of_counted:\t", KEY
													if(T > 0)
														printf "%s\n", total/T*100
													else
														printf "0\n"
													printf "%s_pcnt_sense:\t", KEY
													if(total+AS_total > 0)
														printf "%s\n", total/(total+AS_total)*100
													else 
														printf "0\n"
													
												}
												if(KEY=="IGR")
												{
													printf "%s_pcnt_of_counted:\t", KEY
													if(T > 0)
														printf "%s\n", (total+AS_total)/T*100
													else 
														printf "0\n"
												}
											} 
								
								}' >> $OUTPUT_FILE

	done

}




ALL_COUNT_FILES=`ls -lrt $COUNT_DIR*".counts" | grep -v ".s.counts" | grep -v ".as.counts"  | awk '{print $9}' | sort `
 

ALL_OUTPUT=""


for m in $ALL_COUNT_FILES
do
	
	echo "Working on analysis of"$m"...."
	INPUT_FILE=$m
	SAMPLE_ID=`echo $INPUT_FILE | awk '{y=split($1, ar, "/"); print ar[y]}' | sed 's/.counts//g' | sed 's/_All//g'`
	SUMMARY_FILE=`echo $INPUT_FILE | sed 's/.counts/_s.counts.summary/g'`

#####	
	
	typeset -i NUM_ACC
	NUM_ACC=`cat $INPUT_FILE | sed 1d | sed 1d | awk '{print $2}' | sed 's/;/ /g' |  tr ' ' '\n' | sort | uniq | wc -l`

	ALL_ACC=`cat $INPUT_FILE | sed 1d | sed 1d | awk '{print $2}' | sed 's/;/ /g' |  tr ' ' '\n' | sort | uniq `
	
	if [ $NUM_ACC -gt 1 ];then
	
		ALL_ACC="ALL "$ALL_ACC
	fi
	
	for ACC in $ALL_ACC
	do
		ACC_TEMP_FILE=$TEMP_DIR$SAMPLE_ID"_"$ACC"_temp.counts"
		TEMP_FILE=$TEMP_DIR$SAMPLE_ID"_temp.counts"
		
		cat $INPUT_FILE | sed 's/AS:/AS_/g' | awk -v ACC=$ACC '{if(ACC=="ALL" || $2==ACC) print $0}' > $ACC_TEMP_FILE
		cat $INPUT_FILE | sed 's/AS:/AS_/g'  > $TEMP_FILE
		
		OUTPUT_FILE=$TEMP_DIR$SAMPLE_ID"_"$ACC"_metrics.txt"
		echo "" | sed 1d > $OUTPUT_FILE

		typeset -i TOTAL_COUNTS
		
		TOTAL_FRAGS=`cat $SUMMARY_FILE | awk -v total=0 '{total=total+$2; print total}' | tail -1`
		
		if [ -s $ACC_TEMP_FILE ];then
			TOTAL_COUNTS=`cat $ACC_TEMP_FILE | sed 1d | awk -v total=0 '{total=total+$NF; print total}' | tail -1`
		else
			TOTAL_COUNTS=0
			echo "No counts for replicon "$ACC
			ls -lrt $ACC_TEMP_FILE
		fi
		
		if [ $TOTAL_COUNTS != 0 ];then
			
			echo "...for replicon "$ACC
			
			SENSE_READS=`cat $TEMP_FILE | sed 1d | awk -v total=0 '{split($1, ar, ":");split(ar[1], pr, "_");if(ar[1] !~ /IGR/ && pr[1]!="AS")total=total+$NF; print total}' | tail -1`
			ANTISENSE_READS=`cat $TEMP_FILE | sed 1d | awk -v total=0 '{split($1, ar, ":");split(ar[1], pr, "_");if(ar[1] !~ /IGR/ && pr[1]=="AS")total=total+$NF; print total}' | tail -1`
			PCNT_ALIGNED=`echo $TOTAL_FRAGS $TOTAL_COUNTS | awk '{if($1>0) print $2/$1 * 100; if($1==0) print "0"}'`
			PCNT_SENSE=`echo $TOTAL_COUNTS $SENSE_READS $ANTISENSE_READS | awk '{if($2+$3>0)print $2/($2+$3)*100; if($2+$3==0)print "0"}'`
					
			echo "Sample:	"$SAMPLE_ID"_"$ACC >> $OUTPUT_FILE
			echo "total_fragments:	"$TOTAL_FRAGS >> $OUTPUT_FILE
			echo "total_counted:	"$TOTAL_COUNTS >> $OUTPUT_FILE
			echo "pcnt_counted:	"$PCNT_ALIGNED >> $OUTPUT_FILE
			echo "pcnt_sense:	"$PCNT_SENSE >> $OUTPUT_FILE
		
			#echo "READ_COUNT" $ACC_TEMP_FILE $OUTPUT_FILE $TOTAL_COUNTS $ACC" CDS rRNA misc_RNA tRNA IGR"
			READ_COUNT $ACC_TEMP_FILE $OUTPUT_FILE $TOTAL_COUNTS $ACC CDS rRNA misc_RNA ncRNA tRNA IGR
			
			ALL_OUTPUT=$ALL_OUTPUT" "$OUTPUT_FILE
		fi
	
	
	done

######	


done


i=0

typeset -i NL
NL=`cat $OUTPUT_FILE | wc -l`

ALL_METRICS_FILE=$RESULT_DIR$PROJ_ID"_metrics.txt"



echo $ALL_METRICS_FILE 
echo $ALL_OUTPUT 


for m in $ALL_OUTPUT
do

	echo $i, $m
	if [ $i == 0 ];then
		cat $m | awk '{printf "%s\t", $1}' | sed 's/://g' > $ALL_METRICS_FILE  && echo "" >> $ALL_METRICS_FILE
		cat $m | awk '{printf "%s\t", $2}' | sed 's/://g' >> $ALL_METRICS_FILE && echo "" >> $ALL_METRICS_FILE
	fi
	if [ $i != 0 ];then
		cat $m | awk -v NL=$NL '{printf "%s\t", $2}' | sed 's/://g' >> $ALL_METRICS_FILE && echo "" >> $ALL_METRICS_FILE
	fi
	
	i=`expr $i + 1`
done 


