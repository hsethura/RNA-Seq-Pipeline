#!/usr/bin/env bash
#
# using getopts
#
uflag=
sflag=
pflag=
tflag=
rflag=
while getopts 'us:p:t:r:' OPTION
do
    case $OPTION in
        u) uflag=1
           ;;
        s) sflag=1
           SUBDIR="$OPTARG"
           ;;
        p) pflag=1
           PROJ_ID="$OPTARG"
           ;;
        t) tflag=1
           TEMP_PATH="$OPTARG"
           ;;
        r) rflag=1
           RESULT_PATH="$OPTARG"
           ;;
        ?) printf "Usage: %s: [-u] [-s subdir] -p proj_id -t temp_path -r result_path args\n" ${0##*/} >&2
           exit 2
           ;;
    esac
done

shift $(($OPTIND - 1))

if [ "$uflag" ]
then
    printf "Option -u for umi specified\n"
fi

if [ "$sflag" ]
then
    printf 'Option -s "%s" for SUBDIR specified\n' "$SUBDIR"
fi

if [ "$pflag" ]
then
    printf 'Option -p "%s" for PROJ_ID specified\n' "$PROJ_ID"
fi

if [ "$tflag" ]
then
    printf 'Option -t "%s" for TEMP_PATH specified\n' "$TEMP_PATH"
fi

if [ "$rflag" ]
then
    printf 'Option -r "%s" for RESULT_PATH specified\n' "$RESULT_PATH"
fi

#TEMP_PATH="/broad/hptmp/MOC/"

USERID=`id | cut -d"(" -f2 | cut -d")" -f1`

#TEMP_DIR=$TEMP_PATH"/"$PROJ_ID"/temp_files/"
#TEMP_DIR=$TEMP_PATH"/"$PROJ_ID"/"$USERID"/RPG_temp_files/"
#COUNT_DIR=$TEMP_PATH"/"$PROJ_ID"/patho_result/"$PROJ_ID"/"
COUNT_DIR=$TEMP_PATH"/patho_result/"$PROJ_ID"/"
RESULT_DIR=$RESULT_PATH"/"$PROJ_ID"/"
if [ "$uflag" ]; then
    COUNT_DIR=$COUNT_DIR"/umidir/"
    RESULT_DIR=$RESULT_DIR"/umidir/"
fi

if [ "$sflag" ]; then 
    COUNT_DIR=$COUNT_DIR"/"$SUBDIR"/"
    RESULT_DIR=$RESULT_DIR"/"$SUBDIR"/"
fi

TEMP_DIR=$COUNT_DIR"/RPG_temp_files/"
MOC_SCRIPT_PATH="/broad/IDP-Dx_storage/MOC/scripts/"


printf 'Count_dir: "%s"\n' "$COUNT_DIR"

mkdir -p $TEMP_DIR
mkdir -p $RESULT_DIR

ls -lrt $TEMP_DIR
ls -lrt $RESULT_DIR


echo $COUNT_DIR
echo $TEMP_DIR
echo $RESULT_DIR



#mkdir $TEMP_DIR

READ_COUNT()
{

	READ_FILE=$1
	shift
	OUTPUT_FILE=$1
	shift
	TOTAL_ACC_COUNTS=$1
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
		
		cat $READ_FILE | awk -v total=0 -v AS_total=0 -v KEY=$KEY -v T=$TOTAL_ACC_COUNTS -v NL=$NL -v ACC=$ACC '{
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
													printf "%s_total_counts_for_replicon:\t%s\n", KEY, total
													
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



ALL_COUNT_FILES=`ls -lrt $COUNT_DIR*".counts" | awk '{print $9}' | sort `
 
for INPUT_FILE in $ALL_COUNT_FILES
do

	ls -lrt $INPUT_FILE

done


ALL_OUTPUT=""


for INPUT_FILE in $ALL_COUNT_FILES
do
	
	echo "Working on analysis of"$INPUT_FILE"...."

	SAMPLE_ID=`echo $INPUT_FILE | awk '{y=split($1, ar, "/"); print ar[y]}' | sed 's/.counts//g'`
	MET_FILE=`echo $INPUT_FILE | sed 's/.counts/.mets/g'`

#####	
	
	typeset -i NUM_ACC
	NUM_ACC=`cat $INPUT_FILE | sed 1d  | awk '{print $2}' | sed 's/;/ /g' |  tr ' ' '\n' | sort | uniq | wc -l`

	ALL_ACC=`cat $INPUT_FILE | sed 1d | awk '{print $2}' | sed 's/;/ /g' |  tr ' ' '\n' | sort | uniq `
	
	
	
	echo  $INPUT_FILE 
	echo $ALL_ACC
	
	if [ $NUM_ACC -gt 1 ];then
	
		ALL_ACC="ALL "$ALL_ACC
	fi

	for ACC in $ALL_ACC
	do
		echo $MET_FILE
		
		ACC_TEMP_FILE=$TEMP_DIR$SAMPLE_ID"_"$ACC"_temp.counts"
		TEMP_FILE=$TEMP_DIR$SAMPLE_ID"_temp.counts"
		
		cat $INPUT_FILE | sed 's/AS:/AS_/g' | awk -v ACC=$ACC '{if(ACC=="ALL" || $2==ACC) print $0}' > $ACC_TEMP_FILE
		cat $INPUT_FILE | sed 's/AS:/AS_/g'  > $TEMP_FILE
		
		OUTPUT_FILE=$TEMP_DIR$SAMPLE_ID"_"$ACC"_metrics.txt"
		echo "" | sed 1d > $OUTPUT_FILE

		typeset -f TOTAL_ACC_COUNTS
		
		TOTAL_READS=`cat $MET_FILE | grep "Total_reads:" | awk '{print $2}' | sed 's/\.//g'`
		PCNT_PAIRS=`cat $MET_FILE | grep _reads_properly_mapped_in_pairs | awk '{print $2}' `
		PCNT_MAPPED=`cat $MET_FILE | grep mapped_to_any_reference | awk '{print $2}' `
		INST_LEN=`cat $MET_FILE | grep Avge_insert_length | awk '{print $2}' `
		
		if [ $NUM_ACC -gt 1 ] && [ $ACC != "ALL" ];then
			INST_LEN="-"
		fi
		
		echo $INST_LEN
		echo $ACC
				
		
		echo $TOTAL_FRAGS
		
		if [ -s $ACC_TEMP_FILE ];then
			TOTAL_ACC_COUNTS=`cat $ACC_TEMP_FILE | awk -v total=0 '{total=total+$3;print total}' | tail -1`
		else
			TOTAL_ACC_COUNTS=0
			echo "No counts for replicon "$ACC
			ls -lrt $ACC_TEMP_FILE
		fi
		 
		echo $TOTAL_ACC_COUNTS
		
		
		if [ $TOTAL_ACC_COUNTS != 0 ];then
			
			echo "...for replicon "$ACC
			
			if [ $ACC == "ALL" ];then
				TOTAL_ACC_ALGN=`cat $MET_FILE | grep Total_reads_aligned_to_  | awk -v total=0 '{total=total+$2;print total}' | sed 's/\.//g' | tail -1`						
			else
				TOTAL_ACC_ALGN=`cat $MET_FILE | grep Total_reads_aligned_to_ | grep $ACC | awk '{print $2}' | sed 's/\.//g'`						
			fi

			TOTAL_ACC_COUNTED=`cat $ACC_TEMP_FILE | sed 1d | awk -v total=0 '{total=total+$3;print total}' | tail -1`
			SENSE_READS=`cat $ACC_TEMP_FILE | sed 1d | awk -v total=0 '{split($1, ar, ":");split(ar[1], pr, "_");if(ar[1] !~ /IGR/ && pr[1]!="AS")total=total+$NF; print total}' | tail -1`
			ANTISENSE_READS=`cat $ACC_TEMP_FILE | sed 1d | awk -v total=0 '{split($1, ar, ":");split(ar[1], pr, "_");if(ar[1] !~ /IGR/ && pr[1]=="AS")total=total+$NF; print total}' | tail -1`
			PCNT_ALIGNED=`echo $TOTAL_READS $TOTAL_ACC_ALGN | awk '{if($1>0) print $2/$1 * 100; if($1==0) print "0"}'`
			PCNT_SENSE=`echo $TOTAL_ACC_COUNTS $SENSE_READS $ANTISENSE_READS | awk '{if($2+$3>0)print $2/($2+$3)*100; if($2+$3==0)print "0"}'`
			PCNT_COUNTED=`echo $TOTAL_READS $TOTAL_ACC_ALGN | awk '{if($1>0) print $2/$1 * 100; if($1==0) print "0"}'`
			
			echo $MET_FILE
			
			echo $TOTAL_READS " TOTAL_READS"
			echo $TOTAL_ACC_ALGN " TOTAL_ACC_ALGN"
			echo $PCNT_ALIGNED " PCNT_ALIGNED"
			echo $TOTAL_ACC_COUNTED " TOTAL_ACC_COUNTED"
			
			
			echo "Sample:	"$SAMPLE_ID >> $OUTPUT_FILE
			echo "Total_reads:	"$TOTAL_READS >> $OUTPUT_FILE
			echo "For_replicon...	"$ACC >> $OUTPUT_FILE
			echo "pcnt_aligned:	"$PCNT_ALIGNED >> $OUTPUT_FILE
			echo "pcnt_properly_mapped_pairs:	"$PCNT_PAIRS >> $OUTPUT_FILE
			echo "average_insert_len:	"$INST_LEN >> $OUTPUT_FILE

			echo "Total_frags_counted	"$TOTAL_ACC_COUNTED >> $OUTPUT_FILE
			echo "pcnt_sense:	"$PCNT_SENSE >> $OUTPUT_FILE
		
			#echo "READ_COUNT" $ACC_TEMP_FILE $OUTPUT_FILE $TOTAL_ACC_COUNTS $ACC" CDS rRNA misc_RNA tRNA IGR"
			READ_COUNT $ACC_TEMP_FILE $OUTPUT_FILE $TOTAL_ACC_COUNTED $ACC CDS rRNA misc_RNA ncRNA tRNA IGR
			
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

	echo $i, $m, $NL
	if [ $i == 0 ];then
		cat $m | awk '{printf "%s\t", $1}' | sed 's/://g' > $ALL_METRICS_FILE  && echo "" >> $ALL_METRICS_FILE
		cat $m | awk '{printf "%s\t", $2}' | sed 's/://g' >> $ALL_METRICS_FILE && echo "" >> $ALL_METRICS_FILE
	fi
	if [ $i != 0 ];then
		cat $m | awk -v NL=$NL '{printf "%s\t", $2}' | sed 's/://g' >> $ALL_METRICS_FILE && echo "" >> $ALL_METRICS_FILE
	fi
	 
	i=`expr $i + 1`
done 

echo $ALL_METRICS_FILE


#sh $MOC_SCRIPT_PATH"MOC_KeyMets_Join.sh" $PROJ_ID $RESULT_PATH Sample

