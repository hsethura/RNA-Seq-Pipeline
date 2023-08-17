#!/bin/sh


#### SAM MUST BE SORTED BY KEY, NOT COORD

SAM=$1
FEATURE_FILE=$2
SCRIPT_DIR=$3"/"
TEMP_PATH=$4"/"
COUNT_FILE=$5

########  function for extracting options from command line
 
		extract_option ()
		{
			option_name=$1 
			shift
			default=$1
			shift 
			echo $@ | awk -v def=$default -v name=$option_name '{
		
				for(i=1; i < NF+1; i++)
				{
					if(i==1)			# intitalize to $default in case no option found
						print def
					if($i==name)
						print $(i+1)
				}
		
			}' | tail -1
		}

###############################  extract and checkpoint options  ##############################

options=$@

STRAND_REV=`extract_option -STRAND_REV N $options`
UNIQ_READS=`extract_option -UNIQ_READS N $options`

MET_FILE=`echo $COUNT_FILE | sed 's/.counts/.mets/g'`

ROOT=`echo $SAM | rev | cut -d"/" -f1 | rev | sed 's/.sam//g'`
COORD_FILE=$TEMP_PATH"/"$ROOT"_coords.txt"
COORD_UNIQ_FILE=$TEMP_PATH"/"$ROOT"_coords_uniq.txt"
COUNT_TEMP_FILE=$TEMP_PATH"/"$ROOT"_counts_temp.txt"
ACC_FILE=$TEMP_PATH"/"$ROOT"_ACC.txt"
FEATURE_FILE_PARSED=$TEMP_PATH"/"$ROOT"_parsed.gff"


cat $FEATURE_FILE | grep -v "#" |  awk '{print $1}' | sort | uniq > $ACC_FILE

echo "Pulling out coordinates and metrics..."
echo "$SCRIPT_DIR"SAM_PARSE" $SAM $COORD_FILE $MET_FILE $ACC_FILE $STRAND_REV"
$SCRIPT_DIR"SAM_PARSE" $SAM $COORD_FILE $MET_FILE $ACC_FILE $STRAND_REV 
# echo $COORD_FILE
# 
 cat $FEATURE_FILE | grep -v "#" | sed 's/ /_/g' > $FEATURE_FILE_PARSED

if [ $UNIQ_READS == "Y" ];then

	echo "Finding unique coordinates..."
	$COORD_FILE | awk '{split($NF, ar, ":");$NF=ar[1]":a"r[2]":"ar[3]; print $0}' | sort | uniq  > $COORD_UNIQ_FILE
	cat $COORD_UNIQ_FILE > $COORD_FILE

fi

echo $FEATURE_FILE_PARSED

echo "Generating counts for features..."
echo "$SCRIPT_DIR"FEATURE_COUNT_COORDS" $COORD_FILE $FEATURE_FILE_PARSED $MET_FILE $ACC_FILE > $COUNT_TEMP_FILE"
 $SCRIPT_DIR"FEATURE_COUNT_COORDS" $COORD_FILE $FEATURE_FILE_PARSED $MET_FILE $ACC_FILE > $COUNT_TEMP_FILE


cat $COUNT_TEMP_FILE | grep -v "#" | awk '{

									y=split($1, ar, ";")
									h=split(ar[y], tr, ":")
									acc=tr[h]
									
									for(i=1; i < y+1; i++)
									{
										split(ar[i], pr, "=")
										if(pr[1]=="full_tag")
										{
											if($NF=="AS")
												pr[2]="AS_"pr[2]
											
# 											print pr[1]
# 											print pr[2]
											printf "%s:%s\t%s\t%.0f\n" ,pr[2], acc, acc, $2
										}
									
									}
									
									
								}'  > $COUNT_FILE
								
								
echo $COUNT_FILE
echo $MET_FILE
echo $COUNT_TEMP_FILE
   
