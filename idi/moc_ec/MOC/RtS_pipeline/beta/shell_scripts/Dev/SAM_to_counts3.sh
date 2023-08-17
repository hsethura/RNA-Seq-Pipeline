#!/bin/sh


#### SAM MUST BE SORTED BY KEY, NOT COORD

SAM=$1
FEATURE_FILE=$2
SCRIPT_DIR=$3"/"
TEMP_PATH=$4"/"
COUNT_FILE=$5

### source all functions 
source "/broad/IDP-Dx_storage/MOC/scripts/MOC_functions.sh"


#######   feature_parser  ######

feature_parser ()
{
	FEATURE_FILE=$1
	OUTFILE=$2
	STRAND=$3

	cat $FEATURE_FILE | grep -v "#" | sed 's/ /_/g' | awk '{
															
																y=split($9, ar, ";")
																for(i=1; i<y+1; i++)
																{
																	if(ar[i] ~ /type=/)
																	{
																		split(ar[i], pr, "=")
																		split(pr[2], qr, ":")
																		$3=qr[1]
																	}
																}
																# print sense
																print $0													
															
																if($7=="+")
																	strand="-"
																if($7=="-")
																	strand="+"
															
																$7=strand
																$3="AS_"$3
																# print antisense
																print $0
														
															}' | awk -v STRAND=$STRAND '{if($7==STRAND) print $0}' | sort -k4n  > $OUTFILE

}
################################

options=$@

STRAND_REV=`extract_option -STRAND_REV N 1 $options`
UNIQ_READS=`extract_option -UNIQ_READS N 1 $options`

echo $STRAND_REV


MET_FILE=`echo $COUNT_FILE | sed 's/.counts/.mets/g'`

ROOT=`echo $SAM | rev | cut -d"/" -f1 | rev | sed 's/.sam//g'`
COORD_FILE=$TEMP_PATH"/"$ROOT"_coords.txt"
COORD_UNIQ_FILE=$TEMP_PATH"/"$ROOT"_coords_uniq.txt"
COUNT_TEMP_FILE=$TEMP_PATH"/"$ROOT"_counts_temp.txt"
ACC_FILE=$TEMP_PATH"/"$ROOT"_ACC.txt"
FEATURE_PLUS_PARSED=$TEMP_PATH"/"$ROOT"_parsed_+.gff"
FEATURE_MINUS_PARSED=$TEMP_PATH"/"$ROOT"_parsed_-.gff"

cat $FEATURE_FILE | grep -v "#" |  awk '{print $1}' | sort | uniq > $ACC_FILE

echo "Pulling out coordinates and metrics..."
echo "$SCRIPT_DIR"SAM_PARSE" $SAM $COORD_FILE $MET_FILE $ACC_FILE $STRAND_REV"
$SCRIPT_DIR"SAM_PARSE" $SAM $COORD_FILE $MET_FILE $ACC_FILE $STRAND_REV 
# echo $COORD_FILE
# 
 
feature_parser $FEATURE_FILE $FEATURE_PLUS_PARSED "+"
feature_parser $FEATURE_FILE $FEATURE_MINUS_PARSED "-"

echo $FEATURE_PLUS_PARSED
echo $FEATURE_MINUS_PARSED
echo $FEATURE_FILE

if [ $UNIQ_READS == "Y" ];then

	echo "Finding unique coordinates..."
	$COORD_FILE | awk '{split($NF, ar, ":");$NF=ar[1]":a"r[2]":"ar[3]; print $0}' | sort | uniq  > $COORD_UNIQ_FILE
	cat $COORD_UNIQ_FILE > $COORD_FILE

fi

echo "Generating counts for features..."
$SCRIPT_DIR"FEATURE_COUNTS_COORDS_2" $COORD_FILE $FEATURE_PLUS_PARSED $MET_FILE $ACC_FILE > $COUNT_TEMP_FILE
$SCRIPT_DIR"FEATURE_COUNTS_COORDS_2" $COORD_FILE $FEATURE_MINUS_PARSED $MET_FILE $ACC_FILE >> $COUNT_TEMP_FILE


cat $COUNT_TEMP_FILE | grep -v "#" | awk '{

									y=split($2, ar, ";")
									h=split(ar[y], tr, ":")
									acc=tr[h]
									
									for(i=1; i < y+1; i++)
									{
										split(ar[i], pr, "=")
										if(pr[1]=="full_tag")
										{	
											if($1 !~ /AS_/)	
												printf "%s:%s\t%s\t%.0f\n" ,pr[2], acc, acc, $3
											if($1 ~ /AS_/)	
												printf "AS_%s:%s\t%s\t%.0f\n" ,pr[2], acc, acc, $3
										}
									}
									
									
								}' | sort > $COUNT_FILE
								
								
echo $COUNT_FILE
echo $MET_FILE
echo $COUNT_TEMP_FILE
   
