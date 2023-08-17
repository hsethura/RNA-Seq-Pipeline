#!/bin/sh

IN_DIR=$1
OUTDIR_DIR=$2
ID=$3

RESULT_DIR=$IN_DIR"/"
TEMP_DIR=$IN_DIR"/temp_files/"

mkdir $OUTDIR_DIR
FINAL_FILE=$OUTDIR_DIR"/"$ID"_bcLogfile.txt"
 
mkdir $TEMP_DIR
mkdir $RESULT_DIR

index=`ls -lrt $IN_DIR/*logfile.txt | awk '{
										y=split($9, ar, "/")
										z=split(ar[y], pr, ".")
										if(z==6)
											print pr[3]"."pr[4]"."pr[5]
										if(z==4)
											print pr[3]
								}' | sed 's/_frequency/ /g' | awk '{print $1}' | sort | uniq`
								
files=`ls -lrt $IN_DIR/*logfile.txt | awk '{print $9}'`
		

ALL_TEMP_FILES=""
ALL_COMBINED_FILES="" 

echo "Parsing logfiles...."

for m in $files
do	
	echo "Working on file..."
	echo $m
	
	ROOT=`basename $m .txt`
	TEMP_FILE=$TEMP_DIR$ROOT"_temp.txt"
	ALL_INLINES=`cat $m | grep "Barcode:" | awk '{print $2}' | sort -u`
	
	echo "" | sed 1d > $TEMP_FILE
	
	for INLINE in $ALL_INLINES
	do	
		ZMM=`cat $m | grep -A2 $INLINE | awk -F":" '{print $NF}'| tail -1`
		TOT=`cat $m | grep -A4 $INLINE | tail -1 | cut -d":" -f2 | cut -d"(" -f1`
		PCNT=`cat $m | grep -A4 $INLINE| tail -1 | awk '{print $NF}' | sed 's/)//g'`
		
		echo $INLINE	$PCNT	$TOT	$ZMM 
	
	done | sed 's/\%//g' | sort -k2nr > $TEMP_FILE
	
	cat $m | grep "Total non-match reads:" | cut -d":" -f2 | awk '{print "non-match\t"$2"\t"$1"\t-"}' | sed -e 's/(//g' -e 's/)//g'  | sed 's/\%//g'  >> $TEMP_FILE
	cat $m | grep "Total ambiguous reads:" | cut -d":" -f2 | awk '{print "ambiguous\t"$2"\t"$1"\t-"}' | sed -e 's/(//g' -e 's/)//g'  | sed 's/\%//g'  >> $TEMP_FILE
	
	ls -lrt $TEMP_FILE
done



for p in $index
do
	
	echo "Working on index..."
	echo $p
	
	TEMP_FILE=$TEMP_DIR"temp.txt"	
	TEMP_COMBINED_FILE=$TEMP_DIR$p"_logfile_combined_temp.txt"	
	COMBINED_FILE=$RESULT_DIR$p"_logfile_combined.txt"	
	echo "" | sed 1d > $TEMP_COMBINED_FILE
	
	files=`ls -lrt $TEMP_DIR*$p*frequency_logfile_temp.txt | awk '{print $9}'`
	echo $files
	cat $files > $TEMP_FILE
	ALL_INLINE=`cat $TEMP_FILE | awk '{print $1}' | sort | uniq`

	echo $TEMP_FILE
	
	for INLINE in $ALL_INLINE
	do
		
		total_reads=`cat $TEMP_FILE | awk -v total=0 '{total=total+$3; print total}' | tail -1`
		total_perfect=`cat $TEMP_FILE | awk -v total=0 '{if($4 != 0) perfect_reads=$3*($4/100); total=total+perfect_reads; print total}' | tail -1`

		cat $TEMP_FILE | awk -v INLINE=$INLINE -v i=0 -v total_reads=$total_reads -v total_perfect=$total_perfect -v total_pcnt_all=0 -v total_inline=0 -v total_inline_perfect=0 '{
																										
																										if($1==INLINE) 
																										{
																											i++
																											for(i=1; i < NF+1; i++)
																												if($i ~ /nan/)
																													$i=0
																											
																											reads=$3
																											if($4 != 0) 
																												perfect_reads=$3*($4/100)
																											
																											total_inline=total_inline+reads
																											total_inline_perfect=total_inline_perfect+perfect_reads
																										
																											if( total_inline > 0 && total_reads > 0  && total_inline > 0 )
																												print INLINE"	"total_inline/total_reads*100"	"total_inline"	"total_inline_perfect/total_inline*100
																										}
																								
																									
																									}' | tail -1 >> $TEMP_COMBINED_FILE
																									
	done
	
	
	cat $TEMP_COMBINED_FILE | sort -k2nr > $COMBINED_FILE
	echo $COMBINED_FILE
	ALL_COMBINED_FILES=$ALL_COMBINED_FILES" "$COMBINED_FILE
done


echo "Barcode:" > $FINAL_FILE
paste $COMBINED_FILE | sort -k1 | awk '{print $1}' >> $FINAL_FILE
echo "Total_Reads:" >> $FINAL_FILE

barcodes=`cat $FINAL_FILE | sed 1d | grep -v Total_Reads:`

echo $barcodes

for m in $ALL_COMBINED_FILES
do


	name=`echo $m | awk '{y=split($1, ar, "/"); print ar[y]}' | sed 's/_logfile_combined.txt//g'`
	echo $name > $TEMP_DIR"temp.txt"

	total=`cat $m | awk -v total=0 '{total=total+$3; print total}' | tail -1`

	for p in $barcodes
	do
		pcnt=`cat $m | awk -v bc=$p '{if($1==bc) print $2; print "0"}' | sort -k1n | tail -1`
		perf=`cat $m | awk -v bc=$p '{if($1==bc) print $4; print "0"}' | sort -k1n | tail -1`
	
		echo $pcnt >>  $TEMP_DIR"temp.txt"

	done
	
		echo $total >> $TEMP_DIR"temp.txt"
	cat $FINAL_FILE > $TEMP_DIR"temp2.txt"

	paste $TEMP_DIR"temp2.txt" $TEMP_DIR"temp.txt" > $FINAL_FILE


done

echo $FINAL_FILE

echo $FINAL_FILE
