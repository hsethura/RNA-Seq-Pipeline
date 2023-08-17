#!/bin/sh

STRAIN=$1
shift
GENUS=$1
shift
num_reads=$1
shift
FASTQ_FILE=$1
shift
ALN_DIR=$1
shift
TEMP_DIR=$1
shift
NCBI_PATH=$1
shift
SCRIPT_DIR=$1
shift

sample_and_strip ()
{

	INPUT=$1
	OUTFILE=$2
	num_reads=$3
	
	echo "sample_and_strip"
	GZ=`file $FASTQ_FILE | grep gz | wc -l`
	if [ $GZ == 0 ];then
		cat $INPUT | awk -v num_reads=$num_reads '{
															if(NR % 4 == 2 || NR % 4 == 0)
															{	
																len=length($1)
																t=substr($1, 7, len+1)
																print t
															
															}
															else
																print $1

															
															if(NR>num_reads*4) 
																exit; 
														
													}'  > $OUTFILE
	else
		zcat $INPUT | awk -v num_reads=$num_reads '{
															if(NR % 4 == 2 || NR % 4 == 0)
															{	
																len=length($1)
																t=substr($1, 7, len+1)
																print t
															}
															else
																print $1
															if(NR>num_reads*4) 
																exit; 
														
													}'  > $OUTFILE
	fi


}

DB_FILE="/broad/IDP-Dx_storage/BactRAP/scripts/all_NCBI_DB.txt"


OUTPUT=$ALN_DIR$STRAIN"_"$GENUS"_aligntest.txt"

ROOT=`echo $FASTQ_FILE | awk '{y=split($1, ar, "/"); print ar[y]}' | sed -e 's/.gz//g' -e 's/.fastq//g'`
down_file=$TEMP_DIR$ROOT"downsample_file.txt"
sam=$TEMP_DIR$ROOT".sam"
sampled_sam=$TEMP_DIR$ROOT"_mapped_"$num_reads".sam"
sampled_fastq1=$TEMP_DIR$ROOT"_mapped_"$num_reads".fastq1"

echo $ROOT

GZ=`echo $FASTQ_FILE | grep gz | wc -l`

echo $GZ, $FASTQ_FILE, $num_reads

sample_and_strip $FASTQ_FILE  $sampled_fastq1 $num_reads

ls -lrt $sampled_fastq1
ls -lrt $FASTQ_FILE

echo $sampled_fastq1

f=`cat $DB_FILE | grep $GENUS | grep chromo | awk '{print $2}'`

echo $f

echo "cat $DB_FILE | grep $GENUS | grep chromo | awk '{print $2}'"

TEMP_FILE=$TEMP_DIR"temp_sampled_results.txt"

echo "" | sed 1d > $TEMP_FILE
for m in $f
do
	if [ -s $NCBI_PATH$m".fna" ];then
		
		ALN_SAM=$TEMP_DIR$ROOT"_"$m"_mapped_"$num_reads".aligned.sam"
		num=`cat  $ALN_SAM | grep $m | wc -l | awk -v num_reads=$num_reads '{print $1/num_reads}' `
		name=`cat $DB_FILE | grep $m | awk '{print $1}' `
		echo "Aligning $sampled_fastq1 to generate  $TEMP_DIR$ROOT"_"$m"_mapped_"$num_reads".aligned.sam:""
		echo "sh $SCRIPT_DIR"FASTQ_align.sh" $sampled_fastq1 $NCBI_PATH$m".fna" $TEMP_DIR$ROOT"_"$m"_mapped_"$num_reads"
		sh $SCRIPT_DIR"FASTQ_align.sh" $sampled_fastq1 $NCBI_PATH$m".fna" $TEMP_DIR$ROOT"_"$m"_mapped_"$num_reads
		echo $name $m $num >> $TEMP_FILE
	
	else
		echo "Can't find "$ACC_PATH$m".fna"
	
	fi

done

echo  $TEMP_FILE


echo $STRAIN > $OUTPUT
cat $TEMP_FILE | sort -k3n >> $OUTPUT

echo $OUTPUT