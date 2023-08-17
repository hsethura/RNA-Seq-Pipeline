#!/bin/sh



STRAIN=$1
GENUSES=$2
num_reads=$3
FASTQ_FILE=$4
SAM_PATH=$5
TEMP_PATH=$6
NCBI_PATH=$7
SCRIPT_DIR=$8

source /idi/moc_ec/MOC/scripts/bash_header

### source all functions and set variables 
source "/idi/moc_ec/MOC/scripts/MOC_functions.sh"
paths_and_headers "none" $@

NUM_STRAINS=`extract_option -num_strains 2 1 $@`
ALIGN=`extract_option -align Y 1 $@`


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


align ()
{		
			
			
			FNA_FILE=`ls -lrt $GB_REF_PATH"/"$SPEC"/"$STRAIN"/"* | grep -v rna | grep fna | awk '{print $NF}'` 
			GFF_FILE=`ls -lrt $GB_REF_PATH"/"$SPEC"/"$STRAIN"/"* | grep -v rna | grep gff | awk '{print $NF}'`
			
			#if [ ! -z $FNA_FILE ] && [ ! -z $GFF_FILE ];then
			if [ ! -z $FNA_FILE ];then

				FOUND=`expr $FOUND + 1`
				echo "#############"
				echo $SPEC $STRAIN
				echo $FNA_FILE
				echo $GFF_FILE
				echo "#############"
				echo ""
				SAM_DIR=$SAM_PATH"/"$ROOT"/"$GENUS"/"$SPEC"/"$STRAIN"/"
				mkdir -p $SAM_DIR
				FNA_NAME=`basename $FNA_FILE .fna.gz`
				
				UNZIP_FNA=$TEMP_DIR$FNA_NAME".fna"
				
				zcat $FNA_FILE > $UNZIP_FNA
				
				echo $UNZIP_FNA
				ls -lrt $UNZIP_FNA
				
				TAG=$ROOT"_"$STRAIN"_mapped_"$num_reads
				ALN_SAM=$SAM_DIR"/"$TAG".aligned.sam"
# 				echo $OUT_DIR
# 				echo $ALN_SAM
# 				
				echo "Aligning $SAMPLED_FASTQ to generate $ALN_SAM"
				echo "sh $SCRIPT_DIR"FASTQ_align_GB.sh" $SAMPLED_FASTQ $UNZIP_FNA $TAG $ALN_SAM $TEMP_DIR"
				sh $SCRIPT_DIR"FASTQ_align_GB.sh" $SAMPLED_FASTQ $UNZIP_FNA $TAG $ALN_SAM $TEMP_DIR
				ls -lrt  $ALN_SAM
				

				echo $SPEC $FOUND "******************************"
			fi
}



ROOT=`echo $FASTQ_FILE | awk '{y=split($1, ar, "/"); print ar[y]}' | sed -e 's/.gz//g' -e 's/.fastq//g'`
SAMPLED_FASTQ=$TEMP_DIR$ROOT"_mapped_"$num_reads".fastq1"
SAMPLED_FASTA=`echo $SAMPLED_FASTQ | sed 's/fastq/fasta/g'`
echo "sample_and_strip $FASTQ_FILE $SAMPLED_FASTQ $num_reads"
sample_and_strip $FASTQ_FILE $SAMPLED_FASTQ $num_reads
ls -lrt $SAMPLED_FASTQ

cat $SAMPLED_FASTQ | awk '{if(NR % 4==1) print ">"$1; if(NR % 4==2) print $0}' > $SAMPLED_FASTA
echo $SAMPLED_FASTQ
echo $SAMPLED_FASTA

ALL_GENUSES=`echo $GENUSES | sed 's/,/ /g'`

TEMP_DIR=$TEMP_PATH"/sampled_alignments/"$ROOT"/"
mkdir -p $TEMP_DIR

TEMP_FILE=$TEMP_DIR$ROOT"_temp_sampled_results.txt"
echo $TEMP_FILE

rm $TEMP_FILE
touch $TEMP_FILE


if [ $ALIGN == "Y" ];then
	for GENUS in $ALL_GENUSES
	do
		#echo $GB_REF_PATH
		ALL_SPEC=`ls -lrt $GB_REF_PATH | grep $GENUS | awk '{print $NF}'`
		#echo $ALL_SPEC
		for SPEC in $ALL_SPEC
		do
			FOUND=0
			#echo $SPEC
			ALL_STRAINS=`ls -lrt $GB_REF_PATH"/"$SPEC"/" | awk '{if($5 > 0)print $NF}'`
			#echo $ALL_STRAINS
			for STRAIN in $ALL_STRAINS
			do
			
				if [ $FOUND -lt $NUM_STRAINS ];then
					align
				else
					break
				fi
			
			done
		
		done
	done
fi

if [ $COUNT == "Y" ];then
	for GENUS in $ALL_GENUSES
	do


		ALL_SAMS=`ls -lrt $SAM_PATH"/"$ROOT"/"$GENUS"/"*"/"*"/"*sam | awk '{print $NF}'`
	
		for SAM in $ALL_SAMS
		do
			SPEC=`echo $SAM | rev | cut -d"/" -f3 | rev `
			STRAIN=`echo $SAM | rev | cut -d"/" -f2 | rev `
			FNA_FILE=`ls -lrt $GB_REF_PATH"/"$SPEC"/"$STRAIN"/"* | grep -v rna | grep fna | awk '{print $NF}'` 
			GFF_FILE=`ls -lrt $GB_REF_PATH"/"$SPEC"/"$STRAIN"/"* | grep -v rna | grep gff | awk '{print $NF}'`
			num=`samtools view -F4 $SAM | wc -l | awk -v num_reads=$num_reads '{print $1/num_reads*100}' `
			echo $num $SPEC $STRAIN $SAM $FNA_FILE $GFF_FILE 
		
		
		
		done

	done | sort -k1n
fi


# for GENUS in $ALL_GENUSES
# do
# 	ALL_SAMS=`ls -lrt $SAM_PATH"/"$ROOT"/"$GENUS"/"*"/"*"/"*sam | awk '{print $NF}'`
# 
# 
# 	for SAM in $ALL_SAMS
# 	do
# 		SPEC=`echo $SAM | rev | cut -d"/" -f3 | rev `
# 		STRAIN=`echo $SAM | rev | cut -d"/" -f2 | rev `
# 		FNA_FILE=`ls -lrt $GB_REF_PATH"/"$SPEC"/"$STRAIN"/"* | grep -v rna | grep fna | awk '{print $NF}'` 
# 		GFF_FILE=`ls -lrt $GB_REF_PATH"/"$SPEC"/"$STRAIN"/"* | grep -v rna | grep gff | awk '{print $NF}'`
# 		echo $SAM
# 		echo $GFF_FILE
# 		ls -lrt $FNA_FILE
# 
# 		
# 		
# # 		num=`samtools view -F4 $SAM | wc -l | awk -v num_reads=$num_reads '{print $1/num_reads*100}' `
# # 		echo $SPEC $STRAIN $SAM $num 
# 	
# 	
# 	
# 	done
# 
# done


# SPEC=`cat $TEMP_FILE | sort -k3n | tail -1 | awk '{print $1}'`
# STRAIN=`cat $TEMP_FILE | sort -k3n | tail -1 | awk '{print $2}'`
# 
# echo $SPEC
# 
# echo $STRAIN
# echo $TEMP_FILE
# 
# FNA_FILE=`ls -lrt $GB_REF_PATH"/"$SPEC"/"$STRAIN"/"* | grep -v rna | grep fna | awk '{print $NF}'` 
# GFF_FILE=`ls -lrt $GB_REF_PATH"/"$SPEC"/"$STRAIN"/"* | grep -v rna | grep gff | awk '{print $NF}'`
# 
# echo $FNA_FILE
# 
# ALN_SAM=$TEMP_DIR$ROOT"_"$m"_mapped_"$num_reads".aligned.sam"
# ALN_SAM_ROOT=$TEMP_DIR$ROOT"_"$m"_mapped_"$num_reads
# 
# echo "Aligning $SAMPLED_FASTQ to generate $ALN_SAM"
# echo "sh $SCRIPT_DIR"FASTQ_align_GB.sh" $SAMPLED_FASTQ $FNA_FILE $ALN_SAM_ROOT Y"
# sh $SCRIPT_DIR"FASTQ_align_GB.sh" $SAMPLED_FASTQ $FNA_FILE $ALN_SAM_ROOT Y
# ls -lrt  $ALN_SAM
# 				
# ACC_FILE=$TEMP_PATH"/"$ROOT"_ACC.txt"
# 
# zcat $GFF_FILE | grep -v "#" |  awk '{print $1}' | sort | uniq > $ACC_FILE
# zcat $GFF_FILE | grep -v "#" |  sed 's/ /_/g' |  awk -F"\t" '{
# 	
# 							if($3 != "gene")
# 							{
# 								$NF="type="$3";"$NF
# 								print $0
# 							}
# 						}' > temp.gff
# 
# /idi/moc_ec/MOC/RtS_pipeline/beta/shell_scripts/SAM_PARSE $ALN_SAM coords.txt met.txt $ACC_FILE N
# echo coords.txt met.txt
# 	/broad/IDP-Dx_work/nirmalya/pipeline/beta/shell_scripts/FEATURE_COUNT_COORDS coords.txt temp.gff $OUT_DIR"test.mets" $ACC_FILE N > temp_counts.txt
# 
# cat temp_counts.txt | awk -v total=0 '{total=total+$(NF-1); print total}' | tail -1
# cat temp_counts.txt | grep "type=CDS" | awk -v total=0 '{total=total+$(NF-1); print total}' | tail -1
# cat temp_counts.txt | grep "type=rRNA" | awk -v total=0 '{total=total+$(NF-1); print total}' | tail -1
# 


# 	
# 	
# 	TEMP_FILE=$TEMP_DIR$ROOT"_temp_sampled_results.txt"
# 	OUTPUT=$ALN_DIR$STRAIN"_"$GENUS"_aligntest.txt"
# 
# 	echo "" | sed 1d > $TEMP_FILE
# 	for m in $f
# 	do
# 		if [ -s $NCBI_PATH$m".fna" ];then
# 			
# 			ALN_SAM=$TEMP_DIR$ROOT"_"$m"_mapped_"$num_reads".aligned.sam"
# 			ALN_SAM_ROOT=$TEMP_DIR$ROOT"_"$m"_mapped_"$num_reads
# 			
# 			echo "Aligning $SAMPLED_FASTQ to generate $ALN_SAM"
# 			echo "sh $SCRIPT_DIR"FASTQ_align.sh" $SAMPLED_FASTQ $NCBI_PATH$m".fna" $TEMP_DIR$ROOT"_"$m"_mapped_"$num_reads"
# 			sh $SCRIPT_DIR"FASTQ_align.sh" $SAMPLED_FASTQ $NCBI_PATH$m".fna" $ALN_SAM_ROOT
# 			num=`cat  $ALN_SAM | grep $m | wc -l | awk -v num_reads=$num_reads '{print $1/num_reads}' `
# 			name=`cat $DB_FILE | grep $m | awk '{print $1}' `
# 			echo $name $m $num >> $TEMP_FILE
# 	
# 		else
# 			echo "Can't find "$ACC_PATH$m".fna"
# 		
# 		fi
# 	
# 	
# 	done
# 	
# 			echo  $TEMP_FILE
# 	
# 	echo $STRAIN > $OUTPUT
# 	cat $TEMP_FILE | sort -k3n >> $OUTPUT
# 	
# 	echo $OUTPUT
