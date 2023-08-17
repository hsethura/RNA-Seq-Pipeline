#!/bin/sh


source "/broad/software/scripts/useuse"
reuse -q UGER
reuse -q Java-1.8
reuse -q Samtools
reuse -q BWA
reuse -q R-3.2
reuse -q GCC-5.2

EXPID=$1
shift
FASTQ1=$1
shift

FASTQ2=`echo $FASTQ1 | sed 's/_R1_/_R2_/'g`
FASTQ2=`echo $FASTQ1 | sed 's/.1.fastq/.2.fastq/'g`

SCRIPT_PATH="/broad/IDP-Dx_work/nirmalya/research/BarcodeSplitter/"
ALLBC_FILE=$SCRIPT_PATH"AllSeq_96bc_v1.txt"
ALLBC_DICT=$SCRIPT_PATH"alldict.dat"
OUT_DIR="/broad/hptmp/nirmalya/"$EXPID"/"

mkdir $OUT_DIR

########  function for extracting command from command line
 
	extract_option ()
	{
		option_name=$1
		shift 
		default=$1
		shift 
		
		if [ $# == 0 ];then
			
			echo $default
		fi
		echo $@ | awk -v def=$default -v name=$option_name '{
			for(i=1; i < NF+1; i++)
			{
				if(i==1)			# intitalize to $default in case no option found
					print def
				if($i==name)
				{
					print $(i+1)
				}
			}
	
		}' | tail -1
	}


######################################
DICT_OW=`extract_option -DICT_OW N $@`


if [ ! -s $ALLBC_DICT ] || [ $DICT_OW == "Y" ];then
	echo $SCRIPT_PATH"dict_builder" -i $ALLBC_FILE -o $ALLBC_DICT
	$SCRIPT_PATH"dict_builder" -i $ALLBC_FILE -o $ALLBC_DICT
fi


echo $SCRIPT_PATH"bc_splitter" -t allseq -d $ALLBC_DICT --file1 $FASTQ1 --file2 $FASTQ2 -o $OUT_DIR
$SCRIPT_PATH"bc_splitter" -t allseq -d $ALLBC_DICT --file1 $FASTQ1 --file2 $FASTQ2 -o $OUT_DIR

cd $OUT_DIR 

Rscript $SCRIPT_PATH"plot_dist.R" $OUT_DIR"frequency_logfile.txt" $OUT_DIR"freq_plot.pdf" $EXPID 
